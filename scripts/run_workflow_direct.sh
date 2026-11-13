#!/usr/bin/env bash
# Direct command-line workflow for the top10 MeSS metagenomic simulation.
# Tools are invoked directly; Python is used only for TSV post-processing
# and for the final integrated summary.

set -Eeuo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd -- "$SCRIPT_DIR/.." && pwd)"
MESS_ROOT="${MESS_ROOT:-$HOME/szakdolgozat/MeSS}"
SAMPLE_ID="top10"
MODE="dry-run"

usage() {
    cat <<USAGE
Usage: $(basename "$0") [--dry-run | --execute] [--data-root PATH]

Options:
  --dry-run       Print actions without executing them (default).
  --execute       Execute missing workflow steps.
  --data-root     MeSS project/data directory (default: ~/szakdolgozat/MeSS).
  -h, --help      Show this help.

The workflow reuses completed outputs and refuses to overwrite partial outputs.
USAGE
}

while (($#)); do
    case "$1" in
        --dry-run) MODE="dry-run"; shift ;;
        --execute) MODE="execute"; shift ;;
        --data-root)
            [[ $# -ge 2 ]] || { echo "ERROR: --data-root needs a path." >&2; exit 2; }
            MESS_ROOT="$2"; shift 2 ;;
        -h|--help) usage; exit 0 ;;
        *) echo "ERROR: Unknown argument: $1" >&2; usage >&2; exit 2 ;;
    esac
done

MESS_ROOT="$(cd -- "$MESS_ROOT" 2>/dev/null && pwd)" || {
    echo "ERROR: MeSS data directory does not exist: $MESS_ROOT" >&2
    exit 1
}

# Input and output paths. These match the current top10 project configuration.
COMMUNITY_FILE="$MESS_ROOT/communities/top10_community.tsv"
SIM_DIR="$MESS_ROOT/results/one_sample"
READS="$SIM_DIR/fastq/top10.fq.gz"

FLYE_DIR="$MESS_ROOT/results/assembly/top10/flye"
ASSEMBLY="$FLYE_DIR/assembly.fasta"

MAPPING_DIR="$MESS_ROOT/results/assembly/top10/mapping"
BAM="$MAPPING_DIR/top10_vs_assembly.bam"
BAI="$BAM.bai"
FLAGSTAT="$MAPPING_DIR/flagstat.txt"
COVERAGE="$MAPPING_DIR/coverage.tsv"

BINS_DIR="$MESS_ROOT/results/binning/top10/semibin2/output_bins"
CHECKM2_DIR="$MESS_ROOT/results/qc/top10/checkm2_lowmem_test2"
CHECKM2_REPORT="$CHECKM2_DIR/quality_report.tsv"
CHECKM2_DB="$HOME/szakdolgozat/checkm2-db/CheckM2_database/uniref100.KO.1.dmnd"

REFERENCE_DIR="$SIM_DIR/assembly_finder/download"
FASTANI_DIR="$MESS_ROOT/results/qc/top10/fastani"
FASTANI_FORWARD="$FASTANI_DIR/bin_vs_reference.tsv"
FASTANI_REVERSE="$FASTANI_DIR/reference_vs_bins.tsv"
FASTANI_SUMMARY="$FASTANI_DIR/reference_recovery_summary.tsv"

QUAST_DIR="$MESS_ROOT/results/qc/top10/quast"
QUAST_REPORTS="$QUAST_DIR/reports"
QUAST_SUMMARY="$QUAST_DIR/quast_summary.tsv"

INTEGRATED_SUMMARY="$REPO_ROOT/results/semibin2/10_reference_summary.tsv"

# Tool parameters currently used in the tested top10 workflow.
MESS_THREADS=4
FLYE_THREADS=4
MAPPING_THREADS=4
SEMIBIN_ENVIRONMENT="human_gut"
SEMIBIN_SEQUENCING_TYPE="long_read"
CHECKM2_THREADS=1
FASTANI_THREADS=4
MINIMUM_ANI=95
MINIMUM_REFERENCE_FRACTION=50
QUAST_THREADS=2

log() { printf '\n[%s] %s\n' "$(date '+%F %T')" "$*"; }
fail() { echo "ERROR: $*" >&2; exit 1; }
require_file() { [[ -s "$1" ]] || fail "Required file is missing or empty: $1"; }
require_dir() { [[ -d "$1" ]] || fail "Required directory is missing: $1"; }

# Initialise micromamba activation for non-interactive bash scripts.
if [[ "$MODE" == "execute" ]]; then
    command -v micromamba >/dev/null 2>&1 || fail "micromamba was not found in PATH."
    eval "$(micromamba shell hook --shell bash)"
fi

activate_env() {
    local env_name="$1"
    log "Activating micromamba environment: $env_name"
    micromamba activate "$env_name"
}

deactivate_env() {
    micromamba deactivate
}

show_or_run() {
    if [[ "$MODE" == "dry-run" ]]; then
        printf 'DRY RUN: '
        printf '%q ' "$@"
        printf '\n'
    else
        "$@"
    fi
}

# Step 1: MeSS simulation. Reuse the existing FASTQ when it is already present.
if [[ -s "$READS" ]]; then
    log "MeSS FASTQ already exists; reusing it: $READS"
else
    [[ -s "$COMMUNITY_FILE" ]] || fail "Community table not found: $COMMUNITY_FILE"
    if [[ -d "$SIM_DIR" ]] && find "$SIM_DIR" -mindepth 1 -maxdepth 1 -print -quit | grep -q .; then
        fail "MeSS output directory exists but the expected FASTQ is missing. Inspect it before rerunning: $SIM_DIR"
    fi
    log "MeSS simulation"
    if [[ "$MODE" == "dry-run" ]]; then
        cat <<CMD
DRY RUN: micromamba activate mess
mess run \\
  -i "$COMMUNITY_FILE" --accession --tech ont --model QSHMM-ONT-HQ \\
  --mean-len 15000 --accuracy 0.95 --ratio 39:24:36 --seed 42 \\
  --threads $MESS_THREADS --bam -o "$SIM_DIR"
micromamba deactivate
CMD
    else
        activate_env mess
        mess run \
            -i "$COMMUNITY_FILE" \
            --accession \
            --tech ont \
            --model QSHMM-ONT-HQ \
            --mean-len 15000 \
            --accuracy 0.95 \
            --ratio 39:24:36 \
            --seed 42 \
            --threads "$MESS_THREADS" \
            --bam \
            -o "$SIM_DIR"
        deactivate_env
        require_file "$READS"
    fi
fi

# Step 2: Flye metagenomic assembly.
if [[ -s "$ASSEMBLY" ]]; then
    log "Flye assembly already exists; skipping: $ASSEMBLY"
else
    if [[ "$MODE" == "execute" ]]; then require_file "$READS"; fi
    if [[ -d "$FLYE_DIR" ]] && find "$FLYE_DIR" -mindepth 1 -maxdepth 1 -print -quit | grep -q .; then
        fail "Flye output directory is non-empty but assembly.fasta is missing: $FLYE_DIR"
    fi
    log "Flye assembly"
    if [[ "$MODE" == "dry-run" ]]; then
        cat <<CMD
DRY RUN: micromamba activate flye
flye --nano-hq "$READS" --meta --read-error 0.05 --threads $FLYE_THREADS --out-dir "$FLYE_DIR"
micromamba deactivate
CMD
    else
        mkdir -p "$(dirname "$FLYE_DIR")"
        activate_env flye
        flye --nano-hq "$READS" --meta --read-error 0.05 --threads "$FLYE_THREADS" --out-dir "$FLYE_DIR"
        deactivate_env
        require_file "$ASSEMBLY"
    fi
fi

# Step 3: Map reads to assembly, sort/index BAM, and create mapping statistics.
if [[ -s "$BAM" && -s "$BAI" && -s "$FLAGSTAT" && -s "$COVERAGE" ]]; then
    log "Mapping BAM and statistics already exist; skipping."
else
    if [[ "$MODE" == "execute" ]]; then
        require_file "$READS"
        require_file "$ASSEMBLY"
    fi
    for f in "$BAM" "$BAI" "$FLAGSTAT" "$COVERAGE"; do
        [[ ! -e "$f" ]] || fail "Partial mapping output exists; refusing to overwrite: $f"
    done
    log "minimap2 mapping and samtools processing"
    if [[ "$MODE" == "dry-run" ]]; then
        cat <<CMD
DRY RUN: micromamba activate flye
minimap2 -ax map-ont "$ASSEMBLY" "$READS" | samtools sort -@ $MAPPING_THREADS -o "$BAM"
samtools index "$BAM"
samtools flagstat "$BAM" > "$FLAGSTAT"
samtools coverage "$BAM" > "$COVERAGE"
micromamba deactivate
CMD
    else
        mkdir -p "$MAPPING_DIR"
        activate_env flye
        minimap2 -ax map-ont "$ASSEMBLY" "$READS" | samtools sort -@ "$MAPPING_THREADS" -o "$BAM"
        samtools index "$BAM"
        samtools flagstat "$BAM" > "$FLAGSTAT"
        samtools coverage "$BAM" > "$COVERAGE"
        deactivate_env
        require_file "$BAM"
        require_file "$BAI"
        require_file "$FLAGSTAT"
        require_file "$COVERAGE"
    fi
fi

# Step 4: SemiBin2 single-sample long-read binning.
shopt -s nullglob
BIN_FILES=("$BINS_DIR"/SemiBin_*.fa.gz)
shopt -u nullglob
if ((${#BIN_FILES[@]} == 8)); then
    log "Found 8 SemiBin2 bins; skipping binning."
elif ((${#BIN_FILES[@]} > 0)); then
    fail "SemiBin2 output contains ${#BIN_FILES[@]} bin files, not the expected 8. Inspect before rerunning: $BINS_DIR"
else
    if [[ "$MODE" == "execute" ]]; then
        require_file "$ASSEMBLY"
        require_file "$BAM"
    fi
    [[ ! -e "$BINS_DIR" ]] || {
        [[ -d "$BINS_DIR" ]] && [[ -z "$(find "$BINS_DIR" -mindepth 1 -maxdepth 1 -print -quit)" ]] || fail "SemiBin2 output directory already contains data: $BINS_DIR"
    }
    log "SemiBin2 binning"
    if [[ "$MODE" == "dry-run" ]]; then
        cat <<CMD
DRY RUN: micromamba activate semibin
SemiBin2 single_easy_bin --environment "$SEMIBIN_ENVIRONMENT" --sequencing-type "$SEMIBIN_SEQUENCING_TYPE" -i "$ASSEMBLY" -b "$BAM" -o "$BINS_DIR"
micromamba deactivate
CMD
    else
        mkdir -p "$(dirname "$BINS_DIR")"
        activate_env semibin
        SemiBin2 single_easy_bin \
            --environment "$SEMIBIN_ENVIRONMENT" \
            --sequencing-type "$SEMIBIN_SEQUENCING_TYPE" \
            -i "$ASSEMBLY" \
            -b "$BAM" \
            -o "$BINS_DIR"
        deactivate_env
        shopt -s nullglob
        BIN_FILES=("$BINS_DIR"/SemiBin_*.fa.gz)
        shopt -u nullglob
        ((${#BIN_FILES[@]} > 0)) || fail "SemiBin2 completed but no SemiBin_*.fa.gz files were found."
    fi
fi

# Step 5: CheckM2 quality assessment in low-memory mode.
if [[ -s "$CHECKM2_REPORT" ]]; then
    log "CheckM2 quality_report.tsv already exists; skipping."
else
    if [[ "$MODE" == "execute" ]]; then
        require_dir "$BINS_DIR"
        compgen -G "$BINS_DIR/*.fa.gz" >/dev/null || fail "No .fa.gz bins found in $BINS_DIR"
        require_file "$CHECKM2_DB"
    fi
    if [[ -d "$CHECKM2_DIR" ]] && find "$CHECKM2_DIR" -mindepth 1 -maxdepth 1 -print -quit | grep -q .; then
        fail "CheckM2 output directory is non-empty but quality_report.tsv is missing: $CHECKM2_DIR"
    fi
    log "CheckM2 quality assessment"
    if [[ "$MODE" == "dry-run" ]]; then
        cat <<CMD
DRY RUN: micromamba activate checkm2
checkm2 predict --input "$BINS_DIR" --output-directory "$CHECKM2_DIR" --threads $CHECKM2_THREADS --database_path "$CHECKM2_DB" --extension gz --lowmem
micromamba deactivate
CMD
    else
        mkdir -p "$(dirname "$CHECKM2_DIR")"
        activate_env checkm2
        checkm2 predict \
            --input "$BINS_DIR" \
            --output-directory "$CHECKM2_DIR" \
            --threads "$CHECKM2_THREADS" \
            --database_path "$CHECKM2_DB" \
            --extension gz \
            --lowmem
        deactivate_env
        require_file "$CHECKM2_REPORT"
    fi
fi

# Step 6: Run FastANI in both directions and create the reference recovery table.
if [[ -s "$FASTANI_FORWARD" && -s "$FASTANI_REVERSE" && -s "$FASTANI_SUMMARY" ]]; then
    log "FastANI outputs and recovery summary already exist; skipping."
else
    for f in "$FASTANI_FORWARD" "$FASTANI_REVERSE" "$FASTANI_SUMMARY"; do
        [[ ! -e "$f" ]] || fail "Partial FastANI output exists; refusing to overwrite: $f"
    done
    if [[ "$MODE" == "execute" ]]; then
        require_dir "$BINS_DIR"
        require_dir "$REFERENCE_DIR"
    fi
    BIN_PATHS=()
    REF_PATHS=()
    if [[ -d "$BINS_DIR" ]]; then
        mapfile -t BIN_PATHS < <(find "$BINS_DIR" -maxdepth 1 -type f -name 'SemiBin_*.fa.gz' -print | sort)
    fi
    if [[ -d "$REFERENCE_DIR" ]]; then
        mapfile -t REF_PATHS < <(find "$REFERENCE_DIR" -type f -name '*.fna.gz' -print | sort)
    fi
    if [[ "$MODE" == "execute" ]]; then
        ((${#BIN_PATHS[@]} > 0)) || fail "No SemiBin_*.fa.gz bin files found."
        ((${#REF_PATHS[@]} > 0)) || fail "No .fna.gz reference genomes found under $REFERENCE_DIR"
    fi
    log "FastANI: bins against references and references against bins"
    if [[ "$MODE" == "dry-run" ]]; then
        printf 'DRY RUN: expected inputs: %s bin file(s), %s reference genome file(s).\n' \
            "${#BIN_PATHS[@]}" "${#REF_PATHS[@]}"
        cat <<CMD
DRY RUN: micromamba activate fastani
fastANI --ql "$FASTANI_DIR/query_bins.txt" --rl "$FASTANI_DIR/reference_genomes.txt" --threads $FASTANI_THREADS -o "$FASTANI_FORWARD"
fastANI --ql "$FASTANI_DIR/reference_genomes.txt" --rl "$FASTANI_DIR/query_bins.txt" --threads $FASTANI_THREADS -o "$FASTANI_REVERSE"
# A Python post-processing block will create reference_recovery_summary.tsv using reciprocal-best-hit,
# ANI >= $MINIMUM_ANI and reference fragment fraction >= $MINIMUM_REFERENCE_FRACTION%.
micromamba deactivate
CMD
    else
        if [[ -d "$FASTANI_DIR" ]] && find "$FASTANI_DIR" -mindepth 1 -maxdepth 1 -print -quit | grep -q .; then
            fail "FastANI output directory is non-empty but one or more expected outputs are missing: $FASTANI_DIR"
        fi
        mkdir -p "$FASTANI_DIR"
        printf '%s\n' "${BIN_PATHS[@]}" > "$FASTANI_DIR/query_bins.txt"
        printf '%s\n' "${REF_PATHS[@]}" > "$FASTANI_DIR/reference_genomes.txt"
        activate_env fastani
        fastANI --ql "$FASTANI_DIR/query_bins.txt" --rl "$FASTANI_DIR/reference_genomes.txt" --threads "$FASTANI_THREADS" -o "$FASTANI_FORWARD"
        fastANI --ql "$FASTANI_DIR/reference_genomes.txt" --rl "$FASTANI_DIR/query_bins.txt" --threads "$FASTANI_THREADS" -o "$FASTANI_REVERSE"
        deactivate_env
        require_file "$FASTANI_FORWARD"
        require_file "$FASTANI_REVERSE"
        python3 - "$FASTANI_DIR" "$MINIMUM_ANI" "$MINIMUM_REFERENCE_FRACTION" <<'PY'
import csv
import re
import sys
from pathlib import Path

outdir = Path(sys.argv[1])
minimum_ani = float(sys.argv[2])
minimum_fraction = float(sys.argv[3])
forward_file = outdir / "bin_vs_reference.tsv"
reverse_file = outdir / "reference_vs_bins.tsv"
refs_file = outdir / "reference_genomes.txt"


def read_hits(path):
    hits = []
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            parts = line.rstrip("\r\n").split("\t")
            if len(parts) < 5:
                continue
            query, reference, ani, matched, total = parts[:5]
            fraction = float(parts[5]) if len(parts) > 5 else (float(matched) / float(total) if float(total) else 0.0)
            if fraction <= 1.0:
                fraction *= 100.0
            hits.append({
                "query": query, "reference": reference, "ani": float(ani),
                "matched": int(matched), "total": int(total), "fraction": fraction,
            })
    return hits


def best_by(hits, field):
    result = {}
    for hit in hits:
        key = hit[field]
        score = (hit["ani"], hit["fraction"], hit["matched"])
        if key not in result or score > result[key][0]:
            result[key] = (score, hit)
    return {key: value[1] for key, value in result.items()}


def accession(path):
    match = re.search(r"GCF_\d+\.\d+", Path(path).name)
    return match.group(0) if match else Path(path).name

forward = read_hits(forward_file)  # query=bin, reference=genome
reverse = read_hits(reverse_file)  # query=genome, reference=bin
forward_best = best_by(forward, "query")
reverse_best = best_by(reverse, "query")
with refs_file.open(encoding="utf-8") as handle:
    references = [line.strip() for line in handle if line.strip()]

fields = [
    "Reference", "Best_bin", "ANI_percent", "Matched_fragments",
    "Total_fragments", "Reference_fragment_fraction_percent",
    "Reciprocal_best_hit", "Identifiable_0_1",
]
out = outdir / "reference_recovery_summary.tsv"
with out.open("w", newline="", encoding="utf-8") as handle:
    writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    for ref in references:
        hit = reverse_best.get(ref)
        if hit is None:
            writer.writerow({
                "Reference": accession(ref), "Best_bin": "NA", "ANI_percent": "NA",
                "Matched_fragments": "NA", "Total_fragments": "NA",
                "Reference_fragment_fraction_percent": "NA",
                "Reciprocal_best_hit": 0, "Identifiable_0_1": 0,
            })
            continue
        bin_path = hit["reference"]
        bin_name = Path(bin_path).name
        if bin_name.endswith(".fa.gz"):
            bin_name = bin_name[:-6]
        reciprocal_hit = forward_best.get(bin_path)
        reciprocal = int(reciprocal_hit is not None and reciprocal_hit["reference"] == ref)
        identifiable = int(reciprocal and hit["ani"] >= minimum_ani and hit["fraction"] >= minimum_fraction)
        writer.writerow({
            "Reference": accession(ref), "Best_bin": bin_name,
            "ANI_percent": f"{hit['ani']:.4f}",
            "Matched_fragments": hit["matched"], "Total_fragments": hit["total"],
            "Reference_fragment_fraction_percent": f"{hit['fraction']:.2f}",
            "Reciprocal_best_hit": reciprocal, "Identifiable_0_1": identifiable,
        })
print(f"Created {out}")
PY
        require_file "$FASTANI_SUMMARY"
    fi
fi

# Step 7: QUAST for reference/bin pairs passing the FastANI identification criteria.
if [[ -s "$QUAST_SUMMARY" ]]; then
    log "QUAST summary already exists; skipping."
else
    if [[ "$MODE" == "dry-run" ]]; then
        log "QUAST evaluation (runs only for FastANI-identifiable reference/bin pairs)"
        if [[ -s "$FASTANI_SUMMARY" && -s "$FASTANI_DIR/query_bins.txt" && -s "$FASTANI_DIR/reference_genomes.txt" ]]; then
            PAIR_TMP="$(mktemp)"
            python3 - "$FASTANI_SUMMARY" "$FASTANI_DIR/query_bins.txt" "$FASTANI_DIR/reference_genomes.txt" <<'PY' > "$PAIR_TMP"
import csv
import re
import sys
from pathlib import Path
summary_file, bins_file, refs_file = map(Path, sys.argv[1:])
bins = [line.strip() for line in bins_file.open(encoding="utf-8") if line.strip()]
refs = [line.strip() for line in refs_file.open(encoding="utf-8") if line.strip()]
def bin_id(path):
    name = Path(path).name
    return name[:-6] if name.endswith(".fa.gz") else Path(name).stem
def ref_id(path):
    match = re.search(r"GCF_\d+\.\d+", Path(path).name)
    return match.group(0) if match else Path(path).name
bin_map = {bin_id(path): path for path in bins}
ref_map = {ref_id(path): path for path in refs}
print("reference\tbin\taccession\tbin_name\tani\tref_fraction\treciprocal\tidentifiable")
with summary_file.open(newline="", encoding="utf-8") as handle:
    for row in csv.DictReader(handle, delimiter="\t"):
        if str(row.get("Identifiable_0_1", "0")).strip() != "1":
            continue
        accession, bin_name = row["Reference"], row["Best_bin"]
        if accession not in ref_map or bin_name not in bin_map:
            raise SystemExit(f"Could not resolve reference/bin pair: {accession}, {bin_name}")
        print("\t".join([ref_map[accession], bin_map[bin_name], accession, bin_name,
            row.get("ANI_percent", "NA"), row.get("Reference_fragment_fraction_percent", "NA"),
            row.get("Reciprocal_best_hit", "0"), row.get("Identifiable_0_1", "0")]))
PY
            while IFS=$'\t' read -r ref_file bin_file accession bin_name ani frac reciprocal identifiable; do
                [[ "$ref_file" == "reference" ]] && continue
                printf 'quast.py -r %q %q --threads %s -o %q\n' "$ref_file" "$bin_file" "$QUAST_THREADS" "$QUAST_REPORTS/${bin_name}_vs_${accession}"
            done < "$PAIR_TMP"
            rm -f "$PAIR_TMP"
        else
            echo "DRY RUN: FastANI output is not available yet; QUAST input pairs will be derived after FastANI completes."
            echo "quast.py -r <identified_reference.fna.gz> <identified_bin.fa.gz> --threads $QUAST_THREADS -o <reports/bin_vs_accession>"
        fi
        echo "DRY RUN: the workflow will then create quast_summary.tsv and the integrated summary."
    else
        require_file "$FASTANI_SUMMARY"
        require_file "$FASTANI_DIR/query_bins.txt"
        require_file "$FASTANI_DIR/reference_genomes.txt"
        PAIRS_FILE="$QUAST_DIR/identified_pairs.tsv"
        PAIR_TMP="$(mktemp)"
        python3 - "$FASTANI_SUMMARY" "$FASTANI_DIR/query_bins.txt" "$FASTANI_DIR/reference_genomes.txt" <<'PY' > "$PAIR_TMP"
import csv
import re
import sys
from pathlib import Path
summary_file, bins_file, refs_file = map(Path, sys.argv[1:])
bins = [line.strip() for line in bins_file.open(encoding="utf-8") if line.strip()]
refs = [line.strip() for line in refs_file.open(encoding="utf-8") if line.strip()]
def bin_id(path):
    name = Path(path).name
    return name[:-6] if name.endswith(".fa.gz") else Path(name).stem
def ref_id(path):
    match = re.search(r"GCF_\d+\.\d+", Path(path).name)
    return match.group(0) if match else Path(path).name
bin_map = {bin_id(path): path for path in bins}
ref_map = {ref_id(path): path for path in refs}
print("reference\tbin\taccession\tbin_name\tani\tref_fraction\treciprocal\tidentifiable")
with summary_file.open(newline="", encoding="utf-8") as handle:
    for row in csv.DictReader(handle, delimiter="\t"):
        if str(row.get("Identifiable_0_1", "0")).strip() != "1":
            continue
        accession, bin_name = row["Reference"], row["Best_bin"]
        if accession not in ref_map or bin_name not in bin_map:
            raise SystemExit(f"Could not resolve reference/bin pair: {accession}, {bin_name}")
        print("\t".join([ref_map[accession], bin_map[bin_name], accession, bin_name,
            row.get("ANI_percent", "NA"), row.get("Reference_fragment_fraction_percent", "NA"),
            row.get("Reciprocal_best_hit", "0"), row.get("Identifiable_0_1", "0")]))
PY
        mkdir -p "$QUAST_REPORTS"
        cp "$PAIR_TMP" "$PAIRS_FILE"
        rm -f "$PAIR_TMP"
        activate_env quast311
        while IFS=$'\t' read -r ref_file bin_file accession bin_name ani frac reciprocal identifiable; do
            [[ "$ref_file" == "reference" ]] && continue
            pair_out="$QUAST_REPORTS/${bin_name}_vs_${accession}"
            if [[ -s "$pair_out/report.tsv" ]]; then
                log "QUAST report already exists; skipping $bin_name vs $accession"
            elif [[ -d "$pair_out" ]] && find "$pair_out" -mindepth 1 -maxdepth 1 -print -quit | grep -q .; then
                fail "Partial QUAST output exists; refusing to overwrite: $pair_out"
            else
                quast.py -r "$ref_file" "$bin_file" --threads "$QUAST_THREADS" -o "$pair_out"
                require_file "$pair_out/report.tsv"
            fi
        done < "$PAIRS_FILE"
        deactivate_env

        python3 - "$FASTANI_SUMMARY" "$QUAST_REPORTS" "$QUAST_SUMMARY" <<'PY'
import csv
import sys
from pathlib import Path
fastani_summary, reports_dir, out_file = map(Path, sys.argv[1:])
def report_values(report_file):
    values = {}
    with report_file.open(newline="", encoding="utf-8") as handle:
        for row in csv.reader(handle, delimiter="\t"):
            if len(row) < 2:
                continue
            key, value = row[0].strip(), row[1].strip()
            low = key.lower()
            if low.startswith("genome fraction"):
                values["Genome_fraction_percent"] = value
            elif low == "total length":
                values["Total_length"] = value
            elif low.startswith("# contigs") and "(" not in low:
                values["Contig_count"] = value
            elif low == "n50":
                values["N50"] = value
            elif low == "ng50":
                values["NG50"] = value
            elif low.startswith("# misassemblies") or low == "misassemblies":
                values["Misassemblies"] = value
    return values
fields = ["Reference", "Bin", "ANI", "Reference_fragment_fraction_percent", "Reciprocal_best_hit",
    "Identifiable_0_1", "Genome_fraction_percent", "Total_length", "Contig_count", "N50", "NG50", "Misassemblies"]
with fastani_summary.open(newline="", encoding="utf-8") as source, out_file.open("w", newline="", encoding="utf-8") as target:
    reader = csv.DictReader(source, delimiter="\t")
    writer = csv.DictWriter(target, fieldnames=fields, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    for row in reader:
        accession = row["Reference"]
        bin_name = row.get("Best_bin", "NA")
        metrics = {}
        if row.get("Identifiable_0_1", "0") == "1":
            report_file = reports_dir / f"{bin_name}_vs_{accession}" / "report.tsv"
            if report_file.is_file():
                metrics = report_values(report_file)
        writer.writerow({
            "Reference": accession, "Bin": bin_name,
            "ANI": row.get("ANI_percent", "NA"),
            "Reference_fragment_fraction_percent": row.get("Reference_fragment_fraction_percent", "NA"),
            "Reciprocal_best_hit": row.get("Reciprocal_best_hit", "0"),
            "Identifiable_0_1": row.get("Identifiable_0_1", "0"),
            **{field: metrics.get(field, "NA") for field in fields[6:]},
        })
print(f"Created {out_file}")
PY
        require_file "$QUAST_SUMMARY"
    fi
fi

# Step 8: Build the integrated 10-reference summary table.
if [[ -s "$INTEGRATED_SUMMARY" ]]; then
    log "Integrated summary already exists; skipping: $INTEGRATED_SUMMARY"
else
    log "Building integrated reference summary"
    if [[ "$MODE" == "dry-run" ]]; then
        printf 'DRY RUN: python3 %q --data-root %q --output %q\n' \
            "$REPO_ROOT/scripts/build_10_reference_summary_v2.py" \
            "$MESS_ROOT" "$INTEGRATED_SUMMARY"
    else
        mkdir -p "$(dirname "$INTEGRATED_SUMMARY")"
        python3 "$REPO_ROOT/scripts/build_10_reference_summary_v2.py" \
            --data-root "$MESS_ROOT" \
            --output "$INTEGRATED_SUMMARY"
        require_file "$INTEGRATED_SUMMARY"
    fi
fi

log "Workflow finished (mode: $MODE)."
if [[ "$MODE" == "dry-run" ]]; then
    echo "No tools were executed. Review the planned steps, then rerun with --execute."
fi
