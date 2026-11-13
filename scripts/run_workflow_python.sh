#!/usr/bin/env bash
# Orchestrate the existing Python wrappers in the order of the top10 workflow.
# Completed outputs are reused; partial outputs are left untouched and cause a stop.

set -Eeuo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd -- "$SCRIPT_DIR/.." && pwd)"
MESS_ROOT="${MESS_ROOT:-$HOME/szakdolgozat/MeSS}"
MODE="dry-run"

usage() {
    cat <<USAGE
Usage: $(basename "$0") [--dry-run | --execute] [--data-root PATH]

Options:
  --dry-run       Print wrapper commands without executing them (default).
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

PYTHON="${PYTHON:-python3}"
command -v "$PYTHON" >/dev/null 2>&1 || { echo "ERROR: Python not found: $PYTHON" >&2; exit 1; }

READS="$MESS_ROOT/results/one_sample/fastq/top10.fq.gz"
ASSEMBLY="$MESS_ROOT/results/assembly/top10/flye/assembly.fasta"
MAPPING_DIR="$MESS_ROOT/results/assembly/top10/mapping"
BAM="$MAPPING_DIR/top10_vs_assembly.bam"
BAI="$BAM.bai"
FLAGSTAT="$MAPPING_DIR/flagstat.txt"
COVERAGE="$MAPPING_DIR/coverage.tsv"
BINS_DIR="$MESS_ROOT/results/binning/top10/semibin2/output_bins"
CHECKM2_REPORT="$MESS_ROOT/results/qc/top10/checkm2_lowmem_test2/quality_report.tsv"
FASTANI_DIR="$MESS_ROOT/results/qc/top10/fastani"
FASTANI_FORWARD="$FASTANI_DIR/bin_vs_reference.tsv"
FASTANI_REVERSE="$FASTANI_DIR/reference_vs_bins.tsv"
FASTANI_SUMMARY="$FASTANI_DIR/reference_recovery_summary.tsv"
QUAST_SUMMARY="$MESS_ROOT/results/qc/top10/quast/quast_summary.tsv"
INTEGRATED_SUMMARY="$REPO_ROOT/results/semibin2/10_reference_summary.tsv"

log() { printf '\n[%s] %s\n' "$(date '+%F %T')" "$*"; }
fail() { echo "ERROR: $*" >&2; exit 1; }
require_file() { [[ -s "$1" ]] || fail "Required file is missing or empty: $1"; }

run_wrapper() {
    local script="$1"
    shift
    local command=("$PYTHON" "$REPO_ROOT/scripts/$script" "$@")
    if [[ "$MODE" == "dry-run" ]]; then
        printf 'DRY RUN: '
        printf '%q ' "${command[@]}"
        printf '\n'
    else
        "${command[@]}"
    fi
}

# MeSS: current run_mess.py uses its configured/default MeSS root.
if [[ -s "$READS" ]]; then
    log "MeSS FASTQ already exists; reusing it: $READS"
else
    log "MeSS simulation"
    if [[ "$MODE" == "dry-run" ]]; then
        echo "DRY RUN: $PYTHON scripts/run_mess.py --execute"
    else
        # run_mess.py's default data root is ~/szakdolgozat/MeSS.
        # This is the configured project root used in the current setup.
        [[ "$MESS_ROOT" == "$HOME/szakdolgozat/MeSS" ]] || fail "run_mess.py uses its configured default MeSS root. For a custom --data-root, update run_mess.py's root setting before running this workflow."
        run_wrapper run_mess.py --execute
        require_file "$READS"
    fi
fi

# Flye.
if [[ -s "$ASSEMBLY" ]]; then
    log "Flye assembly already exists; skipping."
else
    if [[ "$MODE" == "execute" ]]; then require_file "$READS"; fi
    log "Flye assembly"
    run_wrapper run_flye.py --data-root "$MESS_ROOT" --execute
    [[ "$MODE" == "dry-run" ]] || require_file "$ASSEMBLY"
fi

# Mapping.
if [[ -s "$BAM" && -s "$BAI" && -s "$FLAGSTAT" && -s "$COVERAGE" ]]; then
    log "Mapping outputs already exist; skipping."
else
    if [[ "$MODE" == "execute" ]]; then
        require_file "$READS"
        require_file "$ASSEMBLY"
    fi
    log "Read mapping and samtools processing"
    run_wrapper run_mapping.py --data-root "$MESS_ROOT" --execute
    if [[ "$MODE" == "execute" ]]; then
        require_file "$BAM"; require_file "$BAI"; require_file "$FLAGSTAT"; require_file "$COVERAGE"
    fi
fi

# SemiBin2.
shopt -s nullglob
BIN_FILES=("$BINS_DIR"/SemiBin_*.fa.gz)
shopt -u nullglob
if ((${#BIN_FILES[@]} == 8)); then
    log "SemiBin2 bins already exist; skipping."
else
    ((${#BIN_FILES[@]} == 0)) || fail "SemiBin2 output has ${#BIN_FILES[@]} bin files; inspect the partial output before continuing."
    if [[ "$MODE" == "execute" ]]; then
        require_file "$ASSEMBLY"; require_file "$BAM"
    fi
    log "SemiBin2 binning"
    run_wrapper run_semibin2.py --data-root "$MESS_ROOT" --execute
    if [[ "$MODE" == "execute" ]]; then
        shopt -s nullglob; BIN_FILES=("$BINS_DIR"/SemiBin_*.fa.gz); shopt -u nullglob
        ((${#BIN_FILES[@]} > 0)) || fail "SemiBin2 did not produce bin files."
    fi
fi

# CheckM2.
if [[ -s "$CHECKM2_REPORT" ]]; then
    log "CheckM2 report already exists; skipping."
else
    if [[ "$MODE" == "execute" ]]; then
        compgen -G "$BINS_DIR/SemiBin_*.fa.gz" >/dev/null || fail "No SemiBin2 bin files found in $BINS_DIR"
    fi
    log "CheckM2 quality assessment"
    run_wrapper run_checkm2.py --data-root "$MESS_ROOT" --execute
    [[ "$MODE" == "dry-run" ]] || require_file "$CHECKM2_REPORT"
fi

# FastANI.
if [[ -s "$FASTANI_FORWARD" && -s "$FASTANI_REVERSE" && -s "$FASTANI_SUMMARY" ]]; then
    log "FastANI outputs already exist; skipping."
else
    log "FastANI identification"
    run_wrapper run_fastani.py --data-root "$MESS_ROOT" --execute
    if [[ "$MODE" == "execute" ]]; then
        require_file "$FASTANI_FORWARD"; require_file "$FASTANI_REVERSE"; require_file "$FASTANI_SUMMARY"
    fi
fi

# QUAST.
if [[ -s "$QUAST_SUMMARY" ]]; then
    log "QUAST summary already exists; skipping."
else
    log "QUAST evaluation"
    if [[ "$MODE" == "dry-run" && ! -s "$FASTANI_SUMMARY" ]]; then
        echo "DRY RUN: run_quast.py will be called after FastANI has produced its pairwise results."
    else
        if [[ "$MODE" == "execute" ]]; then require_file "$FASTANI_SUMMARY"; fi
        run_wrapper run_quast.py --data-root "$MESS_ROOT" --execute
        [[ "$MODE" == "dry-run" ]] || require_file "$QUAST_SUMMARY"
    fi
fi

# Final integrated summary.
if [[ -s "$INTEGRATED_SUMMARY" ]]; then
    log "Integrated summary already exists; skipping."
else
    log "Building integrated reference summary"
    if [[ "$MODE" == "dry-run" ]]; then
        printf 'DRY RUN: %q %q --data-root %q --output %q\n' \
            "$PYTHON" "$REPO_ROOT/scripts/build_10_reference_summary_v2.py" \
            "$MESS_ROOT" "$INTEGRATED_SUMMARY"
    else
        mkdir -p "$(dirname "$INTEGRATED_SUMMARY")"
        "$PYTHON" "$REPO_ROOT/scripts/build_10_reference_summary_v2.py" \
            --data-root "$MESS_ROOT" \
            --output "$INTEGRATED_SUMMARY"
        require_file "$INTEGRATED_SUMMARY"
    fi
fi

log "Workflow finished (mode: $MODE)."
if [[ "$MODE" == "dry-run" ]]; then
    echo "No wrappers were executed. Review the plan, then rerun with --execute."
fi
