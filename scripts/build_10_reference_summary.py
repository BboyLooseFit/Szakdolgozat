#!/usr/bin/env python3
"""Generate a summary table for the ten simulated reference genomes."""

import argparse
import csv
import re
from pathlib import Path


def read_tsv(path):
    """Read a tab-separated file and return its rows."""
    if not path.is_file():
        raise FileNotFoundError(f"Input file not found: {path}")

    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def get_accession(value):
    """Extract an NCBI GCF accession from a string."""
    match = re.search(r"GCF_\d+\.\d+", str(value))
    return match.group(0) if match else None


def get_bin_name(value):
    """Extract a bin name from a filename or path."""
    name = Path(str(value)).name
    return re.sub(r"\.fa(?:sta)?(?:\.gz)?$", "", name)


def find_name(row):
    """Find an organism name in metadata with varying column names."""
    candidates = (
        "name",
        "organism_name",
        "scientific_name",
        "species_name",
        "biosample.description.organism_name",
    )
    for column in candidates:
        if row.get(column):
            return row[column]
    return "NA"


def number(value, decimals=2):
    """Format a numeric value, preserving unavailable values as NA."""
    try:
        return f"{float(value):.{decimals}f}"
    except (TypeError, ValueError):
        return "NA"


def read_quast_report(path):
    """Read the key-value format used by QUAST report.tsv."""
    metrics = {}
    if not path.is_file():
        return metrics

    with path.open(newline="", encoding="utf-8") as handle:
        for row in csv.reader(handle, delimiter="\t"):
            if len(row) >= 2:
                metrics[row[0]] = row[1]
    return metrics


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Combine MeSS metadata, FastANI, CheckM2 and QUAST "
            "results into a ten-reference summary."
        )
    )
    parser.add_argument(
        "--mess-root",
        type=Path,
        required=True,
        help="Root of the MeSS working directory containing communities/ and results/.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("results/10_reference_summary.tsv"),
        help="Output TSV path (default: results/10_reference_summary.tsv).",
    )
    args = parser.parse_args()

    root = args.mess_root.expanduser().resolve()
    output = args.output

    community_file = root / "communities/top10_community.tsv"
    sample = root / "results/one_sample"
    qc = root / "results/qc"

    community_rows = read_tsv(community_file)
    taxonomy_rows = read_tsv(
        sample / "assembly_finder/taxonomy.tsv"
    )
    coverage_rows = read_tsv(sample / "coverages.tsv")
    assembly_rows = read_tsv(
        sample / "assembly_finder/assembly_summary.tsv"
    )
    read_rows = read_tsv(
        qc / "fastani/simulated_reads_by_reference.tsv"
    )
    checkm_rows = read_tsv(qc / "checkm2/quality_report.tsv")

    # Index taxonomy metadata by accession.
    taxonomy = {}
    for row in taxonomy_rows:
        accession = get_accession(
            row.get("accession", "")
        )
        if accession:
            taxonomy[accession] = row

    # Index MeSS coverage metadata by the reference accession.
    coverage = {}
    for row in coverage_rows:
        accession = get_accession(row.get("fasta", ""))
        if accession:
            coverage[accession] = row

    # Index reference assembly metadata.
    assemblies = {}
    for row in assembly_rows:
        accession = get_accession(
            row.get("accession", "")
            or row.get("assembly_accession", "")
            or row.get("assembly_name", "")
        )
        if accession:
            assemblies[accession] = row

    read_stats = {
        row["reference_accession"]: row
        for row in read_rows
    }

    # Index CheckM2 metrics by bin name.
    checkm = {}
    for row in checkm_rows:
        bin_name = get_bin_name(row.get("Name", ""))
        if bin_name:
            checkm[bin_name] = row

    # For each bin, identify its highest-ANI reference.
    forward_best = {}
    for row in read_tsv(qc / "fastani/bin_vs_reference.tsv"):
        bin_name = get_bin_name(row.get("query", ""))
        reference = get_accession(row.get("reference", ""))
        try:
            ani = float(row["ANI"])
        except (KeyError, TypeError, ValueError):
            # FastANI's standard output has no header; handled below if needed.
            continue
        if bin_name and reference:
            current = forward_best.get(bin_name)
            if current is None or ani > current["ani"]:
                forward_best[bin_name] = {
                    "reference": reference,
                    "ani": ani,
                }

    # Parse standard headerless FastANI output if necessary.
    forward_best = {}
    with (qc / "fastani/bin_vs_reference.tsv").open(
        newline="", encoding="utf-8"
    ) as handle:
        for row in csv.reader(handle, delimiter="\t"):
            if len(row) < 5:
                continue
            bin_name = get_bin_name(row[0])
            reference = get_accession(row[1])
            try:
                ani = float(row[2])
            except ValueError:
                continue
            if not reference:
                continue
            current = forward_best.get(bin_name)
            if current is None or ani > current["ani"]:
                forward_best[bin_name] = {
                    "reference": reference,
                    "ani": ani,
                }

    # For each reference, find its highest-ANI bin.
    reverse_best = {}
    with (qc / "fastani/reference_vs_bins.tsv").open(
        newline="", encoding="utf-8"
    ) as handle:
        for row in csv.reader(handle, delimiter="\t"):
            if len(row) < 5:
                continue
            reference = get_accession(row[0])
            bin_name = get_bin_name(row[1])
            try:
                ani = float(row[2])
                matched = int(row[3])
                total = int(row[4])
            except ValueError:
                continue
            if not reference:
                continue

            fraction = 100 * matched / total if total else 0
            current = reverse_best.get(reference)
            if current is None or ani > current["ani"]:
                reverse_best[reference] = {
                    "bin": bin_name,
                    "ani": ani,
                    "fraction": fraction,
                }

    columns = [
        "sample_id",
        "bacterium_name",
        "reference_accession",
        "relative_abundance_fraction",
        "relative_abundance_percent",
        "MeSS_reads_estimate",
        "simulated_reads_from_MeSS_BAM",
        "sample_total_reads",
        "MeSS_bases",
        "binning_tool",
        "best_candidate_bin",
        "FastANI_ANI_percent",
        "FastANI_reference_fragment_fraction_percent",
        "identifiable_0_1",
        "CheckM2_completeness_percent",
        "CheckM2_contamination_percent",
        "CheckM2_strain_heterogeneity",
        "QUAST_genome_fraction_percent",
        "QUAST_contig_count",
        "QUAST_total_assembly_length_bp",
        "QUAST_N50_bp",
        "QUAST_NG50_bp",
        "reference_genome_size_bp",
        "reference_GC_percent",
        "sequencing_technology",
        "simulation_model",
        "target_mean_read_length_bp",
        "simulation_accuracy",
        "simulation_ratio",
        "global_seed",
    ]

    output_rows = []
    for community in community_rows:
        accession = community["accession"]
        meta = assemblies.get(accession, {})
        tax = taxonomy.get(accession, {})
        cov = coverage.get(accession, {})
        reads = read_stats.get(accession, {})
        hit = reverse_best.get(accession, {})
        candidate = hit.get("bin", "NA")

        reciprocal = (
            candidate != "NA"
            and forward_best.get(candidate, {}).get("reference")
            == accession
        )
        identifiable = int(
            reciprocal
            and hit.get("ani", 0) >= 95
            and hit.get("fraction", 0) >= 50
        )

        # Only attach MAG metrics when the reference passes
        # the project's reciprocal-best-hit identification rule.
        cm = checkm.get(candidate, {}) if identifiable else {}
        quast = {}
        if identifiable:
            quast = read_quast_report(
                qc / "quast" / "reports" / candidate / "report.tsv"
            )

        abundance = float(community["seq_abundance"])
        reference_size = (
            meta.get("total_sequence_length")
            or cov.get("seq_len")
            or "NA"
        )

        output_rows.append({
            "sample_id": community.get("sample", "top10"),
            "bacterium_name": find_name(tax) if tax else find_name(meta),
            "reference_accession": accession,
            "relative_abundance_fraction": f"{abundance:.6f}",
            "relative_abundance_percent": f"{100 * abundance:.3f}",
            "MeSS_reads_estimate": cov.get("reads", "NA"),
            "simulated_reads_from_MeSS_BAM": reads.get(
                "primary_mapped_reads", "NA"
            ),
            "sample_total_reads": "67341",
            "MeSS_bases": cov.get("bases", "NA"),
            "binning_tool": "SemiBin2 2.5.0",
            "best_candidate_bin": candidate,
            "FastANI_ANI_percent": number(hit.get("ani")),
            "FastANI_reference_fragment_fraction_percent": number(
                hit.get("fraction")
            ),
            "identifiable_0_1": str(identifiable),
            "CheckM2_completeness_percent": number(
                cm.get("Completeness")
            ),
            "CheckM2_contamination_percent": number(
                cm.get("Contamination")
            ),
            "CheckM2_strain_heterogeneity": "NA",
            "QUAST_genome_fraction_percent": quast.get(
                "Genome fraction (%)", "NA"
            ),
            "QUAST_contig_count": quast.get("# contigs", "NA"),
            "QUAST_total_assembly_length_bp": quast.get(
                "Total length", "NA"
            ),
            "QUAST_N50_bp": quast.get("N50", "NA"),
            "QUAST_NG50_bp": quast.get("NG50", "NA"),
            "reference_genome_size_bp": reference_size,
            "reference_GC_percent": meta.get("gc_percent", "NA"),
            "sequencing_technology": "ONT",
            "simulation_model": "QSHMM-ONT-HQ",
            "target_mean_read_length_bp": "15000",
            "simulation_accuracy": "0.95",
            "simulation_ratio": "39:24:36",
            "global_seed": "42",
        })

    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=columns, delimiter="\t"
        )
        writer.writeheader()
        writer.writerows(output_rows)

    print(f"Output: {output}")
    print(f"Reference rows: {len(output_rows)}")
    print(
        "Identifiable references under project-specific rule: "
        f"{sum(int(row['identifiable_0_1']) for row in output_rows)}/10"
    )


if __name__ == "__main__":
    main()
