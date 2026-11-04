#!/usr/bin/env python3

"""Build a ten-reference summary from the metagenomic workflow outputs."""

import argparse
import csv
import gzip
import re
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]


def read_tsv(path):
    if not path.is_file():
        return []

    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def find_sample_fastq(sample_dir, sample_id):
    """Find the FASTQ belonging to the selected sample."""
    fastq_dir = sample_dir / "fastq"

    candidates = [
        fastq_dir / f"{sample_id}.fq.gz",
        fastq_dir / f"{sample_id}.fastq.gz",
        fastq_dir / f"{sample_id}.fq",
        fastq_dir / f"{sample_id}.fastq",
    ]

    for candidate in candidates:
        if candidate.is_file():
            return candidate

    raise FileNotFoundError(
        f"FASTQ file for sample {sample_id} not found in {fastq_dir}"
    )


def fastq_statistics(path):
    """Calculate descriptive read statistics from a FASTQ file."""
    opener = gzip.open if str(path).endswith(".gz") else open

    count = 0
    total_bases = 0
    min_length = None
    max_length = 0

    with opener(path, "rt", encoding="utf-8") as handle:
        while True:
            header = handle.readline()

            if not header:
                break

            sequence = handle.readline().rstrip("\r\n")
            separator = handle.readline()
            quality = handle.readline()

            if not sequence or not separator or not quality:
                raise ValueError(f"Incomplete FASTQ record in {path}")

            if not header.startswith("@") or not separator.startswith("+"):
                raise ValueError(f"Invalid FASTQ format in {path}")

            length = len(sequence)
            count += 1
            total_bases += length

            min_length = (
                length if min_length is None
                else min(min_length, length)
            )
            max_length = max(max_length, length)

    return {
        "read_count": count,
        "total_bases": total_bases,
        "mean_length": total_bases / count if count else 0,
        "min_length": min_length if count else "NA",
        "max_length": max_length if count else "NA",
    }


def read_parameters(config_file, sample_id=None):
    rows = read_tsv(config_file)

    if sample_id:
        rows = [r for r in rows if r.get("sample_id") == sample_id]

    if not rows:
        raise ValueError(f"No matching configuration in {config_file}")

    return rows[0]


def resolve_path(value, data_root):
    path = Path(value).expanduser()

    if not path.is_absolute():
        path = data_root / path

    return path.resolve()


def get_accession(value):
    match = re.search(r"GCF_\d+\.\d+", str(value))
    return match.group(0) if match else None


def get_bin_name(value):
    name = Path(str(value)).name

    for suffix in (
        ".fasta.gz", ".fna.gz", ".fa.gz",
        ".fasta", ".fna", ".fa",
    ):
        if name.endswith(suffix):
            return name[:-len(suffix)]

    return name


def number(value, decimals=2):
    try:
        return f"{float(value):.{decimals}f}"
    except (TypeError, ValueError):
        return "NA"


def index_by_accession(rows, fields):
    result = {}

    for row in rows:
        accession = None

        for field in fields:
            accession = get_accession(row.get(field, ""))
            if accession:
                break

        if accession:
            result[accession] = row

    return result


def parse_forward_fastani(path):
    """Find the best reference for each bin."""
    best = {}

    for row in read_tsv(path):
        query = row.get("query", "")
        reference = row.get("reference", "")

        # FastANI output is normally headerless.
        if not query:
            continue

        try:
            ani = float(row.get("ANI", ""))
        except (TypeError, ValueError):
            continue

        bin_name = get_bin_name(query)
        accession = get_accession(reference)

        if not bin_name or not accession:
            continue

        current = best.get(bin_name)

        if current is None or ani > current["ani"]:
            best[bin_name] = {
                "reference": accession,
                "ani": ani,
            }

    # Read the original headerless FastANI format.
    if not best and path.is_file():
        with path.open(encoding="utf-8") as handle:
            for line in handle:
                fields = line.rstrip("\n").split("\t")

                if len(fields) < 5:
                    continue

                try:
                    ani = float(fields[2])
                except ValueError:
                    continue

                bin_name = get_bin_name(fields[0])
                accession = get_accession(fields[1])

                if not bin_name or not accession:
                    continue

                current = best.get(bin_name)

                if current is None or ani > current["ani"]:
                    best[bin_name] = {
                        "reference": accession,
                        "ani": ani,
                    }

    return best


def parse_reverse_fastani(path):
    """Find the best bin for each reference."""
    best = {}

    if not path.is_file():
        return best

    with path.open(encoding="utf-8") as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")

            if len(fields) < 5:
                continue

            accession = get_accession(fields[0])
            bin_name = get_bin_name(fields[1])

            try:
                ani = float(fields[2])
                matched = int(fields[3])
                total = int(fields[4])
            except ValueError:
                continue

            if not accession or total == 0:
                continue

            fraction = 100 * matched / total
            current = best.get(accession)

            if current is None or ani > current["ani"]:
                best[accession] = {
                    "bin": bin_name,
                    "ani": ani,
                    "matched": matched,
                    "total": total,
                    "fraction": fraction,
                }

    return best


def read_quast_report(path):
    metrics = {}

    if not path.is_file():
        return metrics

    with path.open(newline="", encoding="utf-8") as handle:
        for row in csv.reader(handle, delimiter="\t"):
            if len(row) >= 2:
                metrics[row[0].strip()] = row[1].strip()

    return metrics


def find_quast_report(new_root, data_root, bin_name, accession):
    candidates = [
        new_root / "reports" / f"{bin_name}_vs_{accession}" / "report.tsv",
        data_root / "results/qc/quast/reports" /
        f"{bin_name}_vs_{accession}" / "report.tsv",
        data_root / "results/qc/quast/reports" /
        bin_name / "report.tsv",
    ]

    for path in candidates:
        if path.is_file():
            return path

    return None


def main():
    parser = argparse.ArgumentParser(
        description="Build a ten-reference metagenomic quality summary."
    )

    parser.add_argument(
        "--data-root",
        type=Path,
        default=Path.home() / "szakdolgozat" / "MeSS",
        help="Root directory of the MeSS project.",
    )

    parser.add_argument(
        "--sample-id",
        default="top10",
        help="Sample identifier.",
    )

    parser.add_argument(
        "--output",
        type=Path,
        default=REPO_ROOT / "results/pilot/10_reference_summary_test.tsv",
        help="Output TSV file.",
    )

    args = parser.parse_args()

    data_root = args.data_root.expanduser().resolve()
    if not data_root.is_dir():
        raise NotADirectoryError(data_root)

    config_dir = REPO_ROOT / "config"

    # Load the configuration files.
    simulation = read_parameters(
        config_dir / "simulation_parameters.tsv",
        args.sample_id,
    )

    fastani_config = read_parameters(
        config_dir / "fastani_parameters.tsv",
        args.sample_id,
    )

    checkm_config = read_parameters(
        config_dir / "checkm2_parameters.tsv",
        args.sample_id,
    )

    quast_config = read_parameters(
        config_dir / "quast_parameters.tsv",
        args.sample_id,
    )

    semibin_config = read_parameters(
        config_dir / "semibin2_parameters.tsv",
        args.sample_id,
    )

    # MeSS input and metadata.
    community_file = resolve_path(
        simulation["community_file"], data_root
    )

    sample_dir = resolve_path(
        simulation["output_dir"], data_root
    )

    community_rows = read_tsv(community_file)
    coverage_rows = read_tsv(sample_dir / "coverages.tsv")

    taxonomy_rows = read_tsv(
        sample_dir / "assembly_finder/taxonomy.tsv"
    )

    assembly_rows = read_tsv(
        sample_dir / "assembly_finder/assembly_summary.tsv"
    )

    if len(community_rows) != 10:
        print(
            f"WARNING: expected 10 community rows, "
            f"found {len(community_rows)}."
        )

    taxonomy = index_by_accession(
        taxonomy_rows, ["accession"]
    )

    assemblies = index_by_accession(
        assembly_rows, ["accession", "assembly_accession",
                        "assembly_name"]
    )

    coverage = index_by_accession(
        coverage_rows, ["fasta", "accession"]
    )

    # MeSS read statistics: use the existing pilot file if available.
    read_statistics_candidates = [
        data_root / "results/qc/fastani/simulated_reads_by_reference.tsv",
        resolve_path(
            fastani_config["output_dir"], data_root
        ) / "simulated_reads_by_reference.tsv",
    ]

    read_statistics_file = next(
        (p for p in read_statistics_candidates if p.is_file()),
        None,
    )

    read_statistics = (
        index_by_accession(
            read_tsv(read_statistics_file),
            ["reference_accession"],
        )
        if read_statistics_file else {}
    )

    # CheckM2: prefer the new output directory, then the pilot output.
    checkm_dir = resolve_path(
        checkm_config["output_dir"], data_root
    )

    checkm_candidates = [
        checkm_dir / "quality_report.tsv",
        data_root / "results/qc/checkm2/quality_report.tsv",
    ]

    checkm_file = next(
        (p for p in checkm_candidates if p.is_file()),
        None,
    )

    checkm = {}

    if checkm_file:
        for row in read_tsv(checkm_file):
            bin_name = get_bin_name(row.get("Name", ""))
            if bin_name:
                checkm[bin_name] = row

    # FastANI: new output directory first, old pilot directory second.
    fastani_dir = resolve_path(
        fastani_config["output_dir"], data_root
    )

    fastani_directories = [
        fastani_dir,
        data_root / "results/qc/fastani",
    ]

    forward_file = None
    reverse_file = None

    for directory in fastani_directories:
        forward_candidate = directory / "bin_vs_reference.tsv"
        reverse_candidate = directory / "reference_vs_bins.tsv"

        if forward_candidate.is_file() and reverse_candidate.is_file():
            forward_file = forward_candidate
            reverse_file = reverse_candidate
            break

    if not forward_file or not reverse_file:
        raise FileNotFoundError(
            "The two FastANI comparison files could not be found."
        )

    forward_best = parse_forward_fastani(forward_file)
    reverse_best = parse_reverse_fastani(reverse_file)

    minimum_ani = float(fastani_config["minimum_ani"])
    minimum_fraction = float(
        fastani_config["minimum_reference_fragment_fraction"]
    )

    # QUAST output directory.
    quast_dir = resolve_path(
        quast_config["output_dir"], data_root
    )

    # Count the actual number of reads in the simulated FASTQ.
    # MeSS per-reference read estimates are kept separate.
    fastq_file = find_sample_fastq(sample_dir, args.sample_id)
    fastq_stats = fastq_statistics(fastq_file)
    sample_total_reads = str(fastq_stats["read_count"])

    # Independent cross-check using primary read counts from the BAM summary.
    bam_primary_total = 0

    for row in read_statistics.values():
        try:
            bam_primary_total += int(float(row.get("primary_mapped_reads", 0)))
        except (TypeError, ValueError):
            pass

    # Binning tool information.
    binning_tool = "SemiBin2 2.5.0"

    output_columns = [
        "sample_id",
        "bacterium_name",
        "reference_accession",
        "relative_abundance_fraction",
        "relative_abundance_percent",
        "MeSS_reads_estimate",
        "simulated_reads_from_MeSS_BAM",
        "sample_total_reads",
        "sample_total_bases_bp",
        "sample_mean_read_length_bp",
        "sample_min_read_length_bp",
        "sample_max_read_length_bp",
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
        accession = get_accession(community.get("accession", ""))

        if not accession:
            continue

        meta = assemblies.get(accession, {})
        tax = taxonomy.get(accession, {})
        cov = coverage.get(accession, {})
        read = read_statistics.get(accession, {})

        name = (
            tax.get("name")
            or meta.get("biosample.description.organism_name")
            or "NA"
        )

        abundance = community.get("seq_abundance", "NA")

        try:
            abundance_value = float(abundance)
            abundance_percent = f"{abundance_value * 100:.3f}"
            abundance_fraction = f"{abundance_value:.6f}"
        except (TypeError, ValueError):
            abundance_percent = "NA"
            abundance_fraction = "NA"

        hit = reverse_best.get(accession, {})
        candidate = hit.get("bin", "NA")

        reciprocal = (
            candidate != "NA"
            and forward_best.get(candidate, {}).get("reference")
            == accession
        )

        ani = hit.get("ani", 0)
        fraction = hit.get("fraction", 0)

        identifiable = int(
            reciprocal
            and ani >= minimum_ani
            and fraction >= minimum_fraction
        )

        # MAG metrics are attached only to references passing the rule.
        cm = checkm.get(candidate, {}) if identifiable else {}

        quast = {}

        if identifiable:
            report_file = find_quast_report(
                quast_dir, data_root, candidate, accession
            )

            if report_file:
                quast = read_quast_report(report_file)

        reference_size = (
            meta.get("total_sequence_length")
            or cov.get("seq_len")
            or "NA"
        )

        output_rows.append({
            "sample_id": community.get(
                "sample", args.sample_id
            ),
            "bacterium_name": name,
            "reference_accession": accession,
            "relative_abundance_fraction": abundance_fraction,
            "relative_abundance_percent": abundance_percent,
            "MeSS_reads_estimate": number(
                read.get("MeSS_reads_estimate")
                or cov.get("reads")
            ),
            "simulated_reads_from_MeSS_BAM": (
                read.get("primary_mapped_reads", "NA")
            ),
            "sample_total_reads": sample_total_reads,
            "sample_total_bases_bp": str(fastq_stats["total_bases"]),
            "sample_mean_read_length_bp": number(
                fastq_stats["mean_length"], 1
            ),
            "sample_min_read_length_bp": str(fastq_stats["min_length"]),
            "sample_max_read_length_bp": str(fastq_stats["max_length"]),
            "MeSS_bases": (
                read.get("MeSS_bases")
                or cov.get("bases")
                or "NA"
            ),
            "binning_tool": binning_tool,
            "best_candidate_bin": candidate,
            "FastANI_ANI_percent": number(hit.get("ani")),
            "FastANI_reference_fragment_fraction_percent":
                number(hit.get("fraction")),
            "identifiable_0_1": str(identifiable),
            "CheckM2_completeness_percent":
                number(cm.get("Completeness")),
            "CheckM2_contamination_percent":
                number(cm.get("Contamination")),
            "CheckM2_strain_heterogeneity":
                cm.get("Strain heterogeneity", "NA"),
            "QUAST_genome_fraction_percent":
                quast.get("Genome fraction (%)", "NA"),
            "QUAST_contig_count":
                quast.get("# contigs", "NA"),
            "QUAST_total_assembly_length_bp":
                quast.get("Total length", "NA"),
            "QUAST_N50_bp": quast.get("N50", "NA"),
            "QUAST_NG50_bp": quast.get("NG50", "NA"),
            "reference_genome_size_bp": reference_size,
            "reference_GC_percent":
                meta.get("gc_percent", "NA"),
            "sequencing_technology":
                simulation.get("technology", "NA").upper(),
            "simulation_model":
                simulation.get("model", "NA"),
            "target_mean_read_length_bp":
                simulation.get("mean_length_bp", "NA"),
            "simulation_accuracy":
                simulation.get("accuracy", "NA"),
            "simulation_ratio":
                simulation.get("ratio", "NA"),
            "global_seed":
                simulation.get("seed", "NA"),
        })

    output = args.output
    if not output.is_absolute():
        output = REPO_ROOT / output

    output.parent.mkdir(parents=True, exist_ok=True)

    with output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=output_columns,
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(output_rows)

    identified_count = sum(
        int(row["identifiable_0_1"]) for row in output_rows
    )

    print(f"Output: {output}")
    print(f"Reference rows: {len(output_rows)}")
    print(f"Identifiable references: {identified_count}/{len(output_rows)}")
    print(f"FASTQ source: {fastq_file}")
    print(f"Actual FASTQ read count: {sample_total_reads}")
    print(f"Total FASTQ bases: {fastq_stats['total_bases']}")
    print(f"Mean read length: {fastq_stats['mean_length']:.1f} bp")
    print(f"Minimum read length: {fastq_stats['min_length']} bp")
    print(f"Maximum read length: {fastq_stats['max_length']} bp")
    print(f"Sum of primary BAM-mapped reads: {bam_primary_total}")
    print(f"CheckM2 source: {checkm_file or 'not found'}")
    print(f"FastANI source: {forward_file.parent}")
    print("Summary generation completed.")


if __name__ == "__main__":
    main()
