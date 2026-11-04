#!/usr/bin/env python3

import argparse
import csv
import re
import shlex
import shutil
import subprocess
from pathlib import Path


def read_tsv(path):
    if not path.is_file():
        raise FileNotFoundError(f"Input file not found: {path}")

    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


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


def find_bin_file(bins_dir, bin_name):
    suffixes = (
        ".fa.gz", ".fasta.gz", ".fna.gz",
        ".fa", ".fasta", ".fna",
    )

    for suffix in suffixes:
        candidate = bins_dir / f"{bin_name}{suffix}"
        if candidate.is_file():
            return candidate

    return bins_dir / f"{bin_name}.fa.gz"


def collect_references(references_dir):
    suffixes = (".fna.gz", ".fa.gz", ".fasta.gz", ".fna", ".fa", ".fasta")
    references = {}

    if not references_dir.is_dir():
        return references

    for path in references_dir.rglob("*"):
        if not path.is_file():
            continue

        if not path.name.endswith(suffixes):
            continue

        accession = get_accession(path.name)

        if accession:
            if accession in references:
                raise ValueError(
                    f"Multiple reference files found for {accession}"
                )
            references[accession] = path

    return references


def parse_forward_results(path):
    """Best reference for each bin."""
    best = {}

    with path.open(encoding="utf-8") as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")

            if len(fields) < 5:
                continue

            bin_name = get_bin_name(fields[0])
            accession = get_accession(fields[1])

            try:
                ani = float(fields[2])
            except ValueError:
                continue

            if not accession:
                continue

            current = best.get(bin_name)

            if current is None or ani > current["ani"]:
                best[bin_name] = {
                    "reference": accession,
                    "ani": ani,
                }

    return best


def parse_reverse_results(path):
    """Best bin for each reference, including reference fragment fraction."""
    best = {}

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

            if not accession:
                continue

            fraction = 100 * matched / total if total else 0

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


def read_quast_report(report_file):
    metrics = {}

    with report_file.open(newline="", encoding="utf-8") as handle:
        for row in csv.reader(handle, delimiter="\t"):
            if len(row) >= 2:
                metrics[row[0].strip()] = row[1].strip()

    return metrics


def main():
    parser = argparse.ArgumentParser(
        description="Run QUAST using reciprocal FastANI matches."
    )

    parser.add_argument(
        "--data-root",
        type=Path,
        default=Path.home() / "szakdolgozat" / "MeSS",
        help="Root directory containing simulation and reference data.",
    )

    parser.add_argument(
        "--sample-id",
        help="Run only the selected sample.",
    )

    parser.add_argument(
        "--execute",
        action="store_true",
        help="Execute QUAST. Without this flag, only display commands.",
    )

    args = parser.parse_args()

    repo_root = Path(__file__).resolve().parents[1]
    config_file = repo_root / "config/quast_parameters.tsv"
    data_root = args.data_root.expanduser().resolve()

    if not data_root.is_dir():
        raise NotADirectoryError(
            f"Data directory not found: {data_root}"
        )

    rows = read_tsv(config_file)

    if args.sample_id:
        rows = [
            row for row in rows
            if row["sample_id"] == args.sample_id
        ]

        if not rows:
            raise ValueError(
                f"Sample ID not found: {args.sample_id}"
            )

    if args.execute and not shutil.which("micromamba"):
        raise RuntimeError("micromamba was not found in PATH.")

    for row in rows:
        sample_id = row["sample_id"]

        bins_dir = resolve_path(row["bins_dir"], data_root)
        references_dir = resolve_path(
            row["references_dir"], data_root
        )
        forward_file = resolve_path(
            row["forward_results"], data_root
        )
        reverse_file = resolve_path(
            row["reverse_results"], data_root
        )
        output_dir = resolve_path(row["output_dir"], data_root)

        threads = int(row["threads"])
        minimum_ani = float(row["minimum_ani"])
        minimum_fraction = float(
            row["minimum_reference_fragment_fraction"]
        )

        if threads < 1:
            raise ValueError("threads must be at least 1.")

        if not forward_file.is_file():
            raise FileNotFoundError(
                f"Forward FastANI file not found: {forward_file}"
            )

        if not reverse_file.is_file():
            raise FileNotFoundError(
                f"Reverse FastANI file not found: {reverse_file}"
            )

        references = collect_references(references_dir)
        forward_best = parse_forward_results(forward_file)
        reverse_best = parse_reverse_results(reverse_file)

        summary_rows = []
        commands = []

        print(f"\nSample: {sample_id}")
        print(f"Bin directory exists: {bins_dir.is_dir()}")
        print(f"Reference directory exists: {references_dir.is_dir()}")
        print(f"Reference genomes found: {len(references)}")
        print(f"Forward FastANI results: {forward_file}")
        print(f"Reverse FastANI results: {reverse_file}")
        print(f"ANI threshold: {minimum_ani}%")
        print(f"Reference fragment threshold: {minimum_fraction}%")

        if not references:
            raise ValueError("No reference genomes were found.")

        for accession, reference_file in sorted(references.items()):
            best = reverse_best.get(accession)

            result = {
                "Reference": accession,
                "Bin": "NA",
                "FastANI_ANI_percent": "NA",
                "Reference_fragment_fraction_percent": "NA",
                "Reciprocal_best_hit": 0,
                "Identifiable_0_1": 0,
                "Genome_fraction_percent": "NA",
                "Total_length_bp": "NA",
                "Contig_count": "NA",
                "N50_bp": "NA",
                "NG50_bp": "NA",
                "Misassemblies": "NA",
            }

            if best is None:
                print(f"\n{accession}: no FastANI hit.")
                summary_rows.append(result)
                continue

            bin_name = best["bin"]
            forward_hit = forward_best.get(bin_name)

            reciprocal = (
                forward_hit is not None
                and forward_hit["reference"] == accession
            )

            identifiable = (
                reciprocal
                and best["ani"] >= minimum_ani
                and best["fraction"] >= minimum_fraction
            )

            result.update({
                "Bin": bin_name,
                "FastANI_ANI_percent": f"{best['ani']:.4f}",
                "Reference_fragment_fraction_percent":
                    f"{best['fraction']:.2f}",
                "Reciprocal_best_hit": int(reciprocal),
                "Identifiable_0_1": int(identifiable),
            })

            print(
                f"\n{accession}: {bin_name}, "
                f"ANI={best['ani']:.3f}%, "
                f"fragment fraction={best['fraction']:.2f}%, "
                f"reciprocal={reciprocal}, "
                f"identifiable={int(identifiable)}"
            )

            if not identifiable:
                print("QUAST skipped: reference is not identifiable.")
                summary_rows.append(result)
                continue

            bin_file = find_bin_file(bins_dir, bin_name)

            pair_output = (
                output_dir / "reports" /
                f"{bin_name}_vs_{accession}"
            )

            command = [
                "micromamba", "run", "-n", row["environment"],
                "quast.py",
                "-r", str(reference_file),
                str(bin_file),
                "--threads", str(threads),
                "-o", str(pair_output),
            ]

            print(f"Bin file: {bin_file}")
            print(f"Bin exists: {bin_file.is_file()}")
            print(f"Reference file: {reference_file}")

            commands.append({
                "command": command,
                "bin_file": bin_file,
                "reference_file": reference_file,
                "pair_output": pair_output,
                "result": result,
            })

            if args.execute:
                if not bin_file.is_file():
                    raise FileNotFoundError(
                        f"Bin file not found: {bin_file}. Run SemiBin2 first."
                    )

                if pair_output.exists() and any(pair_output.iterdir()):
                    raise FileExistsError(
                        f"QUAST output directory is not empty: {pair_output}"
                    )

                pair_output.parent.mkdir(parents=True, exist_ok=True)

                print("Starting QUAST...")
                subprocess.run(
                    command,
                    cwd=data_root,
                    check=True,
                )

                report_file = pair_output / "report.tsv"

                if not report_file.is_file():
                    raise RuntimeError(
                        f"QUAST report not found: {report_file}"
                    )

                metrics = read_quast_report(report_file)

                result.update({
                    "Genome_fraction_percent":
                        metrics.get("Genome fraction (%)", "NA"),
                    "Total_length_bp":
                        metrics.get("Total length", "NA"),
                    "Contig_count":
                        metrics.get("# contigs", "NA"),
                    "N50_bp": metrics.get("N50", "NA"),
                    "NG50_bp": metrics.get("NG50", "NA"),
                    "Misassemblies":
                        metrics.get("# misassemblies", "NA"),
                })

            summary_rows.append(result)

        identified_count = sum(
            item["Identifiable_0_1"] for item in summary_rows
        )

        print(
            f"\nIdentifiable references under project rule: "
            f"{identified_count}/{len(summary_rows)}"
        )

        if not args.execute:
            print("\nPrepared QUAST commands:")

            for item in commands:
                print(shlex.join(item["command"]))

            print(
                "\nDRY RUN completed. No QUAST analysis was performed."
            )
            continue

        output_dir.mkdir(parents=True, exist_ok=True)
        summary_file = output_dir / "quast_summary.tsv"

        columns = list(summary_rows[0].keys())

        with summary_file.open(
            "w", newline="", encoding="utf-8"
        ) as handle:
            writer = csv.DictWriter(
                handle,
                fieldnames=columns,
                delimiter="\t",
            )
            writer.writeheader()
            writer.writerows(summary_rows)

        print(f"\nQUAST summary: {summary_file}")


if __name__ == "__main__":
    main()
