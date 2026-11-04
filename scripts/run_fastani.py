#!/usr/bin/env python3

import argparse
import csv
import re
import shlex
import shutil
import subprocess
from pathlib import Path


def read_parameters(config_file):
    with config_file.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))

    if not rows:
        raise ValueError("The FastANI parameter file is empty.")

    return rows


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


def read_fastani_output(path):
    """Read the standard five-column FastANI output."""
    rows = []

    if not path.is_file():
        return rows

    with path.open(encoding="utf-8") as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")

            if len(fields) < 5:
                continue

            try:
                ani = float(fields[2])
                matched = int(fields[3])
                total = int(fields[4])
            except ValueError:
                continue

            rows.append({
                "query": fields[0],
                "reference": fields[1],
                "ani": ani,
                "matched_fragments": matched,
                "total_fragments": total,
            })

    return rows


def write_recovery_summary(
    forward_rows,
    reverse_rows,
    reference_files,
    output_file,
    minimum_ani,
    minimum_fragment_fraction,
):
    """Identify reciprocal best bin-reference matches."""

    forward_best = {}

    for row in forward_rows:
        bin_name = get_bin_name(row["query"])
        accession = get_accession(row["reference"])

        if not accession:
            continue

        current = forward_best.get(bin_name)

        if current is None or row["ani"] > current["ani"]:
            forward_best[bin_name] = {
                "reference": accession,
                "ani": row["ani"],
            }

    reverse_best = {}

    for row in reverse_rows:
        accession = get_accession(row["query"])
        bin_name = get_bin_name(row["reference"])

        if not accession:
            continue

        fraction = (
            100 * row["matched_fragments"] /
            row["total_fragments"]
            if row["total_fragments"] else 0
        )

        current = reverse_best.get(accession)

        if current is None or row["ani"] > current["ani"]:
            reverse_best[accession] = {
                "bin": bin_name,
                "ani": row["ani"],
                "matched": row["matched_fragments"],
                "total": row["total_fragments"],
                "fraction": fraction,
            }

    columns = [
        "Reference",
        "Best_bin",
        "ANI_percent",
        "Matched_fragments",
        "Total_fragments",
        "Reference_fragment_fraction_percent",
        "Reciprocal_best_hit",
        "Identifiable_0_1",
    ]

    summary = []

    for reference_file in reference_files:
        accession = get_accession(reference_file.name)

        if not accession:
            continue

        best = reverse_best.get(accession)

        if best is None:
            summary.append({
                "Reference": accession,
                "Best_bin": "NA",
                "ANI_percent": "NA",
                "Matched_fragments": "NA",
                "Total_fragments": "NA",
                "Reference_fragment_fraction_percent": "NA",
                "Reciprocal_best_hit": 0,
                "Identifiable_0_1": 0,
            })
            continue

        candidate = best["bin"]
        forward = forward_best.get(candidate)

        reciprocal = int(
            forward is not None
            and forward["reference"] == accession
        )

        identifiable = int(
            reciprocal == 1
            and best["ani"] >= minimum_ani
            and best["fraction"] >= minimum_fragment_fraction
        )

        summary.append({
            "Reference": accession,
            "Best_bin": candidate,
            "ANI_percent": f"{best['ani']:.4f}",
            "Matched_fragments": best["matched"],
            "Total_fragments": best["total"],
            "Reference_fragment_fraction_percent":
                f"{best['fraction']:.2f}",
            "Reciprocal_best_hit": reciprocal,
            "Identifiable_0_1": identifiable,
        })

    with output_file.open(
        "w", newline="", encoding="utf-8"
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=columns,
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(summary)

    return summary


def main():
    parser = argparse.ArgumentParser(
        description="Run two-way FastANI comparisons from a TSV configuration."
    )

    parser.add_argument(
        "--data-root",
        type=Path,
        default=Path.home() / "szakdolgozat" / "MeSS",
        help="Root directory containing bins and reference genomes.",
    )

    parser.add_argument(
        "--sample-id",
        help="Run only the selected sample.",
    )

    parser.add_argument(
        "--execute",
        action="store_true",
        help="Execute FastANI. Without this option, only display commands.",
    )

    args = parser.parse_args()

    repo_root = Path(__file__).resolve().parents[1]
    config_file = repo_root / "config/fastani_parameters.tsv"
    data_root = args.data_root.expanduser().resolve()

    if not data_root.is_dir():
        raise NotADirectoryError(
            f"Data directory not found: {data_root}"
        )

    rows = read_parameters(config_file)

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
        output_dir = resolve_path(row["output_dir"], data_root)

        threads = int(row["threads"])
        minimum_ani = float(row["minimum_ani"])
        minimum_fraction = float(
            row["minimum_reference_fragment_fraction"]
        )

        if threads < 1:
            raise ValueError("threads must be at least 1.")

        bin_files = (
            sorted(bins_dir.glob(row["bin_pattern"]))
            if bins_dir.is_dir() else []
        )

        reference_files = (
            sorted(references_dir.rglob(row["reference_pattern"]))
            if references_dir.is_dir() else []
        )

        query_list = output_dir / "query_bins.txt"
        reference_list = output_dir / "reference_genomes.txt"

        forward_output = output_dir / "bin_vs_reference.tsv"
        reverse_output = output_dir / "reference_vs_bins.tsv"

        summary_output = output_dir / "reference_recovery_summary.tsv"

        print(f"\nSample: {sample_id}")
        print(f"Bin directory exists: {bins_dir.is_dir()}")
        print(f"Bin files found: {len(bin_files)}")
        print(f"Reference directory exists: {references_dir.is_dir()}")
        print(f"Reference genomes found: {len(reference_files)}")
        print(f"ANI threshold: {minimum_ani}%")
        print(f"Reference fragment threshold: {minimum_fraction}%")

        command_forward = [
            "micromamba", "run", "-n", row["environment"],
            "fastANI",
            "--ql", str(query_list),
            "--rl", str(reference_list),
            "--threads", str(threads),
            "-o", str(forward_output),
        ]

        command_reverse = [
            "micromamba", "run", "-n", row["environment"],
            "fastANI",
            "--ql", str(reference_list),
            "--rl", str(query_list),
            "--threads", str(threads),
            "-o", str(reverse_output),
        ]

        if not args.execute:
            print("\nDRY RUN - FastANI was not started.")
            print("\nForward comparison:")
            print(shlex.join(command_forward))
            print("\nReverse comparison:")
            print(shlex.join(command_reverse))
            print(f"\nExpected summary: {summary_output}")
            continue

        if not bin_files:
            raise FileNotFoundError(
                f"No bin files found in {bins_dir}. Run SemiBin2 first."
            )

        if not reference_files:
            raise FileNotFoundError(
                f"No reference genomes found in {references_dir}."
            )

        if output_dir.exists() and any(output_dir.iterdir()):
            raise FileExistsError(
                f"Output directory is not empty: {output_dir}. "
                "Choose a new output directory before rerunning."
            )

        output_dir.mkdir(parents=True, exist_ok=True)

        query_list.write_text(
            "".join(f"{path}\n" for path in bin_files),
            encoding="utf-8",
        )

        reference_list.write_text(
            "".join(f"{path}\n" for path in reference_files),
            encoding="utf-8",
        )

        print("\nRunning bin-to-reference comparison...")
        subprocess.run(command_forward, check=True)

        print("\nRunning reference-to-bin comparison...")
        subprocess.run(command_reverse, check=True)

        forward_rows = read_fastani_output(forward_output)
        reverse_rows = read_fastani_output(reverse_output)

        summary = write_recovery_summary(
            forward_rows,
            reverse_rows,
            reference_files,
            summary_output,
            minimum_ani,
            minimum_fraction,
        )

        identified = sum(
            row["Identifiable_0_1"] for row in summary
        )

        print(f"\nFastANI completed: {sample_id}")
        print(f"Reference genomes: {len(summary)}")
        print(f"Identifiable references: {identified}/{len(summary)}")
        print(f"Summary: {summary_output}")

    if not args.execute:
        print("\nDry run completed. No ANI analysis was performed.")


if __name__ == "__main__":
    main()
