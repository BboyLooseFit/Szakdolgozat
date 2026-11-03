#!/usr/bin/env python3

import argparse
import csv
import shlex
import shutil
import subprocess
from pathlib import Path


def read_parameters(config_file):
    """Read SemiBin2 parameters from a TSV file."""
    with config_file.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))

    if not rows:
        raise ValueError("The SemiBin2 parameter file is empty.")

    return rows


def resolve_path(value, data_root):
    """Resolve relative paths against the MeSS project directory."""
    path = Path(value)

    if not path.is_absolute():
        path = data_root / path

    return path.resolve()


def main():
    parser = argparse.ArgumentParser(
        description="Run SemiBin2 using a TSV parameter table."
    )

    parser.add_argument(
        "--data-root",
        type=Path,
        default=Path.home() / "szakdolgozat" / "MeSS",
        help="Root directory containing assembly and mapping results.",
    )

    parser.add_argument(
        "--sample-id",
        help="Run only the selected sample.",
    )

    parser.add_argument(
        "--execute",
        action="store_true",
        help="Execute SemiBin2. Without this flag, only display the command.",
    )

    args = parser.parse_args()

    repo_root = Path(__file__).resolve().parents[1]
    config_file = repo_root / "config/semibin2_parameters.tsv"
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

        assembly = resolve_path(
            row["assembly_fasta"], data_root
        )
        bam = resolve_path(
            row["bam_file"], data_root
        )
        output_dir = resolve_path(
            row["output_dir"], data_root
        )

        command = [
            "micromamba", "run", "-n", "semibin",
            "SemiBin2", "single_easy_bin",
            "--environment", row["environment"],
            "--sequencing-type", row["sequencing_type"],
            "-i", str(assembly),
            "-b", str(bam),
            "-o", str(output_dir),
        ]

        print(f"\nSample: {sample_id}")
        print(f"Assembly: {assembly}")
        print(f"Assembly exists: {assembly.is_file()}")
        print(f"BAM: {bam}")
        print(f"BAM exists: {bam.is_file()}")
        print(f"Output: {output_dir}")

        if not args.execute:
            print("\nDRY RUN - binning was not started.")
            print(shlex.join(command))
            continue

        if not assembly.is_file():
            raise FileNotFoundError(
                f"Assembly FASTA not found: {assembly}. "
                "Run Flye first."
            )

        if not bam.is_file():
            raise FileNotFoundError(
                f"BAM file not found: {bam}. "
                "Run read mapping first."
            )

        if output_dir.exists():
            raise FileExistsError(
                f"Output directory already exists: {output_dir}. "
                "Choose a new output directory before rerunning."
            )

        output_dir.parent.mkdir(parents=True, exist_ok=True)

        print("Starting SemiBin2...")

        subprocess.run(
            command,
            cwd=data_root,
            check=True,
        )

        print(f"SemiBin2 completed: {sample_id}")
        print(f"Results: {output_dir}")

    if not args.execute:
        print("\nDry run completed. No binning was performed.")


if __name__ == "__main__":
    main()
