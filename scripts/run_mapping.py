#!/usr/bin/env python3

import argparse
import csv
import shlex
import shutil
import subprocess
from pathlib import Path


def read_parameters(config_file):
    with config_file.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))

    if not rows:
        raise ValueError("The mapping parameter file is empty.")

    return rows


def resolve_path(value, data_root):
    path = Path(value)
    if not path.is_absolute():
        path = data_root / path
    return path.resolve()


def main():
    parser = argparse.ArgumentParser(
        description="Run minimap2 read mapping and samtools processing."
    )

    parser.add_argument(
        "--data-root",
        type=Path,
        default=Path.home() / "szakdolgozat" / "MeSS",
        help="Root directory containing the simulation and assembly data.",
    )

    parser.add_argument(
        "--sample-id",
        help="Run only the selected sample.",
    )

    parser.add_argument(
        "--execute",
        action="store_true",
        help="Execute mapping. Without this option, only display commands.",
    )

    args = parser.parse_args()

    repo_root = Path(__file__).resolve().parents[1]
    config_file = repo_root / "config/mapping_parameters.tsv"
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

    if args.execute:
        for tool in ("micromamba", "samtools"):
            if not shutil.which(tool):
                raise RuntimeError(f"{tool} was not found in PATH.")

    for row in rows:
        sample_id = row["sample_id"]

        reads = resolve_path(row["input_reads"], data_root)
        assembly = resolve_path(row["assembly_fasta"], data_root)
        output_dir = resolve_path(row["output_dir"], data_root)

        threads = int(row["sort_threads"])
        preset = row["preset"]

        if threads < 1:
            raise ValueError("sort_threads must be at least 1.")

        bam_file = output_dir / f"{sample_id}_vs_assembly.bam"
        bai_file = Path(str(bam_file) + ".bai")
        flagstat_file = output_dir / "flagstat.txt"
        coverage_file = output_dir / "coverage.tsv"

        print(f"\nSample: {sample_id}")
        print(f"Reads: {reads} (exists: {reads.is_file()})")
        print(f"Assembly: {assembly} (exists: {assembly.is_file()})")
        print(f"Output BAM: {bam_file}")

        minimap2_command = [
            "micromamba", "run", "-n", "flye",
            "minimap2",
            "-ax", preset,
            str(assembly),
            str(reads),
        ]

        samtools_sort_command = [
            "samtools", "sort",
            "-@", str(threads),
            "-o", str(bam_file),
        ]

        if not args.execute:
            print("DRY RUN - mapping was not started.")
            print("Minimap2:")
            print(shlex.join(minimap2_command))
            print("Samtools sort:")
            print(shlex.join(samtools_sort_command))
            print(f"Samtools index: samtools index {bam_file}")
            print(f"Samtools flagstat: samtools flagstat {bam_file}")
            print(f"Samtools coverage: samtools coverage {bam_file}")
            continue

        if not reads.is_file():
            raise FileNotFoundError(f"FASTQ file not found: {reads}")

        if not assembly.is_file():
            raise FileNotFoundError(
                f"Assembly FASTA not found: {assembly}. "
                "Run the Flye module first."
            )

        if any(path.exists() for path in (
            bam_file, bai_file, flagstat_file, coverage_file
        )):
            raise FileExistsError(
                f"Mapping output already exists for sample {sample_id}. "
                "Choose a new output directory before rerunning."
            )

        output_dir.mkdir(parents=True, exist_ok=True)

        print("Starting minimap2 and samtools pipeline...")

        with subprocess.Popen(
            minimap2_command,
            stdout=subprocess.PIPE,
        ) as mapper:
            sort_result = subprocess.run(
                samtools_sort_command,
                stdin=mapper.stdout,
                check=False,
            )
            mapper.stdout.close()
            mapper_returncode = mapper.wait()

        if mapper_returncode != 0:
            raise subprocess.CalledProcessError(mapper_returncode, minimap2_command)

        if sort_result.returncode != 0:
            raise subprocess.CalledProcessError(
                sort_result.returncode, samtools_sort_command
            )

        subprocess.run(
            ["samtools", "index", str(bam_file)],
            check=True,
        )

        with flagstat_file.open("w") as handle:
            subprocess.run(
                ["samtools", "flagstat", str(bam_file)],
                stdout=handle,
                check=True,
                text=True,
            )

        with coverage_file.open("w") as handle:
            subprocess.run(
                ["samtools", "coverage", str(bam_file)],
                stdout=handle,
                check=True,
                text=True,
            )

        print(f"Mapping completed: {sample_id}")
        print(f"BAM: {bam_file}")
        print(f"Index: {bai_file}")
        print(f"Flagstat: {flagstat_file}")
        print(f"Coverage: {coverage_file}")

    if not args.execute:
        print("\nDry run completed. No mapping was performed.")


if __name__ == "__main__":
    main()
