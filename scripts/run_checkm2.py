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
        raise ValueError("The CheckM2 parameter file is empty.")

    return rows


def resolve_path(value, data_root):
    path = Path(value).expanduser()

    if not path.is_absolute():
        path = data_root / path

    return path.resolve()


def main():
    parser = argparse.ArgumentParser(
        description="Run CheckM2 using a TSV parameter table."
    )

    parser.add_argument(
        "--data-root",
        type=Path,
        default=Path.home() / "szakdolgozat" / "MeSS",
        help="Root directory containing the binning results.",
    )

    parser.add_argument(
        "--sample-id",
        help="Run only the selected sample.",
    )

    parser.add_argument(
        "--execute",
        action="store_true",
        help="Execute CheckM2. Without this flag, only display the command.",
    )

    args = parser.parse_args()

    repo_root = Path(__file__).resolve().parents[1]
    config_file = repo_root / "config/checkm2_parameters.tsv"
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

        bins_dir = resolve_path(
            row["bins_dir"], data_root
        )
        database = resolve_path(
            row["database_path"], data_root
        )
        output_dir = resolve_path(
            row["output_dir"], data_root
        )

        extension = row["extension"].lstrip(".")
        threads = int(row["threads"])

        if threads < 1:
            raise ValueError("threads must be at least 1.")

        if not extension:
            raise ValueError("File extension must be specified.")

        bin_files = (
            list(bins_dir.glob(f"*.{extension}"))
            if bins_dir.is_dir()
            else []
        )

        command = [
            "micromamba", "run",
            "-n", row["environment"],
            "checkm2", "predict",
            "--input", str(bins_dir),
            "--output-directory", str(output_dir),
            "--threads", str(threads),
            "--database_path", str(database),
            "--extension", extension,
        ]

        print(f"\nSample: {sample_id}")
        print(f"Input bins: {bins_dir}")
        print(f"Input directory exists: {bins_dir.is_dir()}")
        print(f"Matching bin files: {len(bin_files)}")
        print(f"Database: {database}")
        print(f"Database exists: {database.is_file()}")
        print(f"Output: {output_dir}")

        if not args.execute:
            print("\nDRY RUN - CheckM2 was not started.")
            print(shlex.join(command))
            continue

        if not bins_dir.is_dir():
            raise FileNotFoundError(
                f"Bin directory not found: {bins_dir}. "
                "Run SemiBin2 first."
            )

        if not bin_files:
            raise FileNotFoundError(
                f"No files with extension .{extension} found in {bins_dir}"
            )

        if not database.is_file():
            raise FileNotFoundError(
                f"CheckM2 database not found: {database}"
            )

        if output_dir.exists() and any(output_dir.iterdir()):
            raise FileExistsError(
                f"Output directory is not empty: {output_dir}. "
                "Choose a new output directory before rerunning."
            )

        output_dir.parent.mkdir(parents=True, exist_ok=True)

        print("Starting CheckM2...")

        subprocess.run(
            command,
            cwd=data_root,
            check=True,
        )

        report = output_dir / "quality_report.tsv"

        if not report.is_file():
            raise RuntimeError(
                f"CheckM2 finished, but quality_report.tsv "
                f"was not found in {output_dir}"
            )

        print(f"CheckM2 completed: {sample_id}")
        print(f"Quality report: {report}")

    if not args.execute:
        print("\nDry run completed. No QC analysis was performed.")


if __name__ == "__main__":
    main()
