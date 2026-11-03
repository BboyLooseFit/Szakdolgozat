#!/usr/bin/env python3

import argparse
import csv
import shlex
import shutil
import subprocess
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]


def read_parameters(config_file):
    with config_file.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))

    if not rows:
        raise ValueError("The assembly parameter file is empty.")

    return rows


def main():
    parser = argparse.ArgumentParser(
        description="Run Flye assembly using a TSV parameter table."
    )

    parser.add_argument(
        "--data-root",
        type=Path,
        default=Path.home() / "szakdolgozat" / "MeSS",
        help="Root directory containing the MeSS simulation results.",
    )

    parser.add_argument(
        "--sample-id",
        help="Run only the selected sample.",
    )

    parser.add_argument(
        "--execute",
        action="store_true",
        help="Execute Flye. Without this flag, only display the command.",
    )

    args = parser.parse_args()

    data_root = args.data_root.expanduser().resolve()
    config_file = REPO_ROOT / "config/assembly_parameters.tsv"

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

        input_reads = Path(row["input_reads"])
        if not input_reads.is_absolute():
            input_reads = data_root / input_reads

        output_dir = Path(row["output_dir"])
        if not output_dir.is_absolute():
            output_dir = data_root / output_dir

        if not input_reads.is_file():
            raise FileNotFoundError(
                f"Input FASTQ not found: {input_reads}"
            )

        read_type = row["read_type"]
        if read_type not in {"nano-hq", "nano-raw"}:
            raise ValueError(
                f"Unsupported read type: {read_type}"
            )

        command = [
            "micromamba", "run", "-n", "flye",
            "flye",
            f"--{read_type}",
            str(input_reads),
        ]

        if row["assembly_mode"].lower() == "meta":
            command.append("--meta")

        command.extend([
            "--read-error", row["read_error"],
            "--threads", row["threads"],
            "--out-dir", str(output_dir),
        ])

        print(f"\nSample: {sample_id}")
        print(f"Input: {input_reads}")
        print(f"Output: {output_dir}")

        if args.execute:
            if output_dir.exists() and any(output_dir.iterdir()):
                raise FileExistsError(
                    f"Output directory is not empty: {output_dir}"
                )

            output_dir.parent.mkdir(parents=True, exist_ok=True)

            print("Starting Flye assembly...")
            subprocess.run(
                command,
                cwd=data_root,
                check=True,
            )
            print(f"Finished: {sample_id}")

        else:
            print("DRY RUN - assembly was not started.")
            print(shlex.join(command))

    if not args.execute:
        print("\nDry run completed. No assembly was generated.")


if __name__ == "__main__":
    main()
