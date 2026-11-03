#!/usr/bin/env python3

import argparse
import csv
import shlex
import shutil
import subprocess
from pathlib import Path


def read_parameters(config_file):
    """Read simulation parameters from a TSV file."""
    with config_file.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = list(reader)

    if not rows:
        raise ValueError("The simulation parameter file is empty.")

    return rows


def parse_bool(value):
    """Convert a TSV boolean value to Python bool."""
    return str(value).strip().lower() in ("true", "1", "yes")


def build_command(row, mess_root):
    """Build the MeSS command from one configuration row."""

    community = Path(row["community_file"])
    if not community.is_absolute():
        community = mess_root / community

    output = Path(row["output_dir"])
    if not output.is_absolute():
        output = mess_root / output

    if not community.is_file():
        raise FileNotFoundError(
            f"Community file not found: {community}"
        )

    command = [
        "micromamba", "run", "-n", "mess",
        "mess", "run",
        "-i", str(community),
        "--accession",
        "--tech", row["technology"],
        "--model", row["model"],
        "--mean-len", str(row["mean_length_bp"]),
        "--accuracy", str(row["accuracy"]),
        "--ratio", row["ratio"],
        "--seed", str(row["seed"]),
        "--threads", str(row["threads"]),
        "-o", str(output),
    ]

    if parse_bool(row["save_bam"]):
        command.append("--bam")

    return command, output


def main():
    parser = argparse.ArgumentParser(
        description="Run MeSS simulations from a TSV parameter table."
    )

    parser.add_argument(
        "--mess-root",
        type=Path,
        default=Path.home() / "szakdolgozat" / "MeSS",
        help="Path to the original MeSS project directory.",
    )

    parser.add_argument(
        "--config",
        type=Path,
        default=Path("config/simulation_parameters.tsv"),
        help="Simulation parameter TSV file.",
    )

    parser.add_argument(
        "--sample-id",
        help="Run only the selected sample.",
    )

    parser.add_argument(
        "--execute",
        action="store_true",
        help="Actually execute simulations. Without this option, only commands are displayed.",
    )

    args = parser.parse_args()

    mess_root = args.mess_root.expanduser().resolve()
    config_file = args.config.expanduser().resolve()

    if not mess_root.is_dir():
        raise NotADirectoryError(
            f"MeSS directory not found: {mess_root}"
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
        print(f"\nSample: {row['sample_id']}")

        command, output = build_command(row, mess_root)

        if args.execute:
            if output.exists() and any(output.iterdir()):
                raise FileExistsError(
                    f"Output directory already contains files: {output}. "
                    "Choose a new output_dir in the parameter table."
                )

            print("Starting MeSS simulation...")
            subprocess.run(
                command,
                cwd=mess_root,
                check=True,
            )
            print(f"Finished: {row['sample_id']}")
        else:
            print("DRY RUN - simulation was not started.")
            print(shlex.join(command))

    if not args.execute:
        print("\nDry run completed. No simulation data were generated.")


if __name__ == "__main__":
    main()
