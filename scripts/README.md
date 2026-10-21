# Scripts

This directory contains custom scripts developed for the thesis project.

## Current scripts

### `build_10_reference_summary.py`

Generates a 10-reference summary table by combining:

- MeSS community composition and simulation metadata
- Reference genome taxonomy and assembly metadata
- FastANI reference-to-bin comparisons
- CheckM2 MAG quality metrics
- QUAST assembly quality metrics
- Read counts derived from the simulated-read BAM file

The output contains one row for each original reference genome.

## Inputs

The script currently expects the MeSS project directory structure and the corresponding simulation, FastANI, CheckM2 and QUAST output files.

Large raw sequencing files, reference databases and intermediate analysis files are not stored in this repository.

## Reproducibility

The current implementation was tested with one simulated ONT sample generated using seed 42.

The identifiability criterion used in the summary is project-specific and is documented alongside the generated table.
