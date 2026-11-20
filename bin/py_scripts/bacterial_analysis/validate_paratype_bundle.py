#!/usr/bin/env python3
"""Validate a user-supplied, read-only Paratype v1.1-style bundle."""

from __future__ import annotations

import argparse
import os
from pathlib import Path


VERSION = "1.0.0"
REQUIRED_FILES = (
    "paratype.py",
    "SParatyphiAKU12601.fasta",
    "SParatyphiA_genotype_specific_alleles.txt",
    "SParatyphiA_gene_mutation_codons.txt",
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bundle", required=True, type=Path)
    parser.add_argument("--version", action="version", version=VERSION)
    return parser.parse_args()


def validate_bundle(bundle: Path) -> Path:
    """Check that a staged bundle contains the exact files Paratype needs."""
    resolved = bundle.resolve()
    if not resolved.is_dir():
        raise ValueError(f"Paratype bundle is not a directory: {resolved}")
    for filename in REQUIRED_FILES:
        candidate = resolved / filename
        if not candidate.is_file() or not os.access(candidate, os.R_OK):
            raise ValueError(
                f"Paratype bundle is missing a readable required file '{filename}': {resolved}"
            )
        if candidate.stat().st_size == 0:
            raise ValueError(f"Paratype bundle file is empty: {candidate}")
    return resolved


def main() -> int:
    args = parse_args()
    try:
        validate_bundle(args.bundle)
    except ValueError as error:
        raise SystemExit(f"ERROR: {error}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
