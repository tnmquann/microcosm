#!/usr/bin/env python3
"""Validate and normalize the standalone bacterial-genome input CSV."""

from __future__ import annotations

import argparse
import csv
import os
import re
import sys
from pathlib import Path
from typing import Iterable

VERSION = "0.2.0"
REQUIRED_HEADERS = ("sample_id", "fasta")
OPTIONAL_HEADERS = ("run_id", "group", "assembler", "lineage")
OUTPUT_HEADERS = (
    "sample_id",
    "run_id",
    "group",
    "assembler",
    "fasta",
    "lineage",
    "domain",
    "phylum",
    "class",
    "order",
    "family",
    "genus",
    "species",
    "lineage_source",
    "classification_key",
    "source_route",
    "row_number",
)
RANK_ORDER = ("d", "p", "c", "o", "f", "g", "s")
SUPPORTED_FASTA_SUFFIXES = (
    ".fa",
    ".fasta",
    ".fna",
    ".fas",
    ".fa.gz",
    ".fasta.gz",
    ".fna.gz",
    ".fas.gz",
)
SAFE_SAMPLE_ID = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]*$")


class ValidationError(ValueError):
    """Raised when a user-facing input validation check fails."""


def _unchanged_trimmed(value: str, row_number: int, field: str, sample: str) -> str:
    """Reject surrounding whitespace instead of silently changing user metadata."""
    if value != value.strip():
        raise ValidationError(
            f"row {row_number} sample '{sample or '<unknown>'}': field '{field}' "
            "contains leading or trailing whitespace"
        )
    return value


def parse_lineage(lineage: str, context: str) -> tuple[str, ...]:
    """Validate ordered GTDB-style rank tokens and return routing ranks."""
    if "\t" in lineage or "\n" in lineage or "\r" in lineage:
        raise ValidationError(f"{context}: lineage contains a tab or newline")

    tokens = lineage.split(";")
    if not tokens or any(token == "" for token in tokens):
        raise ValidationError(f"{context}: lineage contains an empty rank token")
    if len(tokens) > len(RANK_ORDER):
        raise ValidationError(
            f"{context}: lineage has {len(tokens)} ranks; at most 7 are supported"
        )

    parsed: dict[str, str] = {}
    previous_index = -1
    for token in tokens:
        match = re.fullmatch(r"([dpcofgs])__(.*)", token)
        if not match:
            raise ValidationError(
                f"{context}: invalid lineage token '{token}'; expected GTDB form such as d__Bacteria"
            )
        rank, value = match.groups()
        rank_index = RANK_ORDER.index(rank)
        if rank_index <= previous_index:
            raise ValidationError(
                f"{context}: lineage ranks are duplicated or out of canonical order at '{token}'"
            )
        if rank in parsed:
            raise ValidationError(f"{context}: duplicate lineage rank '{rank}__'")
        parsed[rank] = f"{rank}__{value}"
        previous_index = rank_index

    domain = parsed.get("d", "")
    if not domain or domain == "d__":
        raise ValidationError(f"{context}: lineage must contain a non-empty domain rank")
    return tuple(parsed.get(rank, "") for rank in RANK_ORDER)


def fasta_stem(path: Path) -> str:
    """Remove one supported FASTA suffix, including optional gzip."""
    name = path.name
    lower_name = name.lower()
    for suffix in sorted(SUPPORTED_FASTA_SUFFIXES, key=len, reverse=True):
        if lower_name.endswith(suffix):
            return name[: -len(suffix)]
    raise ValidationError(
        f"unsupported FASTA extension for '{path}'; supported extensions: "
        + ", ".join(SUPPORTED_FASTA_SUFFIXES)
    )


def validate_rows(
    input_csv: Path, classification_available: bool
) -> list[dict[str, str]]:
    """Validate all rows and return normalized output records."""
    with input_csv.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None:
            raise ValidationError("input CSV has no header row")

        normalized_headers = [header.strip() for header in reader.fieldnames]
        if len(normalized_headers) != len(set(normalized_headers)):
            raise ValidationError("input CSV contains duplicate header names")
        missing_headers = [
            header for header in REQUIRED_HEADERS if header not in normalized_headers
        ]
        if missing_headers:
            raise ValidationError(
                "input CSV is missing required header(s): " + ", ".join(missing_headers)
            )
        reader.fieldnames = normalized_headers

        output_rows: list[dict[str, str]] = []
        sample_rows: dict[str, int] = {}
        missing_lineage_rows: list[tuple[int, str]] = []
        classification_keys: dict[str, tuple[int, str]] = {}

        for row_number, row in enumerate(reader, start=2):
            if None in row:
                raise ValidationError(
                    f"row {row_number}: contains more values than the header defines"
                )
            sample_id = _unchanged_trimmed(
                row.get("sample_id") or "", row_number, "sample_id", ""
            )
            if not sample_id:
                raise ValidationError(
                    f"row {row_number}: required field 'sample_id' is empty"
                )
            if not SAFE_SAMPLE_ID.fullmatch(sample_id):
                raise ValidationError(
                    f"row {row_number} sample '{sample_id}': sample_id must match "
                    "[A-Za-z0-9][A-Za-z0-9._-]* for collision-safe output paths"
                )
            if sample_id in sample_rows:
                raise ValidationError(
                    f"row {row_number} sample '{sample_id}': duplicate sample_id; "
                    f"first seen on row {sample_rows[sample_id]}"
                )
            sample_rows[sample_id] = row_number

            raw_fasta = _unchanged_trimmed(
                row.get("fasta") or "", row_number, "fasta", sample_id
            )
            if not raw_fasta:
                raise ValidationError(
                    f"row {row_number} sample '{sample_id}': required field 'fasta' is empty"
                )
            fasta_path = Path(raw_fasta).expanduser()
            if not fasta_path.is_absolute():
                fasta_path = input_csv.parent / fasta_path
            fasta_path = fasta_path.resolve()
            classification_key = fasta_stem(fasta_path)
            if not fasta_path.is_file():
                raise ValidationError(
                    f"row {row_number} sample '{sample_id}': FASTA does not exist: {fasta_path}"
                )
            if not os.access(fasta_path, os.R_OK):
                raise ValidationError(
                    f"row {row_number} sample '{sample_id}': FASTA is not readable: {fasta_path}"
                )

            optional_values: dict[str, str] = {}
            for field in OPTIONAL_HEADERS:
                optional_values[field] = _unchanged_trimmed(
                    row.get(field) or "", row_number, field, sample_id
                )

            lineage = optional_values["lineage"]
            if lineage:
                domain, phylum, class_rank, order, family, genus, species = parse_lineage(
                    lineage, f"row {row_number} sample '{sample_id}' field 'lineage'"
                )
                lineage_source = "provided"
            else:
                domain, phylum, class_rank, order, family, genus, species = (
                    "",
                    "",
                    "",
                    "",
                    "",
                    "",
                    "",
                )
                lineage_source = "classified"
                missing_lineage_rows.append((row_number, sample_id))
                if classification_key in classification_keys:
                    first_row, first_sample = classification_keys[classification_key]
                    raise ValidationError(
                        f"row {row_number} sample '{sample_id}': FASTA basename "
                        f"'{classification_key}' collides with row {first_row} sample "
                        f"'{first_sample}' among rows requiring classification"
                    )
                classification_keys[classification_key] = (row_number, sample_id)

            output_rows.append(
                {
                    "sample_id": sample_id,
                    "run_id": optional_values["run_id"],
                    "group": optional_values["group"],
                    "assembler": optional_values["assembler"],
                    "fasta": str(fasta_path),
                    "lineage": lineage,
                    "domain": domain,
                    "phylum": phylum,
                    "class": class_rank,
                    "order": order,
                    "family": family,
                    "genus": genus,
                    "species": species,
                    "lineage_source": lineage_source,
                    "classification_key": classification_key,
                    "source_route": "standalone",
                    "row_number": str(row_number),
                }
            )

    if not output_rows:
        raise ValidationError("input CSV contains no data rows")
    if missing_lineage_rows and not classification_available:
        details = ", ".join(
            f"row {row_number} sample '{sample}'"
            for row_number, sample in missing_lineage_rows
        )
        raise ValidationError(
            "lineage classification is required for "
            f"{details}, but the bin classifier inputs are missing (--sourmash_database and "
            "--sourmash_taxonomy_csv, or --bin_classifier gtdbtk with --gtdbtk_db)"
        )
    return output_rows


def write_rows(output_csv: Path, rows: Iterable[dict[str, str]]) -> None:
    """Write the deterministic normalized manifest."""
    with output_csv.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=OUTPUT_HEADERS)
        writer.writeheader()
        writer.writerows(rows)


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, type=Path, help="Genome input CSV")
    parser.add_argument(
        "--output", required=True, type=Path, help="Normalized validated CSV"
    )
    parser.add_argument(
        "--classification-available",
        action="store_true",
        help="Confirm that both required classification database inputs are configured",
    )
    parser.add_argument("--version", action="version", version=VERSION)
    return parser.parse_args()


def main() -> int:
    """Run validation and normalization."""
    args = parse_args()
    try:
        input_csv = args.input.expanduser().resolve()
        if not input_csv.is_file():
            raise ValidationError(f"genome input CSV does not exist: {input_csv}")
        rows = validate_rows(input_csv, args.classification_available)
        write_rows(args.output, rows)
    except (OSError, csv.Error, ValidationError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
