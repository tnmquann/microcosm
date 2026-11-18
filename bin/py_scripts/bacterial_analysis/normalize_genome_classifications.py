#!/usr/bin/env python3
"""Normalize Sourmash genome-classification CSV files for lineage routing."""

from __future__ import annotations

import argparse
import csv
import re
import sys
from pathlib import Path

VERSION = "0.2.0"
RANK_ORDER = ("d", "p", "c", "o", "f", "g", "s")
OUTPUT_HEADERS = (
    "query_name",
    "status",
    "lineage",
    "domain",
    "phylum",
    "class",
    "order",
    "family",
    "genus",
    "species",
)


class ClassificationError(ValueError):
    """Raised when a classification result cannot be routed safely."""


def parse_lineage(lineage: str, context: str) -> tuple[str, ...]:
    """Validate ordered GTDB-style rank tokens and return routing ranks."""
    if "\t" in lineage or "\n" in lineage or "\r" in lineage:
        raise ClassificationError(f"{context}: lineage contains a tab or newline")
    tokens = lineage.split(";")
    if not tokens or any(token == "" for token in tokens):
        raise ClassificationError(f"{context}: lineage contains an empty rank token")
    if len(tokens) > len(RANK_ORDER):
        raise ClassificationError(
            f"{context}: lineage has {len(tokens)} ranks; at most 7 are supported"
        )

    parsed: dict[str, str] = {}
    previous_index = -1
    for token in tokens:
        match = re.fullmatch(r"([dpcofgs])__(.*)", token)
        if not match:
            raise ClassificationError(
                f"{context}: invalid lineage token '{token}'; expected GTDB rank form"
            )
        rank, value = match.groups()
        rank_index = RANK_ORDER.index(rank)
        if rank_index <= previous_index or rank in parsed:
            raise ClassificationError(
                f"{context}: lineage ranks are duplicated or out of canonical order"
            )
        parsed[rank] = f"{rank}__{value}"
        previous_index = rank_index

    domain = parsed.get("d", "")
    if not domain or domain == "d__":
        raise ClassificationError(f"{context}: lineage has no non-empty domain")
    return tuple(parsed.get(rank, "") for rank in RANK_ORDER)


def normalize(
    classification_paths: list[Path],
) -> tuple[list[dict[str, str]], list[dict[str, str]]]:
    """Read classifications, validate uniqueness, and split resolved/unresolved rows."""
    resolved: list[dict[str, str]] = []
    unresolved: list[dict[str, str]] = []
    seen_queries: dict[str, Path] = {}

    for classification_path in classification_paths:
        with classification_path.open(
            "r", encoding="utf-8-sig", newline=""
        ) as handle:
            reader = csv.DictReader(handle)
            if reader.fieldnames is None:
                continue
            headers = [header.strip() for header in reader.fieldnames]
            if len(headers) != len(set(headers)):
                raise ClassificationError(
                    f"{classification_path}: contains duplicate header names"
                )
            reader.fieldnames = headers
            required = [field for field in ("query_name", "lineage") if field not in headers]
            if required:
                raise ClassificationError(
                    f"{classification_path}: missing required header(s): "
                    + ", ".join(required)
                )

            for row_number, row in enumerate(reader, start=2):
                if None in row:
                    raise ClassificationError(
                        f"{classification_path} row {row_number}: contains more values "
                        "than the header defines"
                    )
                query_name = (row.get("query_name") or "").strip()
                status = (row.get("status") or "").strip()
                lineage = (row.get("lineage") or "").strip()
                if not query_name:
                    raise ClassificationError(
                        f"{classification_path} row {row_number}: query_name is empty"
                    )
                if query_name in seen_queries:
                    raise ClassificationError(
                        f"duplicate classification for query '{query_name}' in "
                        f"{classification_path} and {seen_queries[query_name]}"
                    )
                seen_queries[query_name] = classification_path

                if not lineage:
                    unresolved.append(
                        {
                            "query_name": query_name,
                            "status": status,
                            "lineage": "",
                            "domain": "",
                            "phylum": "",
                            "class": "",
                            "order": "",
                            "family": "",
                            "genus": "",
                            "species": "",
                        }
                    )
                    continue

                domain, phylum, class_rank, order, family, genus, species = parse_lineage(
                    lineage, f"{classification_path} row {row_number} query '{query_name}'"
                )
                resolved.append(
                    {
                        "query_name": query_name,
                        "status": status,
                        "lineage": lineage,
                        "domain": domain,
                        "phylum": phylum,
                        "class": class_rank,
                        "order": order,
                        "family": family,
                        "genus": genus,
                        "species": species,
                    }
                )
    return resolved, unresolved


def write_tsv(path: Path, rows: list[dict[str, str]]) -> None:
    """Write a deterministic tab-separated routing table."""
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=OUTPUT_HEADERS, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--classification",
        required=True,
        nargs="+",
        type=Path,
        help="One or more Sourmash tax-genome CSV files",
    )
    parser.add_argument(
        "--resolved-output", required=True, type=Path, help="Resolved lineage TSV"
    )
    parser.add_argument(
        "--unresolved-output", required=True, type=Path, help="Unresolved lineage TSV"
    )
    parser.add_argument("--version", action="version", version=VERSION)
    return parser.parse_args()


def main() -> int:
    """Normalize all classification files."""
    args = parse_args()
    try:
        resolved, unresolved = normalize(args.classification)
        write_tsv(args.resolved_output, resolved)
        write_tsv(args.unresolved_output, unresolved)
    except (OSError, csv.Error, ClassificationError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
