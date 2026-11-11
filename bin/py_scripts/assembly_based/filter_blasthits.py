#!/usr/bin/env python3
"""Filter fixed-schema BLASTN TSV results and write deterministic summaries."""

from __future__ import annotations

import argparse
import csv
import math
from collections import defaultdict
from pathlib import Path
from typing import Any


VERSION = "1.0.0"
BLAST_COLUMNS = (
    "qseqid",
    "sseqid",
    "pident",
    "length",
    "mismatch",
    "gapopen",
    "qstart",
    "qend",
    "sstart",
    "send",
    "evalue",
    "bitscore",
    "slen",
)
FILTERED_COLUMNS = (*BLAST_COLUMNS, "subject_coverage")


def percentage(value: float | str, name: str) -> float:
    """Return a finite percentage in the inclusive range 0..100."""
    try:
        parsed = float(value)
    except (TypeError, ValueError) as error:
        raise ValueError(f"{name} must be a number between 0 and 100") from error
    if not math.isfinite(parsed) or not 0 <= parsed <= 100:
        raise ValueError(f"{name} must be a finite number between 0 and 100")
    return parsed


def positive_integer(value: int | str, name: str) -> int:
    """Return a strictly positive integer."""
    try:
        parsed = int(value)
    except (TypeError, ValueError) as error:
        raise ValueError(f"{name} must be a positive integer") from error
    if parsed < 1:
        raise ValueError(f"{name} must be a positive integer")
    return parsed


def intervals_overlap(
    first: tuple[str, str, int, int],
    second: tuple[str, str, int, int],
    minimum_overlap_bp: int,
) -> bool:
    """Return whether two hits to the same subject share enough query bases."""
    if first[:2] != second[:2]:
        return False
    overlap_bp = min(first[3], second[3]) - max(first[2], second[2]) + 1
    return overlap_bp >= minimum_overlap_bp


def parse_candidate(
    row: list[str],
    row_number: int,
    minimum_identity: float,
    minimum_coverage: float,
) -> dict[str, Any] | None:
    """Validate one fixed-schema row and return it when thresholds are met."""
    if len(row) != len(BLAST_COLUMNS):
        raise ValueError(
            f"BLAST row {row_number} has {len(row)} columns; "
            f"expected {len(BLAST_COLUMNS)} ({', '.join(BLAST_COLUMNS)})"
        )
    qseqid, sseqid = row[0], row[1]
    if not qseqid or not sseqid:
        raise ValueError(f"BLAST row {row_number} has an empty query or subject identifier")
    try:
        pident = float(row[2])
        alignment_length = int(row[3])
        qstart = int(row[6])
        qend = int(row[7])
        bitscore = float(row[11])
        subject_length = int(row[12])
    except ValueError as error:
        raise ValueError(f"BLAST row {row_number} contains an invalid numeric value") from error
    if not all(math.isfinite(value) for value in (pident, bitscore)):
        raise ValueError(f"BLAST row {row_number} contains a non-finite numeric value")
    if not 0 <= pident <= 100:
        raise ValueError(f"BLAST row {row_number} has pident outside 0..100")
    if alignment_length < 1 or subject_length < 1:
        raise ValueError(
            f"BLAST row {row_number} requires positive alignment length and slen"
        )
    if qstart < 1 or qend < 1:
        raise ValueError(f"BLAST row {row_number} requires positive query coordinates")

    coverage = min(100.0, alignment_length / subject_length * 100.0)
    if pident < minimum_identity or coverage < minimum_coverage:
        return None
    hit_key = (qseqid, sseqid, min(qstart, qend), max(qstart, qend))
    return {
        "row": row,
        "hit_key": hit_key,
        "pident": pident,
        "coverage": coverage,
        "bitscore": bitscore,
        "alignment_length": alignment_length,
    }


def filter_results(
    blast_results: str | Path,
    minimum_identity: float | str,
    minimum_coverage: float | str,
    minimum_overlap_bp: int | str = 1,
) -> tuple[list[list[str]], dict[str, int]]:
    """Select threshold-passing nonduplicate hits independently of row order."""
    identity = percentage(minimum_identity, "minimum identity")
    coverage = percentage(minimum_coverage, "minimum coverage")
    overlap_bp = positive_integer(minimum_overlap_bp, "minimum overlap")
    candidates: list[dict[str, Any]] = []
    with Path(blast_results).open("r", encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        for row_number, row in enumerate(reader, start=1):
            if not row:
                continue
            candidate = parse_candidate(row, row_number, identity, coverage)
            if candidate is not None:
                candidates.append(candidate)

    candidates.sort(
        key=lambda item: (
            -item["pident"],
            -item["coverage"],
            -item["bitscore"],
            -item["alignment_length"],
            item["hit_key"],
            tuple(item["row"]),
        )
    )
    selected: list[dict[str, Any]] = []
    selected_keys: list[tuple[str, str, int, int]] = []
    for candidate in candidates:
        if any(
            intervals_overlap(existing, candidate["hit_key"], overlap_bp)
            for existing in selected_keys
        ):
            continue
        selected.append(candidate)
        selected_keys.append(candidate["hit_key"])

    selected.sort(key=lambda item: (item["hit_key"], tuple(item["row"])))
    filtered_rows: list[list[str]] = []
    gene_counts: defaultdict[str, int] = defaultdict(int)
    for candidate in selected:
        coverage_text = f'{candidate["coverage"]:.6f}'.rstrip("0").rstrip(".")
        filtered_rows.append([*candidate["row"], coverage_text])
        gene_counts[candidate["row"][1]] += 1
    return filtered_rows, dict(gene_counts)


def write_outputs(
    filtered_rows: list[list[str]], gene_counts: dict[str, int], output: str | Path
) -> None:
    """Write escaped, deterministic detailed and gene-count TSV files."""
    output_prefix = str(output)
    with open(
        output_prefix + ".blastn_filtered_hits.tsv",
        "w",
        encoding="utf-8",
        newline="",
    ) as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(FILTERED_COLUMNS)
        writer.writerows(filtered_rows)

    with open(
        output_prefix + ".blastn_summary.tsv",
        "w",
        encoding="utf-8",
        newline="",
    ) as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(("gene", "count"))
        for gene in sorted(gene_counts):
            writer.writerow((gene, gene_counts[gene]))


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--blast-results", required=True, type=Path)
    parser.add_argument("--min-identity", required=True, type=float)
    parser.add_argument("--min-coverage", required=True, type=float)
    parser.add_argument("--overlap-bp", default=1, type=int)
    parser.add_argument("--output", required=True)
    parser.add_argument("--version", action="version", version=VERSION)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    rows, counts = filter_results(
        args.blast_results,
        args.min_identity,
        args.min_coverage,
        args.overlap_bp,
    )
    write_outputs(rows, counts, args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
