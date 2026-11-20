#!/usr/bin/env python3
"""GTDB-Tk classify_wf summaries -> the sourmash tax genome CSV layout used by microcosm.

The bacterial routing (normalize_genome_classifications.py) and the bin MultiQC tables read
query_name, status and lineage; this keeps one classification contract for both classifiers.
  - query_name = user_genome (GTDB-Tk already removed the --extension)
  - lineage    = classification without empty ranks (g__X;s__ -> ...;g__X), like sourmash
  - "Unclassified Bacteria" / "Unclassified Archaea" (too few markers) -> lineage d__<domain>
  - "Unclassified" -> status nomatch, empty lineage (routed as unresolved)
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

VERSION = "1.0.0"
RANKS = {"d": "superkingdom", "p": "phylum", "c": "class", "o": "order", "f": "family", "g": "genus", "s": "species"}
EXTRA = ("classification_method", "closest_genome_reference", "closest_genome_ani", "note", "warnings")
HEADER = ("query_name", "status", "rank", "lineage") + EXTRA


def convert_row(row: dict) -> dict:
    classification = (row.get("classification") or "").strip()
    if classification.startswith("Unclassified "):
        lineage = "d__" + classification.split(" ", 1)[1].strip()
    elif classification.startswith("d__"):
        lineage = ";".join(token.strip() for token in classification.split(";") if not token.strip().endswith("__"))
    else:
        lineage = ""
    rank = RANKS[lineage.split(";")[-1][0]] if lineage else ""
    out = {"query_name": row["user_genome"].strip(), "status": "match" if lineage else "nomatch", "rank": rank, "lineage": lineage}
    out.update({key: (row.get(key) or "").strip() for key in EXTRA})
    return out


def convert(paths: list[Path]) -> list[dict]:
    rows = []
    for path in paths:
        with path.open(newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if reader.fieldnames is None:
                continue
            if "user_genome" not in reader.fieldnames or "classification" not in reader.fieldnames:
                raise ValueError(f"{path}: not a GTDB-Tk summary (needs user_genome and classification)")
            rows.extend(convert_row(row) for row in reader)
    return sorted(rows, key=lambda row: row["query_name"])


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--summary", type=Path, nargs="+", required=True, help="GTDB-Tk *.bac120.summary.tsv / *.ar53.summary.tsv")
    parser.add_argument("--output", type=Path, required=True, help="CSV in the sourmash tax genome layout")
    parser.add_argument("--version", action="version", version=VERSION)
    args = parser.parse_args(argv)
    try:
        rows = convert(args.summary)
    except (OSError, csv.Error, ValueError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 2
    with args.output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=HEADER, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
