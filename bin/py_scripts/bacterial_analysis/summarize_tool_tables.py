#!/usr/bin/env python3
"""Concatenate one tool's per-genome tables into one summary TSV.

The manifest is a two-column TSV (genome_id, table path, no header).
Every output row starts with genome_id; columns are the union of the input
headers in first-seen order; rows are sorted by genome_id. A leading '#' on
the first header cell (ABRicate '#FILE', DeepARG '#ARG') is dropped.
--no-header handles header-less tables (mlst): the first --columns names are
used as given and any remaining cells are joined with ' ' into the last name.
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

VERSION = "1.0.0"


def read_table(path: Path, no_header: bool, columns: list[str]) -> tuple[list[str], list[dict]]:
    with path.open(newline="") as handle:
        lines = [line for line in csv.reader(handle, delimiter="\t") if any(cell.strip() for cell in line)]
    if not lines:
        return (columns if no_header else []), []
    if no_header:
        fixed = len(columns) - 1
        rows = [
            dict(zip(columns, cells[:fixed] + [" ".join(cells[fixed:])]))
            for cells in lines
        ]
        return columns, rows
    header = [cell.strip() for cell in lines[0]]
    header[0] = header[0].lstrip("#").strip()
    rows = [dict(zip(header, cells + [""] * (len(header) - len(cells)))) for cells in lines[1:]]
    return header, rows


def summarize(manifest: Path, output: Path, no_header: bool, columns: list[str]) -> int:
    if no_header and len(columns) < 2:
        raise SystemExit("--no-header needs at least two --columns")
    entries = []
    with manifest.open(newline="") as handle:
        for line in csv.reader(handle, delimiter="\t"):
            if line:
                entries.append((line[0], Path(line[1])))
    header: list[str] = []
    out_rows: list[dict] = []
    for genome_id, path in sorted(entries, key=lambda entry: entry[0]):
        table_header, rows = read_table(path, no_header, columns)
        header.extend(column for column in table_header if column not in header and column != "genome_id")
        out_rows.extend({**row, "genome_id": genome_id} for row in rows)
    with output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["genome_id"] + header, delimiter="\t",
                                extrasaction="ignore", lineterminator="\n")
        writer.writeheader()
        writer.writerows(out_rows)
    return len(out_rows)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--no-header", action="store_true")
    parser.add_argument("--columns", default="", help="comma-separated names for --no-header tables")
    parser.add_argument("--version", action="version", version=VERSION)
    args = parser.parse_args(argv)
    columns = [name for name in args.columns.split(",") if name]
    count = summarize(args.manifest, args.output, args.no_header, columns)
    print(f"{args.output}: {count} rows", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
