#!/usr/bin/env python3
"""MultiQC custom-content tables from the microcosm processed read-based table.

Reads <prefix>.processed_metaflow.csv (post_processing.py default output) and
writes microcosm_<rank>_mqc.tsv bar graphs (percent of weighted k-mers per
sample, f_weighted_at_rank) for phylum, genus and species, in the spirit of
MultiQC's sylph module: the --top taxa by mean abundance are kept, the rest
are summed into "Other", "unclassified" stays its own category. When the table
holds several YACHT min_coverage levels, each level is a separate sample row.
"""

from __future__ import annotations

import argparse
import csv
import sys
from collections import defaultdict
from pathlib import Path

VERSION = "1.0.0"
RANKS = {"phylum": 1, "genus": 5, "species": 6}


def expand_levels(rows: list[dict]) -> list[tuple[str, dict]]:
    """(sample label, row) pairs; rows without min_coverage (synthetic unclassified) join every level."""
    levels = sorted({row.get("min_coverage", "") for row in rows} - {""})
    if len(levels) <= 1:
        return [(row["sample_name"], row) for row in rows]
    pairs = []
    for row in rows:
        for level in [row["min_coverage"]] if row.get("min_coverage") else levels:
            pairs.append((f"{row['sample_name']} ({level})", row))
    return pairs


def taxon_at_rank(lineage: str, index: int) -> str:
    if not lineage or lineage == "unclassified":
        return "unclassified"
    parts = lineage.split(";")
    if index >= len(parts) or not parts[index].strip() or parts[index].strip().endswith("__"):
        return "unassigned"
    return parts[index].strip()


def build_tables(rows: list[dict], top: int) -> dict[str, dict[str, dict[str, float]]]:
    pairs = expand_levels(rows)
    tables = {}
    for rank, index in RANKS.items():
        values: dict[str, dict[str, float]] = defaultdict(lambda: defaultdict(float))
        for label, row in pairs:
            value = float(row.get("f_weighted_at_rank") or 0) * 100
            values[label][taxon_at_rank(row.get("lineage", ""), index)] += value
        taxa = {taxon for sample in values.values() for taxon in sample} - {"unclassified"}
        ranked = sorted(taxa, key=lambda taxon: (-sum(s.get(taxon, 0) for s in values.values()), taxon))
        kept = ranked[:top]
        table = {}
        for sample, sample_values in values.items():
            out = {taxon: sample_values.get(taxon, 0.0) for taxon in kept}
            out["Other"] = sum(v for t, v in sample_values.items() if t not in kept and t != "unclassified")
            out["unclassified"] = sample_values.get("unclassified", 0.0)
            table[sample] = out
        tables[rank] = table
    return tables


def write_table(path: Path, rank: str, table: dict[str, dict[str, float]]) -> None:
    columns = list(next(iter(table.values())).keys()) if table else []
    header = [
        f"# id: 'microcosm_{rank}'",
        f"# section_name: 'microcosm: {rank}'",
        f"# description: 'microcosm (sourmash + YACHT) read-based profile at {rank} level: percent of weighted k-mers (f_weighted_at_rank), top taxa plus Other and unclassified.'",
        "# plot_type: 'bargraph'",
        "# pconfig:",
        f"#     id: 'microcosm_{rank}_plot'",
        f"#     title: 'microcosm: {rank}'",
        "#     ylab: '% weighted k-mers'",
        "#     cpswitch_c_active: False",
    ]
    with path.open("w", newline="") as handle:
        handle.write("\n".join(header) + "\n")
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["Sample"] + columns)
        for sample in sorted(table):
            writer.writerow([sample] + [f"{table[sample][c]:.4f}" for c in columns])


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--input", type=Path, required=True, help="<prefix>.processed_metaflow.csv")
    parser.add_argument("--outdir", type=Path, default=Path("."))
    parser.add_argument("--top", type=int, default=10)
    parser.add_argument("--version", action="version", version=VERSION)
    args = parser.parse_args(argv)
    with args.input.open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    missing = {"sample_name", "lineage", "f_weighted_at_rank"} - set(rows[0] if rows else {})
    if rows and missing:
        raise SystemExit(f"{args.input}: missing column(s) {sorted(missing)}")
    args.outdir.mkdir(parents=True, exist_ok=True)
    for rank, table in build_tables(rows, args.top).items():
        write_table(args.outdir / f"microcosm_{rank}_mqc.tsv", rank, table)
    print(f"{len(rows)} rows -> {args.outdir}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
