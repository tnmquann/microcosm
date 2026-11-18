#!/usr/bin/env python3
"""MultiQC custom content for the microcosm read-based profile.

Reads <prefix>.processed_metaflow.csv (post_processing.py default output) and
writes microcosm_mqc.json: one "microcosm" section with one "Top taxa" bar
graph that has a tab per rank (Strain = genome, Species ... Domain), laid out
like MultiQC's sylph-tax module. Values are percent of weighted k-mers
(f_weighted_at_rank); per rank the --top taxa by summed abundance are shown,
the rest of the classified part is "Other", and "Unclassified" stays separate.
When the table holds several YACHT min_coverage levels, each level is its own bar. A genome's
min_coverage cell lists every level it passes ("min_coverage0.1 | min_coverage0.05 | ..."), so it
joins each of those bars; Unclassified of a level bar is 100 - its classified share (the
untagged unclassified rows of the table cannot be assigned to a level).
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from collections import defaultdict
from pathlib import Path

VERSION = "2.1.0"
# tab label -> lineage index (None = genome / match_name)
RANKS = [("Strain", None), ("Species", 6), ("Genus", 5), ("Family", 4), ("Order", 3), ("Class", 2), ("Phylum", 1), ("Domain", 0)]
OTHER, UNCLASSIFIED = "zzy_other", "zzz_unclassified"


def row_levels(row: dict) -> list[str]:
    return [level.strip() for level in (row.get("min_coverage") or "").split("|") if level.strip()]


def expand_levels(rows: list[dict]) -> tuple[list[tuple[str, dict]], bool]:
    """(bar label, row) pairs and whether bars are per level (then untagged unclassified rows are dropped)."""
    if len({level for row in rows for level in row_levels(row)}) <= 1:
        return [(row["sample_name"], row) for row in rows], False
    return [(f"{row['sample_name']} ({level})", row) for row in rows for level in row_levels(row)], True


def taxon(row: dict, index: int | None) -> str | None:
    """Taxon of a row at a rank; None = unclassified at that rank."""
    lineage = (row.get("lineage") or "").strip()
    if not lineage or lineage == "unclassified":
        return None
    if index is None:
        return row.get("match_name") or None
    parts = lineage.split(";")
    if index >= len(parts):
        return None
    name = parts[index].strip()
    return name if name and not name.endswith("__") else None


def build(rows: list[dict], top: int) -> tuple[list[dict], list[dict], list[str]]:
    pairs, per_level = expand_levels(rows)
    remainder: dict[str, float] = defaultdict(lambda: 100.0)
    if per_level:
        for bar, row in pairs:
            remainder[bar] -= float(row.get("f_weighted_at_rank") or 0) * 100
    datasets, categories, labels = [], [], []
    for label, index in RANKS:
        values: dict[str, dict[str, float]] = defaultdict(lambda: defaultdict(float))
        for bar, row in pairs:
            name = taxon(row, index)
            values[bar][name or UNCLASSIFIED] += float(row.get("f_weighted_at_rank") or 0) * 100
        for bar in remainder:
            values[bar][UNCLASSIFIED] += max(remainder[bar], 0.0)
        totals: dict[str, float] = defaultdict(float)
        for sample in values.values():
            for name, value in sample.items():
                if name != UNCLASSIFIED:
                    totals[name] += value
        if not totals:
            continue
        kept = [name for name, _ in sorted(totals.items(), key=lambda item: (-item[1], item[0]))[:top]]
        data = {}
        for bar, sample in values.items():
            out = {name: round(sample.get(name, 0.0), 4) for name in kept}
            out[OTHER] = round(sum(v for n, v in sample.items() if n not in kept and n != UNCLASSIFIED), 4)
            out[UNCLASSIFIED] = round(sample.get(UNCLASSIFIED, 0.0), 4)
            data[bar] = out
        cats = {name: {"name": name} for name in kept}
        cats[OTHER] = {"name": "Other", "color": "#cccccc"}
        cats[UNCLASSIFIED] = {"name": "Unclassified", "color": "#8c8c8c"}
        datasets.append(data)
        categories.append(cats)
        labels.append(label)
    return datasets, categories, labels


def write_json(path: Path, datasets: list[dict], categories: list[dict], labels: list[str], top: int) -> None:
    content = {
        "id": "microcosm",
        "section_name": "microcosm",
        "section_anchor": "microcosm",
        "description": (
            "Read-based taxonomic profile from microcosm (sourmash fastmultigather + YACHT, GTDB). "
            f"Top {top} taxa per rank as percent of weighted k-mers; Other = remaining classified taxa, "
            "Unclassified = sequence not assigned at that rank. Full table: "
            "Readbased_Analysis/Post_processing/&lt;prefix&gt;_default/&lt;prefix&gt;.processed_metaflow.csv."
        ),
        "plot_type": "bargraph",
        "pconfig": {
            "id": "microcosm-top-n-plot",
            "title": "microcosm: Top taxa",
            "ylab": "Relative Abundance (% weighted k-mers)",
            "data_labels": labels,
            "cpswitch": False,
        },
        "categories": categories,
        "data": datasets,
    }
    path.write_text(json.dumps(content, indent=1) + "\n")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--input", type=Path, required=True, help="<prefix>.processed_metaflow.csv")
    parser.add_argument("--outdir", type=Path, default=Path("."))
    parser.add_argument("--top", type=int, default=10)
    parser.add_argument("--version", action="version", version=VERSION)
    args = parser.parse_args(argv)
    with args.input.open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    missing = {"sample_name", "lineage", "f_weighted_at_rank", "match_name"} - set(rows[0] if rows else {})
    if rows and missing:
        raise SystemExit(f"{args.input}: missing column(s) {sorted(missing)}")
    args.outdir.mkdir(parents=True, exist_ok=True)
    datasets, categories, labels = build(rows, args.top)
    write_json(args.outdir / "microcosm_mqc.json", datasets, categories, labels, args.top)
    print(f"{len(rows)} rows, ranks {labels} -> {args.outdir / 'microcosm_mqc.json'}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
