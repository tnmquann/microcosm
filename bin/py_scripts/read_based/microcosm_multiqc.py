#!/usr/bin/env python3
"""MultiQC custom content for the microcosm read-based profile.

Reads <prefix>.processed_microcosm.csv (post_processing.py default output) and
writes microcosm_mqc.json: one "microcosm" section with one "Top taxa" bar
graph that has a tab per rank (Strain = genome, Species ... Domain), laid out
like MultiQC's sylph-tax module. Values are percent of weighted k-mers
(f_weighted_at_rank); per rank the --top taxa by summed abundance are shown,
the rest of the classified part is "Other", and "Unclassified" stays separate.
One bar per sample at one YACHT min_coverage level (--min_coverage, the level post-processing
was asked for; the pipeline passes the --min_coverage of --readbased_postprocess_options, default 0.05):
a genome's min_coverage cell lists every level it passes ("min_coverage0.1 | min_coverage0.05"),
so it is kept when that list holds the level. Unclassified = 100 - classified share; with
--remove_unclassified it is left out and the classified taxa are rescaled to 100 %.
"""

from __future__ import annotations

import argparse
import csv
import json
import re
import sys
from collections import defaultdict
from pathlib import Path

VERSION = "3.0.0"
# tab label -> lineage index (None = genome / match_name)
RANKS = [("Strain", None), ("Species", 6), ("Genus", 5), ("Family", 4), ("Order", 3), ("Class", 2), ("Phylum", 1), ("Domain", 0)]
OTHER, UNCLASSIFIED = "zzy_other", "zzz_unclassified"


LEVEL = re.compile(r"min_coverage\s*([0-9eE.+-]+)")


def row_levels(row: dict) -> set[float]:
    return {float(value) for value in LEVEL.findall(row.get("min_coverage") or "")}


def select_rows(rows: list[dict], min_coverage: float | None) -> list[dict]:
    """Classified rows at the chosen level (synthetic unclassified rows are recomputed in build)."""
    classified = [row for row in rows if taxon(row, 0) or taxon(row, None)]
    if min_coverage is None or not any(row_levels(row) for row in classified):
        return classified
    return [row for row in classified if any(abs(level - min_coverage) < 1e-9 for level in row_levels(row))]


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


def build(rows: list[dict], top: int, min_coverage: float | None = None, remove_unclassified: bool = False) -> tuple[list[dict], list[dict], list[str]]:
    selected = select_rows(rows, min_coverage)
    share: dict[str, float] = defaultdict(float)
    for row in selected:
        share[row["sample_name"]] += float(row.get("f_weighted_at_rank") or 0) * 100
    samples = sorted({row["sample_name"] for row in rows})
    datasets, categories, labels = [], [], []
    for label, index in RANKS:
        values: dict[str, dict[str, float]] = {sample: defaultdict(float) for sample in samples}
        for row in selected:
            sample = row["sample_name"]
            value = float(row.get("f_weighted_at_rank") or 0) * 100
            if remove_unclassified and share[sample] > 0:
                value *= 100 / share[sample]
            values[sample][taxon(row, index) or UNCLASSIFIED] += value
        if not remove_unclassified:
            for sample in samples:
                values[sample][UNCLASSIFIED] += max(100 - share[sample], 0.0)
        totals: dict[str, float] = defaultdict(float)
        for sample in values.values():
            for name, value in sample.items():
                if name != UNCLASSIFIED:
                    totals[name] += value
        if not totals:
            continue
        kept = [name for name, _ in sorted(totals.items(), key=lambda item: (-item[1], item[0]))[:top]]
        show_unclassified = any(sample.get(UNCLASSIFIED, 0) > 0 for sample in values.values())
        data = {}
        for bar, sample in values.items():
            out = {name: round(sample.get(name, 0.0), 4) for name in kept}
            out[OTHER] = round(sum(v for n, v in sample.items() if n not in kept and n != UNCLASSIFIED), 4)
            if show_unclassified:
                out[UNCLASSIFIED] = round(sample.get(UNCLASSIFIED, 0.0), 4)
            data[bar] = out
        cats = {name: {"name": name} for name in kept}
        cats[OTHER] = {"name": "Other", "color": "#cccccc"}
        if show_unclassified:
            cats[UNCLASSIFIED] = {"name": "Unclassified", "color": "#8c8c8c"}
        datasets.append(data)
        categories.append(cats)
        labels.append(label)
    return datasets, categories, labels


def write_json(path: Path, datasets: list[dict], categories: list[dict], labels: list[str], top: int,
               min_coverage: float | None = None, remove_unclassified: bool = False) -> None:
    level = f"YACHT min_coverage {min_coverage:g}; " if min_coverage is not None else ""
    scale = ("unclassified sequence removed and taxa rescaled to 100 % (--remove_unclassified). "
             if remove_unclassified else "Unclassified = sequence not assigned at that rank. ")
    content = {
        "id": "microcosm",
        "section_name": "microcosm",
        "section_anchor": "microcosm",
        "description": (
            "Read-based taxonomic profile from microcosm (sourmash fastmultigather + YACHT, GTDB). "
            f"{level}Top {top} taxa per rank as percent of weighted k-mers; Other = remaining classified taxa; "
            f"{scale}Full table: "
            "Readbased_Analysis/Post_processing/&lt;prefix&gt;_default/&lt;prefix&gt;.processed_microcosm.csv."
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
    parser.add_argument("--input", type=Path, required=True, help="<prefix>.processed_microcosm.csv")
    parser.add_argument("--outdir", type=Path, default=Path("."))
    parser.add_argument("--top", type=int, default=10)
    parser.add_argument("--min_coverage", type=float, help="YACHT min_coverage level to plot (default: every row as one bar)")
    parser.add_argument("--remove_unclassified", action="store_true", help="drop Unclassified and rescale taxa to 100 %%")
    parser.add_argument("--version", action="version", version=VERSION)
    args = parser.parse_args(argv)
    with args.input.open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    missing = {"sample_name", "lineage", "f_weighted_at_rank", "match_name"} - set(rows[0] if rows else {})
    if rows and missing:
        raise SystemExit(f"{args.input}: missing column(s) {sorted(missing)}")
    args.outdir.mkdir(parents=True, exist_ok=True)
    datasets, categories, labels = build(rows, args.top, args.min_coverage, args.remove_unclassified)
    if args.min_coverage is not None and not datasets:
        print(f"WARNING: no taxa at min_coverage {args.min_coverage:g} in {args.input}", file=sys.stderr)
    write_json(args.outdir / "microcosm_mqc.json", datasets, categories, labels, args.top, args.min_coverage, args.remove_unclassified)
    print(f"{len(rows)} rows, ranks {labels} -> {args.outdir / 'microcosm_mqc.json'}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
