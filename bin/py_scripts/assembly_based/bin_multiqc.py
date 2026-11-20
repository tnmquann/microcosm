#!/usr/bin/env python3
"""MultiQC files for the assembly-based bins (MultiQC has no module for sourmash / GTDB-Tk bin taxonomy).

Writes, from whatever inputs are given:
  bin_taxonomy_mqc.json  bin classification (sourmash tax genome or GTDB-Tk, --classifier): one bar graph with a tab per rank
                         (Species ... Domain), one bar per bin set (<assembler>-<binner>-<sample>), value =
                         number of bins; top --top taxa per rank, rest "Other", unclassified separate
                         (same layout as the read-based microcosm section).
  bin_summary_mqc.json   one table row per bin: CheckM2 / CheckM / BUSCO completeness and contamination,
                         QUAST size, contigs, N50, GC and the GTDB taxon (like nf-core/mag's bin summary).
  report_bins.tsv        QUAST bin summary transposed back to QUAST's report.tsv layout for MultiQC's quast module.
"""

from __future__ import annotations

import argparse
import csv
import json
import re
import sys
from collections import defaultdict
from pathlib import Path

VERSION = "1.1.0"
RANKS = [("Species", 6), ("Genus", 5), ("Family", 4), ("Order", 3), ("Class", 2), ("Phylum", 1), ("Domain", 0)]
OTHER, UNCLASSIFIED = "zzy_other", "zzz_unclassified"
EXT = re.compile(r"(\.(fa|fasta|fna|fas))?(\.gz)?$", re.IGNORECASE)


# --classifier: section label, method, results folder
CLASSIFIERS = {
    "sourmash": ("sourmash", "sourmash tax genome", "sourmash"),
    "gtdbtk": ("GTDB-Tk", "GTDB-Tk classify_wf", "GTDB-Tk"),
}


def bin_name(name: str) -> str:
    return EXT.sub("", Path(name.strip()).name)


def bin_set(name: str) -> str:
    """Bin set label: bin name without its trailing '.<bin id>'."""
    return name.rsplit(".", 1)[0] if "." in name else name


def read_table(path: Path, delimiter: str) -> list[dict]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter=delimiter))


def number(value: str | None) -> float | None:
    try:
        return float(str(value).rstrip("%"))
    except (TypeError, ValueError):
        return None


def lineage_taxa(row: dict) -> list[str]:
    """Lineage parts of a sourmash 'match'; [] when the bin is not classified."""
    if (row.get("status") or "").strip() != "match":
        return []
    return [part.strip() for part in (row.get("lineage") or "").split(";") if part.strip() and not part.strip().endswith("__")]


def taxonomy_plot(rows: list[dict], top: int, classifier: str = "sourmash") -> dict | None:
    datasets, categories, labels = [], [], []
    for label, index in RANKS:
        counts: dict[str, dict[str, float]] = defaultdict(lambda: defaultdict(float))
        for row in rows:
            taxa = lineage_taxa(row)
            counts[bin_set(row["query_name"])][taxa[index] if index < len(taxa) else UNCLASSIFIED] += 1
        totals: dict[str, float] = defaultdict(float)
        for per_set in counts.values():
            for name, value in per_set.items():
                if name != UNCLASSIFIED:
                    totals[name] += value
        kept = [name for name, _ in sorted(totals.items(), key=lambda item: (-item[1], item[0]))[:top]]
        data = {}
        for set_name, per_set in sorted(counts.items()):
            out = {name: per_set.get(name, 0) for name in kept}
            out[OTHER] = sum(v for n, v in per_set.items() if n not in kept and n != UNCLASSIFIED)
            out[UNCLASSIFIED] = per_set.get(UNCLASSIFIED, 0)
            data[set_name] = out
        cats = {name: {"name": name} for name in kept}
        cats[OTHER] = {"name": "Other", "color": "#cccccc"}
        cats[UNCLASSIFIED] = {"name": "Unclassified", "color": "#8c8c8c"}
        datasets.append(data)
        categories.append(cats)
        labels.append(label)
    if not rows:
        return None
    return {
        "id": "microcosm_bin_taxonomy",
        "section_name": f"Bin taxonomy ({CLASSIFIERS[classifier][0]})",
        "section_anchor": "microcosm_bin_taxonomy",
        "description": (
            f"{CLASSIFIERS[classifier][1]} classification of the bins against GTDB, counted per bin set "
            f"(&lt;assembler&gt;-&lt;binner&gt;-&lt;sample&gt;). Top {top} taxa per rank; Other = remaining classified bins, "
            f"Unclassified = bins without a match at that rank. Full tables: Taxonomy/{CLASSIFIERS[classifier][2]}/."
        ),
        "plot_type": "bargraph",
        "pconfig": {
            "id": "microcosm-bin-taxonomy-plot",
            "title": f"Bins: {CLASSIFIERS[classifier][0]} taxonomy",
            "ylab": "Number of bins",
            "data_labels": labels,
            "cpswitch_c_active": True,
        },
        "categories": categories,
        "data": datasets,
    }


HEADERS = {
    "checkm2_completeness": {"title": "CheckM2 compl.", "description": "CheckM2 completeness", "suffix": "%", "max": 100, "min": 0, "scale": "RdYlGn"},
    "checkm2_contamination": {"title": "CheckM2 contam.", "description": "CheckM2 contamination", "suffix": "%", "min": 0, "scale": "RdYlGn-rev"},
    "checkm_completeness": {"title": "CheckM compl.", "description": "CheckM completeness", "suffix": "%", "max": 100, "min": 0, "scale": "RdYlGn"},
    "checkm_contamination": {"title": "CheckM contam.", "description": "CheckM contamination", "suffix": "%", "min": 0, "scale": "RdYlGn-rev"},
    "busco_complete": {"title": "BUSCO complete", "description": "BUSCO complete (single + duplicated)", "suffix": "%", "max": 100, "min": 0, "scale": "RdYlGn"},
    "busco_dataset": {"title": "BUSCO lineage", "description": "BUSCO lineage dataset"},
    "size": {"title": "Size", "description": "Total length (QUAST)", "suffix": " bp", "format": "{:,.0f}", "scale": "Blues"},
    "contigs": {"title": "Contigs", "description": "Number of contigs (QUAST)", "format": "{:,.0f}", "scale": "Oranges"},
    "n50": {"title": "N50", "description": "N50 (QUAST)", "suffix": " bp", "format": "{:,.0f}", "scale": "Greens"},
    "gc": {"title": "GC", "description": "GC content (QUAST)", "suffix": "%", "format": "{:,.1f}"},
    "gtdb_ani": {"title": "GTDB-Tk ANI", "description": "ANI to the closest GTDB reference genome (GTDB-Tk)", "suffix": "%", "max": 100, "min": 75, "scale": "RdYlGn"},
    "gtdb_method": {"title": "GTDB-Tk method", "description": "GTDB-Tk classification method (ANI screen, topology + ANI, RED ...)"},
    "taxon": {"title": "GTDB taxon", "description": "Lowest GTDB rank assigned by the bin classifier (sourmash tax genome or GTDB-Tk)"},
}


def summary_rows(taxonomy: list[dict], qc: list[Path], quast: list[dict]) -> dict[str, dict]:
    rows: dict[str, dict] = defaultdict(dict)
    for path in qc:
        table = read_table(path, "\t")
        header = set(table[0]) if table else set()
        for row in table:
            if "Completeness_Model_Used" in header or ("Name" in header and "Completeness" in header):
                rows[bin_name(row["Name"])].update(checkm2_completeness=number(row["Completeness"]), checkm2_contamination=number(row["Contamination"]))
            elif "Bin Id" in header:
                rows[bin_name(row["Bin Id"])].update(checkm_completeness=number(row["Completeness"]), checkm_contamination=number(row["Contamination"]))
            elif "Input_file" in header:
                rows[bin_name(row["Input_file"])].update(busco_complete=number(row.get("Complete")), busco_dataset=row.get("Dataset"))
    for row in quast:
        rows[bin_name(row["Assembly"])].update(
            size=number(row.get("Total length (>= 0 bp)") or row.get("Total length")),
            contigs=number(row.get("# contigs (>= 0 bp)") or row.get("# contigs")),
            n50=number(row.get("N50")),
            gc=number(row.get("GC (%)")),
        )
    for row in taxonomy:
        taxa = lineage_taxa(row)
        rows[bin_name(row["query_name"])]["taxon"] = taxa[-1] if taxa else "unclassified"
        # GTDB-Tk (gtdbtk_to_taxgenome.py columns): ANI to the closest GTDB genome and how the taxon was assigned
        if row.get("classification_method"):
            rows[bin_name(row["query_name"])].update(gtdb_ani=number(row.get("closest_genome_ani")), gtdb_method=row["classification_method"])
    return {name: {k: v for k, v in values.items() if v is not None} for name, values in sorted(rows.items())}


def summary_table(rows: dict[str, dict]) -> dict | None:
    if not rows:
        return None
    used = [key for key in HEADERS if any(key in row for row in rows.values())]
    return {
        "id": "microcosm_bin_summary",
        "section_name": "Bin summary",
        "section_anchor": "microcosm_bin_summary",
        "description": "One row per bin: quality (CheckM2 / CheckM / BUSCO), assembly statistics (QUAST) and the GTDB taxon of the bin classifier.",
        "plot_type": "table",
        "pconfig": {"id": "microcosm-bin-summary-table", "title": "Bins: summary", "namespace": "Bins", "col1_header": "Bin"},
        "headers": {key: HEADERS[key] for key in used},
        "data": rows,
    }


def write_quast_report(quast: list[dict], path: Path) -> None:
    """QUAST transposed_report layout (one row per assembly) -> report.tsv layout (one row per metric)."""
    metrics = [key for key in quast[0] if key != "Assembly"]
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["Assembly"] + [row["Assembly"] for row in quast])
        for metric in metrics:
            writer.writerow([metric] + [row.get(metric, "") for row in quast])


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--taxgenome", type=Path, nargs="*", default=[], help="combined_genome_classification CSVs (sourmash tax genome layout; GTDB-Tk converted to it)")
    parser.add_argument("--qc", type=Path, nargs="*", default=[], help="<tool>_summary.tsv from bin QC (checkm2, checkm, busco)")
    parser.add_argument("--quast", type=Path, nargs="*", default=[], help="QUAST bin summary (quast_summary.tsv)")
    parser.add_argument("--outdir", type=Path, default=Path("."))
    parser.add_argument("--top", type=int, default=10)
    parser.add_argument("--classifier", choices=sorted(CLASSIFIERS), default="sourmash", help="bin classifier that wrote --taxgenome")
    parser.add_argument("--version", action="version", version=VERSION)
    args = parser.parse_args(argv)

    taxonomy = [row for path in args.taxgenome for row in read_table(path, ",")]
    quast = [row for path in args.quast for row in read_table(path, "\t")]
    args.outdir.mkdir(parents=True, exist_ok=True)
    written = []
    for name, content in [("bin_taxonomy_mqc.json", taxonomy_plot(taxonomy, args.top, args.classifier)), ("bin_summary_mqc.json", summary_table(summary_rows(taxonomy, args.qc, quast)))]:
        if content:
            (args.outdir / name).write_text(json.dumps(content, indent=1) + "\n")
            written.append(name)
    if quast:
        write_quast_report(quast, args.outdir / "report_bins.tsv")
        written.append("report_bins.tsv")
    print(f"{len(taxonomy)} classified bins, {len(quast)} QUAST rows -> {written}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
