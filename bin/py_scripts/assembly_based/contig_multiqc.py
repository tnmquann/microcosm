#!/usr/bin/env python3
"""MultiQC files for the contig analysis of each assembly (no MultiQC module reads these tables).

Writes, from whatever inputs are given (one row / bar per assembly, <assembler>-<sample>):
  contig_summary_mqc.json  table: contigs and length (SeqKit fx2tab), mapped reads, mean depth and breadth
                           (samtools coverage), predicted genes (Pyrodigal / Prokka .faa), geNomad plasmid and
                           virus contigs (+ plasmid AMR genes, proviruses), contigs with a skani hit
  contig_tiara_mqc.json    Tiara domain classification: bar graph, tabs "Contigs" and "Length (bp)"
  contig_skani_mqc.json    skani search (--qi): best reference per contig (highest ANI, then aligned fraction);
                           bar graph of the top references, tabs "Length (bp)" and "Contigs"
Columns of a skipped tool are left out.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import re
import sys
from collections import defaultdict
from pathlib import Path

VERSION = "1.1.0"
# file name -> assembly name: MEGAHIT-s1.tiara.txt, MEGAHIT-s1.contigs_plasmid_summary.tsv, MEGAHIT-s1.pyrodigal.faa ...
SUFFIX = re.compile(r"(_plasmid_summary|_virus_summary)?\.(txt|tsv|faa)(\.gz)?$")
TOOL = re.compile(r"\.(tiara|seqkit|skani|coverage|pyrodigal|prokka|geNomad|contigs)$")
TIARA_CLASSES = ["bacteria", "archaea", "prokarya", "eukarya", "organelle", "unknown"]
SKANI_OTHER, SKANI_NONE = "Other references", "No skani hit"


def assembly_name(path: Path) -> str:
    name = SUFFIX.sub("", path.name)
    while TOOL.search(name):
        name = TOOL.sub("", name)
    return name


def open_text(path: Path):
    return gzip.open(path, "rt", newline="") if path.name.endswith(".gz") else path.open(newline="")


def read_tsv(path: Path) -> list[dict]:
    with open_text(path) as handle:
        lines = [line for line in handle if line.strip()]
    if lines:
        lines[0] = lines[0].lstrip("#")  # samtools coverage: #rname
    return list(csv.DictReader(lines, delimiter="\t"))


def number(value) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return 0.0


def summary(args) -> tuple[dict, dict]:
    rows: dict[str, dict] = defaultdict(dict)
    lengths: dict[str, dict[str, float]] = defaultdict(dict)
    for path in args.seqkit:
        table = read_tsv(path)
        name = assembly_name(path)
        id_col = next((key for key in ("id", "name") if table and key in table[0]), None)
        for row in table:
            if id_col:
                lengths[name][row[id_col].split()[0]] = number(row.get("length"))
        rows[name]["contigs"] = len(table)
        rows[name]["length"] = sum(number(row.get("length")) for row in table)
    for path in args.coverage:
        table = read_tsv(path)
        size = sum(number(row["endpos"]) - number(row["startpos"]) + 1 for row in table)
        row = rows[assembly_name(path)]
        row["mapped_reads"] = sum(number(r["numreads"]) for r in table)
        row["mean_depth"] = sum(number(r["meandepth"]) * (number(r["endpos"]) - number(r["startpos"]) + 1) for r in table) / size if size else 0.0
        row["breadth"] = 100 * sum(number(r["covbases"]) for r in table) / size if size else 0.0
    for path in args.faa:
        with open_text(path) as handle:
            rows[assembly_name(path)]["genes"] = sum(1 for line in handle if line.startswith(">"))
    for path in args.plasmid:
        table = read_tsv(path)
        row = rows[assembly_name(path)]
        row["plasmids"] = len(table)
        row["plasmid_bp"] = sum(number(r.get("length")) for r in table)
        row["plasmid_amr"] = sum(len([g for g in (r.get("amr_genes") or "").split(";") if g and g != "NA"]) for r in table)
    for path in args.virus:
        table = read_tsv(path)
        row = rows[assembly_name(path)]
        row["viruses"] = len(table)
        row["proviruses"] = sum(1 for r in table if (r.get("topology") or "").lower() == "provirus")
        row["virus_bp"] = sum(number(r.get("length")) for r in table)
    for name, hits in skani_best_hits(args.skani).items():
        row = rows[name]
        row["skani_contigs"] = len(hits)
        total = sum(lengths.get(name, {}).values())
        if total:
            row["skani_bp"] = 100 * sum(lengths[name].get(contig, 0.0) for contig in hits) / total
    return rows, lengths


def reference_label(row: dict) -> str:
    """'NC_016845.1 Klebsiella pneumoniae subsp. pneumoniae HS11286 chromosome, complete genome' -> species/strain text."""
    words = (row.get("Ref_name") or "").split()
    label = " ".join(words[1:]).split(",")[0].strip()
    label = re.sub(r"\s+(chromosome|plasmid|complete genome|genome assembly).*$", "", label)
    return label or SUFFIX.sub("", Path(row.get("Ref_file") or "unknown").name)


def skani_best_hits(paths: list[Path]) -> dict[str, dict[str, str]]:
    """assembly -> {contig: reference label} for the best hit of each contig."""
    best: dict[str, dict[str, tuple]] = defaultdict(dict)
    for path in paths:
        name = assembly_name(path)
        for row in read_tsv(path):
            contig = (row.get("Query_name") or "").split()[0] if row.get("Query_name") else ""
            if not contig:
                continue
            score = (number(row.get("ANI")), number(row.get("Align_fraction_query")))
            if contig not in best[name] or score > best[name][contig][0]:
                best[name][contig] = (score, reference_label(row))
    for path in paths:  # assemblies without any hit still get a row
        best.setdefault(assembly_name(path), {})
    return {name: {contig: hit[1] for contig, hit in hits.items()} for name, hits in best.items()}


def skani_plot(paths: list[Path], lengths: dict, top: int) -> dict | None:
    hits = skani_best_hits(paths)
    if not hits:
        return None
    sizes: dict[str, dict[str, float]] = {}
    counts: dict[str, dict[str, float]] = {}
    for name, contigs in hits.items():
        size, count = defaultdict(float), defaultdict(float)
        for contig, label in contigs.items():
            size[label] += lengths.get(name, {}).get(contig, 0.0)
            count[label] += 1
        if lengths.get(name):
            size[SKANI_NONE] = sum(length for contig, length in lengths[name].items() if contig not in contigs)
            count[SKANI_NONE] = sum(1 for contig in lengths[name] if contig not in contigs)
        sizes[name], counts[name] = size, count
    has_bp = any(sum(value for label, value in size.items() if label != SKANI_NONE) for size in sizes.values())
    ranking = defaultdict(float)
    for name in hits:
        for label, value in (sizes if has_bp else counts)[name].items():
            if label != SKANI_NONE:
                ranking[label] += value
    keep = sorted(ranking, key=lambda label: -ranking[label])[:top]
    categories = keep + [SKANI_OTHER, SKANI_NONE]

    def collapse(table: dict[str, dict[str, float]]) -> dict[str, dict[str, float]]:
        out = {}
        for name, values in table.items():
            row = {label: values.get(label, 0.0) for label in keep}
            row[SKANI_OTHER] = sum(value for label, value in values.items() if label not in keep and label != SKANI_NONE)
            row[SKANI_NONE] = values.get(SKANI_NONE, 0.0)
            out[name] = row
        return out

    data, labels = ([collapse(sizes), collapse(counts)], ["Length (bp)", "Contigs"]) if has_bp else ([collapse(counts)], ["Contigs"])
    return {
        "id": "microcosm_contig_skani",
        "section_name": "Contig ANI search (skani)",
        "section_anchor": "microcosm_contig_skani",
        "description": f"skani search of each contig against --skani_db; every contig counts for its best reference (highest ANI, "
                       f"then aligned fraction). Top {top} references; No skani hit = contigs without a hit (short contigs "
                       "are not searched). Full tables: ContigsAnalysis/skani/.",
        "plot_type": "bargraph",
        "pconfig": {"id": "microcosm-contig-skani-plot", "title": "Contigs: skani best hits", "ylab": labels[0], "data_labels": labels, "cpswitch_c_active": False},
        "categories": categories,
        "data": data,
    }


HEADERS = {
    "contigs": {"title": "Contigs", "description": "Number of contigs (SeqKit)", "format": "{:,.0f}", "scale": "Oranges"},
    "length": {"title": "Length", "description": "Total contig length (SeqKit)", "suffix": " bp", "format": "{:,.0f}", "scale": "Blues"},
    "mapped_reads": {"title": "Mapped reads", "description": "Reads mapped to the contigs (BWA-MEM2, samtools coverage)", "format": "{:,.0f}", "scale": "Purples"},
    "mean_depth": {"title": "Mean depth", "description": "Length-weighted mean read depth over all contigs", "suffix": "X", "format": "{:,.1f}", "scale": "Greens"},
    "breadth": {"title": "Covered", "description": "Contig bases covered by at least one read", "suffix": "%", "max": 100, "min": 0, "scale": "RdYlGn"},
    "genes": {"title": "Genes", "description": "Predicted protein-coding genes (contig annotation)", "format": "{:,.0f}", "scale": "Greens"},
    "plasmids": {"title": "Plasmid contigs", "description": "geNomad plasmid contigs", "format": "{:,.0f}", "scale": "Reds"},
    "plasmid_bp": {"title": "Plasmid bp", "description": "Total length of geNomad plasmid contigs", "suffix": " bp", "format": "{:,.0f}", "scale": "Reds"},
    "plasmid_amr": {"title": "Plasmid AMR genes", "description": "AMR genes on geNomad plasmid contigs (amr_genes)", "format": "{:,.0f}", "scale": "Reds"},
    "skani_contigs": {"title": "skani hits", "description": "Contigs with a skani hit in --skani_db", "format": "{:,.0f}", "scale": "Greens"},
    "skani_bp": {"title": "skani hit bp", "description": "Share of the assembly length in contigs with a skani hit", "suffix": "%", "max": 100, "min": 0, "scale": "Greens"},
    "viruses": {"title": "Virus contigs", "description": "geNomad virus contigs (including proviruses)", "format": "{:,.0f}", "scale": "Blues"},
    "proviruses": {"title": "Proviruses", "description": "geNomad proviruses (virus regions inside host contigs)", "format": "{:,.0f}", "scale": "Blues"},
    "virus_bp": {"title": "Virus bp", "description": "Total length of geNomad virus sequences", "suffix": " bp", "format": "{:,.0f}", "scale": "Blues"},
}


def summary_table(rows: dict) -> dict | None:
    if not rows:
        return None
    used = [key for key in HEADERS if any(key in row for row in rows.values())]
    return {
        "id": "microcosm_contig_summary",
        "section_name": "Contig analysis",
        "section_anchor": "microcosm_contig_summary",
        "description": "One row per assembly: contig statistics (SeqKit), read mapping (BWA-MEM2 + samtools coverage), "
                       "predicted genes, geNomad plasmids / viruses and skani hits. Full tables: ContigsAnalysis/.",
        "plot_type": "table",
        "pconfig": {"id": "microcosm-contig-summary-table", "title": "Contigs: summary", "namespace": "Contigs", "col1_header": "Assembly"},
        "headers": {key: HEADERS[key] for key in used},
        "data": {name: {key: row[key] for key in used if key in row} for name, row in sorted(rows.items())},
    }


def tiara_plot(paths: list[Path], lengths: dict) -> dict | None:
    counts: dict[str, dict[str, float]] = {}
    sizes: dict[str, dict[str, float]] = {}
    for path in paths:
        name = assembly_name(path)
        count, size = defaultdict(float), defaultdict(float)
        for row in read_tsv(path):
            label = (row.get("class_fst_stage") or "unknown").strip().lower()
            label = label if label in TIARA_CLASSES else "unknown"
            count[label] += 1
            size[label] += lengths.get(name, {}).get((row.get("sequence_id") or "").split()[0], 0.0)
        counts[name] = {label: count[label] for label in TIARA_CLASSES}
        sizes[name] = {label: size[label] for label in TIARA_CLASSES}
    if not counts:
        return None
    data, labels = [counts], ["Contigs"]
    if any(sum(value.values()) for value in sizes.values()):
        data.append(sizes)
        labels.append("Length (bp)")
    return {
        "id": "microcosm_contig_tiara",
        "section_name": "Contig domains (Tiara)",
        "section_anchor": "microcosm_contig_tiara",
        "description": "Tiara classification of the contigs (only contigs of at least --tiara_min_length are classified). "
                       "Full tables: Taxonomy/Tiara/.",
        "plot_type": "bargraph",
        "pconfig": {"id": "microcosm-contig-tiara-plot", "title": "Contigs: Tiara domains", "ylab": "Contigs", "data_labels": labels, "cpswitch_c_active": False},
        "categories": TIARA_CLASSES,
        "data": data,
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    for flag, text in [("--seqkit", "SeqKit fx2tab tables"), ("--coverage", "samtools coverage tables"),
                       ("--faa", "predicted proteins (.faa)"), ("--plasmid", "geNomad *_plasmid_summary.tsv"),
                       ("--virus", "geNomad *_virus_summary.tsv"), ("--tiara", "Tiara classification tables"),
                       ("--skani", "skani search tables (--qi)")]:
        parser.add_argument(flag, type=Path, nargs="*", default=[], help=text)
    parser.add_argument("--top", type=int, default=10, help="references shown in the skani plot")
    parser.add_argument("--outdir", type=Path, default=Path("."))
    parser.add_argument("--version", action="version", version=VERSION)
    args = parser.parse_args(argv)

    rows, lengths = summary(args)
    args.outdir.mkdir(parents=True, exist_ok=True)
    written = []
    for name, content in [("contig_summary_mqc.json", summary_table(rows)), ("contig_tiara_mqc.json", tiara_plot(args.tiara, lengths)),
                          ("contig_skani_mqc.json", skani_plot(args.skani, lengths, args.top))]:
        if content:
            (args.outdir / name).write_text(json.dumps(content, indent=1) + "\n")
            written.append(name)
    print(f"{len(rows)} assemblies -> {written}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
