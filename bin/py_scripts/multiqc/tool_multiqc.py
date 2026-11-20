#!/usr/bin/env python3
"""MultiQC custom-content files for tools that MultiQC has no module for.

One call per tool: --tool <name> <files or directories> writes <name>_mqc.json (and, for some tools, a second
plot file). Directories are searched recursively. A tool with no usable rows writes nothing.

  contig_summary  SeqKit fx2tab + samtools coverage: contigs, length, mapped reads, depth, breadth (per assembly)
  pyrodigal       Pyrodigal .faa: genes, complete / partial genes, mean protein length (per assembly)
  genomad         geNomad *_plasmid_summary.tsv / *_virus_summary.tsv (per assembly)
  tiara_contigs   Tiara contig classification (+ SeqKit tables for lengths): bar graph per assembly
  skani           skani search (--qi) best reference per contig (+ SeqKit tables for lengths): bar graph per assembly
  tiara_bins      Tiara bin classification (<prefix>.binclassification.tsv): bins per domain, per bin set
  gunc            GUNC *maxCSS_level.tsv: per-bin chimerism table
  assembly_stats  assembly_stats.py *.assembly_stats.csv: per-bin statistics
  rgi_bwt         RGI bwt (reads) <sample>.gene_mapping_data.txt (+ overall_mapping_stats.txt): table + AMR gene families
  rgi_amr         RGI main (bacterial genomes) <genome>.rgi.txt: table + AMR gene families
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import os
import re
import sys
from collections import defaultdict
from pathlib import Path

VERSION = "1.2.0"
# file name -> assembly name: MEGAHIT-s1.tiara.txt, MEGAHIT-s1.contigs_plasmid_summary.tsv, MEGAHIT-s1.pyrodigal.faa ...
SUFFIX = re.compile(r"(_plasmid_summary|_virus_summary)?\.(txt|tsv|faa|csv)(\.gz)?$")
TOOL = re.compile(r"\.(tiara|seqkit|skani|coverage|pyrodigal|prokka|geNomad|contigs|gene_mapping_data|overall_mapping_stats|rgi|binclassification|assembly_stats|kleborate|seqsero2_result|genotype_predictResults|salmonella_routing)$")
TIARA_CLASSES = ["bacteria", "archaea", "prokarya", "eukarya", "organelle", "unknown"]
OTHER, SKANI_NONE = "Other", "No skani hit"
csv.field_size_limit(sys.maxsize)  # RGI rows carry whole sequences


def sample_name(path: Path) -> str:
    name = SUFFIX.sub("", path.name)
    while TOOL.search(name):
        name = TOOL.sub("", name)
    return name


def open_text(path: Path):
    return gzip.open(path, "rt", newline="") if path.name.endswith(".gz") else path.open(newline="")


def read_table(path: Path, delimiter: str = "\t") -> list[dict]:
    with open_text(path) as handle:
        lines = [line for line in handle if line.strip()]
    if lines:
        lines[0] = lines[0].lstrip("#")  # samtools coverage: #rname
    return list(csv.DictReader(lines, delimiter=delimiter))


def number(value) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return 0.0


def pick(paths: list[Path], *patterns: str) -> list[Path]:
    return [path for path in paths if any(re.search(pattern, path.name) for pattern in patterns)]


def seqkit_lengths(paths: list[Path]) -> dict[str, dict[str, float]]:
    lengths: dict[str, dict[str, float]] = defaultdict(dict)
    for path in pick(paths, r"\.seqkit\.txt$"):
        for row in read_table(path):
            contig = (row.get("id") or row.get("name") or "").split()
            if contig:
                lengths[sample_name(path)][contig[0]] = number(row.get("length"))
    return lengths


def table(anchor: str, name: str, description: str, headers: dict, rows: dict, col1: str) -> dict | None:
    rows = {sample: row for sample, row in rows.items() if row}
    if not rows:
        return None
    used = [key for key in headers if any(key in row for row in rows.values())]
    return {
        "id": anchor, "section_name": name, "section_anchor": anchor, "description": description, "plot_type": "table",
        "pconfig": {"id": f"{anchor.replace('_', '-')}-table", "title": name, "namespace": name, "col1_header": col1},
        "headers": {key: headers[key] for key in used},
        "data": {sample: {key: row[key] for key in used if key in row} for sample, row in sorted(rows.items())},
    }


def bargraph(anchor: str, name: str, description: str, categories: list[str], data: list[dict], labels: list[str]) -> dict | None:
    if not data or not any(any(values.values()) for values in data[0].values()):
        return None
    return {
        "id": anchor, "section_name": name, "section_anchor": anchor, "description": description, "plot_type": "bargraph",
        "pconfig": {"id": f"{anchor.replace('_', '-')}-plot", "title": name, "ylab": labels[0], "data_labels": labels, "cpswitch_c_active": False},
        "categories": categories, "data": data,
    }


def top_categories(values: dict[str, dict[str, float]], top: int, other: str = OTHER, fixed: tuple = ()) -> tuple[list[str], dict]:
    """Keep the top labels over all samples, sum the rest into 'Other'."""
    ranking = defaultdict(float)
    for row in values.values():
        for label, value in row.items():
            if label not in fixed:
                ranking[label] += value
    keep = sorted(ranking, key=lambda label: (-ranking[label], label))[:top]
    out = {}
    for sample, row in values.items():
        new = {label: row.get(label, 0.0) for label in keep}
        new[other] = sum(value for label, value in row.items() if label not in keep and label not in fixed)
        for label in fixed:
            new[label] = row.get(label, 0.0)
        out[sample] = new
    return keep + [other, *fixed], out


def fmt(title, description, scale, suffix="", fmt_="{:,.0f}", **extra):
    header = {"title": title, "description": description, "scale": scale, "format": fmt_}
    if suffix:
        header["suffix"] = suffix
    header.update(extra)
    return header


# ---------- contigs ----------

def contig_summary(paths: list[Path], _top: int) -> list:
    rows: dict[str, dict] = defaultdict(dict)
    for path in pick(paths, r"\.seqkit\.txt$"):
        records = read_table(path)
        rows[sample_name(path)].update(contigs=len(records), length=sum(number(r.get("length")) for r in records))
    for path in pick(paths, r"\.coverage\.txt$"):
        records = read_table(path)
        size = sum(number(r["endpos"]) - number(r["startpos"]) + 1 for r in records)
        rows[sample_name(path)].update(
            mapped_reads=sum(number(r["numreads"]) for r in records),
            mean_depth=sum(number(r["meandepth"]) * (number(r["endpos"]) - number(r["startpos"]) + 1) for r in records) / size if size else 0.0,
            breadth=100 * sum(number(r["covbases"]) for r in records) / size if size else 0.0,
        )
    headers = {
        "contigs": fmt("Contigs", "Number of contigs (SeqKit)", "Oranges"),
        "length": fmt("Length", "Total contig length (SeqKit)", "Blues", " bp"),
        "mapped_reads": fmt("Mapped reads", "Reads mapped to the contigs (BWA-MEM2, samtools coverage)", "Purples"),
        "mean_depth": fmt("Mean depth", "Length-weighted mean read depth over all contigs", "Greens", "X", "{:,.1f}"),
        "breadth": fmt("Covered", "Contig bases covered by at least one read", "RdYlGn", "%", "{:,.1f}", max=100, min=0),
    }
    return [table("microcosm_contig_summary", "Contigs: summary", "One row per assembly: contig statistics (SeqKit) and read "
                  "mapping (BWA-MEM2 + samtools coverage). Full tables: ContigsAnalysis/.", headers, rows, "Assembly")]


def pyrodigal(paths: list[Path], _top: int) -> list:
    rows = {}
    for path in pick(paths, r"\.faa(\.gz)?$"):
        genes = complete = residues = 0
        with open_text(path) as handle:
            for line in handle:
                if line.startswith(">"):
                    genes += 1
                    complete += "partial=00" in line
                else:
                    residues += len(line.strip().rstrip("*"))
        rows[sample_name(path)] = {"genes": genes, "complete": complete, "partial": genes - complete,
                                   "mean_length": residues / genes if genes else 0.0}
    headers = {
        "genes": fmt("Genes", "Predicted protein-coding genes", "Greens"),
        "complete": fmt("Complete", "Genes with both ends inside the contig (partial=00)", "Greens"),
        "partial": fmt("Partial", "Genes running off a contig end", "Oranges"),
        "mean_length": fmt("Mean length", "Mean protein length", "Blues", " aa", "{:,.0f}"),
    }
    return [table("microcosm_pyrodigal", "Contigs: genes (Pyrodigal)", "Pyrodigal gene prediction on the contigs. "
                  "Files: ContigsAnalysis/annotation/pyrodigal/.", headers, rows, "Assembly")]


def genomad(paths: list[Path], _top: int) -> list:
    rows: dict[str, dict] = defaultdict(dict)
    for path in pick(paths, r"_plasmid_summary\.tsv$"):
        records = read_table(path)
        rows[sample_name(path)].update(
            plasmids=len(records), plasmid_bp=sum(number(r.get("length")) for r in records),
            plasmid_conjugation=sum(1 for r in records if (r.get("conjugation_genes") or "NA") not in ("", "NA")),
            plasmid_amr=sum(len([g for g in (r.get("amr_genes") or "").split(";") if g and g != "NA"]) for r in records))
    for path in pick(paths, r"_virus_summary\.tsv$"):
        records = read_table(path)
        rows[sample_name(path)].update(
            viruses=len(records), proviruses=sum(1 for r in records if (r.get("topology") or "").lower() == "provirus"),
            virus_bp=sum(number(r.get("length")) for r in records))
    headers = {
        "plasmids": fmt("Plasmids", "geNomad plasmid contigs", "Reds"),
        "plasmid_bp": fmt("Plasmid bp", "Total length of plasmid contigs", "Reds", " bp"),
        "plasmid_conjugation": fmt("Conjugative", "Plasmid contigs with conjugation genes", "Reds"),
        "plasmid_amr": fmt("Plasmid AMR genes", "AMR genes on plasmid contigs (amr_genes)", "Reds"),
        "viruses": fmt("Viruses", "geNomad virus sequences (including proviruses)", "Blues"),
        "proviruses": fmt("Proviruses", "Virus regions inside host contigs", "Blues"),
        "virus_bp": fmt("Virus bp", "Total length of virus sequences", "Blues", " bp"),
    }
    return [table("microcosm_genomad", "Contigs: plasmids and viruses (geNomad)", "geNomad end-to-end on the contigs. "
                  "Files: ContigsAnalysis/genomad/.", headers, rows, "Assembly")]


def tiara_contigs(paths: list[Path], _top: int) -> list:
    lengths = seqkit_lengths(paths)
    counts, sizes = {}, {}
    for path in pick(paths, r"\.tiara\.txt(\.gz)?$"):
        name = sample_name(path)
        count, size = defaultdict(float), defaultdict(float)
        for row in read_table(path):
            label = (row.get("class_fst_stage") or "unknown").strip().lower()
            label = label if label in TIARA_CLASSES else "unknown"
            count[label] += 1
            size[label] += lengths.get(name, {}).get((row.get("sequence_id") or "").split()[0] if row.get("sequence_id") else "", 0.0)
        counts[name] = {label: count[label] for label in TIARA_CLASSES}
        sizes[name] = {label: size[label] for label in TIARA_CLASSES}
    data, labels = [counts], ["Contigs"]
    if any(sum(row.values()) for row in sizes.values()):
        data, labels = [counts, sizes], ["Contigs", "Length (bp)"]
    return [bargraph("microcosm_contig_tiara", "Contigs: domains (Tiara)", "Tiara classification of the contigs (only contigs "
                     "of at least --tiara_min_length are classified). Files: Taxonomy/Tiara/.", TIARA_CLASSES, data, labels)]


def reference_label(row: dict) -> str:
    """'NC_016845.1 Klebsiella pneumoniae subsp. pneumoniae HS11286 chromosome, complete genome' -> species/strain text."""
    words = (row.get("Ref_name") or "").split()
    label = " ".join(words[1:]).split(",")[0].strip()
    label = re.sub(r"\s+(chromosome|plasmid|complete genome|genome assembly).*$", "", label)
    return label or SUFFIX.sub("", Path(row.get("Ref_file") or "unknown").name)


def skani(paths: list[Path], top: int) -> list:
    lengths = seqkit_lengths(paths)
    sizes, counts = {}, {}
    for path in pick(paths, r"\.skani\.tsv$"):
        name = sample_name(path)
        best: dict[str, tuple] = {}
        for row in read_table(path):
            contig = (row.get("Query_name") or "").split()
            if not contig:
                continue
            score = (number(row.get("ANI")), number(row.get("Align_fraction_query")))
            if contig[0] not in best or score > best[contig[0]][0]:
                best[contig[0]] = (score, reference_label(row))
        size, count = defaultdict(float), defaultdict(float)
        for contig, (_score, label) in best.items():
            size[label] += lengths.get(name, {}).get(contig, 0.0)
            count[label] += 1
        if lengths.get(name):
            size[SKANI_NONE] = sum(length for contig, length in lengths[name].items() if contig not in best)
            count[SKANI_NONE] = sum(1 for contig in lengths[name] if contig not in best)
        sizes[name], counts[name] = dict(size), dict(count)
    has_bp = any(value for row in sizes.values() for label, value in row.items() if label != SKANI_NONE)
    categories, ranked = top_categories(sizes if has_bp else counts, top, fixed=(SKANI_NONE,))
    keep = categories[:-2]

    def collapse(values):
        return {name: {**{label: row.get(label, 0.0) for label in keep},
                       OTHER: sum(v for label, v in row.items() if label not in keep and label != SKANI_NONE),
                       SKANI_NONE: row.get(SKANI_NONE, 0.0)} for name, row in values.items()}

    data, labels = ([ranked, collapse(counts)], ["Length (bp)", "Contigs"]) if has_bp else ([ranked], ["Contigs"])
    return [bargraph("microcosm_contig_skani", "Contigs: ANI search (skani)", f"skani search of each contig against --skani_db; "
                     f"every contig counts for its best reference (highest ANI, then aligned fraction). Top {top} references; "
                     f"{SKANI_NONE} = contigs without a hit (short contigs are not searched). Files: ContigsAnalysis/skani/.",
                     categories, data, labels)]


# ---------- bins ----------

def tiara_bins(paths: list[Path], _top: int) -> list:
    counts = {}
    for path in pick(paths, r"\.binclassification\.tsv$"):
        count = defaultdict(float)
        for row in read_table(path):
            label = (row.get("classification") or "unknown").strip().lower()
            count[label if label in TIARA_CLASSES else "unknown"] += 1
        counts[sample_name(path)] = {label: count[label] for label in TIARA_CLASSES}
    return [bargraph("microcosm_bin_tiara", "Bins: domains (Tiara)", "Domain of each bin from the Tiara classification of its "
                     "contigs (--bin_domain_classification). Only prokaryotic bins go on to bin QC, annotation and classification.",
                     TIARA_CLASSES, [counts], ["Bins"])]


def gunc(paths: list[Path], _top: int) -> list:
    rows = {}
    for path in pick(paths, r"maxCSS_level\.tsv$"):
        for row in read_table(path):
            rows[row.get("genome") or "unknown"] = {
                "n_contigs": number(row.get("n_contigs")), "n_genes_called": number(row.get("n_genes_called")),
                "taxonomic_level": row.get("taxonomic_level") or "",
                "clade_separation_score": number(row.get("clade_separation_score")),
                "contamination_portion": number(row.get("contamination_portion")),
                "reference_representation_score": number(row.get("reference_representation_score")),
                "pass": (row.get("pass.GUNC") or "").strip(),
            }
    headers = {
        "n_contigs": fmt("Contigs", "Contigs in the bin", "Oranges"),
        "n_genes_called": fmt("Genes", "Genes called by GUNC", "Greens"),
        "taxonomic_level": {"title": "Level", "description": "Taxonomic level with the highest clade separation score"},
        "clade_separation_score": fmt("CSS", "Clade separation score (> 0.45 = likely chimeric)", "RdYlGn-rev", fmt_="{:,.2f}", max=1, min=0),
        "contamination_portion": fmt("Contamination", "Share of genes outside the major clade", "RdYlGn-rev", fmt_="{:,.2f}", max=1, min=0),
        "reference_representation_score": fmt("RRS", "Reference representation score (low = poorly represented in the GUNC DB)",
                                              "RdYlGn", fmt_="{:,.2f}", max=1, min=0),
        "pass": {"title": "Pass", "description": "pass.GUNC (False = likely chimeric bin)"},
    }
    return [table("microcosm_gunc", "Bins: chimerism (GUNC)", "GUNC run on each bin (--run_gunc). Files: Binning/QC/GUNC/.",
                  headers, rows, "Bin")]


def assembly_stats(paths: list[Path], _top: int) -> list:
    rows = {}
    for path in pick(paths, r"\.assembly_stats\.csv$"):
        for row in read_table(path, ","):
            rows[row.get("bin_name") or "unknown"] = {
                key: number(row.get(key)) for key in ("assembly_length_bp", "contig_count", "contig_N50_bp", "contig_N90_bp", "GC_perc", "gaps_sum_bp")}
    headers = {
        "assembly_length_bp": fmt("Length", "Total bin length", "Blues", " bp"),
        "contig_count": fmt("Contigs", "Contigs (scaffolds split at N runs)", "Oranges"),
        "contig_N50_bp": fmt("N50", "Contig N50", "Greens", " bp"),
        "contig_N90_bp": fmt("N90", "Contig N90", "Greens", " bp"),
        "GC_perc": fmt("GC", "GC content", "PuBu", "%", "{:,.1f}", max=100, min=0),
        "gaps_sum_bp": fmt("Gap bp", "Total length of N runs", "Reds", " bp"),
    }
    return [table("microcosm_assembly_stats", "Bins: assembly statistics", "assembly_stats.py on each bin. "
                  "Files: Binning/QC/AssemblyStats/.", headers, rows, "Bin")]


# ---------- RGI ----------

def family_plot(anchor: str, name: str, description: str, families: dict, top: int, label: str) -> dict | None:
    categories, data = top_categories(families, top)
    return bargraph(anchor, name, description, categories, [data], [label])


def rgi_bwt(paths: list[Path], top: int) -> list:
    rows: dict[str, dict] = defaultdict(dict)
    families: dict[str, dict] = {}
    for path in pick(paths, r"\.overall_mapping_stats\.txt$"):
        with open_text(path) as handle:
            text = handle.read()
        total = re.search(r"Total reads:\s+(\d+)", text)
        mapped = re.search(r"Mapped reads:\s+(\d+)", text)
        if total and mapped:
            rows[sample_name(path)].update(total_reads=number(total.group(1)), mapped_reads=number(mapped.group(1)),
                                           mapped_pct=100 * number(mapped.group(1)) / number(total.group(1)) if number(total.group(1)) else 0.0)
    for path in pick(paths, r"\.gene_mapping_data\.txt$"):
        records = read_table(path)
        name = sample_name(path)
        family = defaultdict(float)
        for row in records:
            family[row.get("AMR Gene Family") or "unknown"] += number(row.get("All Mapped Reads"))
        classes = {c.strip() for row in records for c in (row.get("Drug Class") or "").split(";") if c.strip()}
        rows[name].update(genes=len(records), arg_reads=sum(number(r.get("All Mapped Reads")) for r in records), drug_classes=len(classes))
        families[name] = dict(family)
    headers = {
        "total_reads": fmt("Reads", "Reads given to RGI bwt", "Blues"),
        "mapped_reads": fmt("Mapped", "Reads mapped to CARD references", "Purples"),
        "mapped_pct": fmt("Mapped %", "Share of reads mapped to CARD references", "Purples", "%", "{:,.2f}", max=100, min=0),
        "genes": fmt("ARGs", "CARD reference genes with mapped reads (gene_mapping_data)", "Reds"),
        "arg_reads": fmt("Reads on ARGs", "Mapped reads summed over the reported genes", "Reds"),
        "drug_classes": fmt("Drug classes", "Distinct drug classes of the reported genes", "Oranges"),
    }
    return [
        table("microcosm_rgi_bwt", "Reads: AMR genes (RGI bwt)", "RGI bwt read mapping to CARD (--enable_rgi_bwt). "
              "Files: Readbased_Analysis/rgi_bwt/.", headers, rows, "Sample"),
        family_plot("microcosm_rgi_bwt_families", "Reads: AMR gene families (RGI bwt)", f"Mapped reads per AMR gene family "
                    f"(top {top}).", families, top, "Mapped reads"),
    ]


def sylph(paths: list[Path], top: int) -> list:
    """sylph profile TSVs (one per sample): works without sylph-tax, which MultiQC's sylphtax module needs."""
    rows, abundance = {}, {}
    for path in pick(paths, r"\.tsv$"):
        records = [r for r in read_table(path) if "Taxonomic_abundance" in r]
        if not records:
            continue
        name = sample_name(path)
        genomes = defaultdict(float)
        for row in records:
            genomes[reference_label({"Ref_name": row.get("Contig_name"), "Ref_file": row.get("Genome_file")})] += number(row["Taxonomic_abundance"])
        abundance[name] = dict(genomes)
        rows[name] = {"genomes": len(records), "mean_ani": sum(number(r.get("Adjusted_ANI")) for r in records) / len(records)}
    headers = {
        "genomes": fmt("Genomes", "Reference genomes detected by sylph profile", "Greens"),
        "mean_ani": fmt("Mean ANI", "Mean adjusted ANI of the detected genomes", "RdYlGn", "%", "{:,.2f}", max=100, min=90),
    }
    return [
        table("microcosm_sylph", "Reads: sylph profile", "sylph profile per sample (--readbased_taxotools sylph). "
              "Files: Readbased_Analysis/sylph/profile/.", headers, rows, "Sample"),
        family_plot("microcosm_sylph_abundance", "Reads: sylph taxonomic abundance", f"Taxonomic abundance (%) per reference "
                    f"genome (top {top}).", abundance, top, "Taxonomic abundance (%)"),
    ]


def rgi_amr(paths: list[Path], top: int) -> list:
    rows, families = {}, {}
    for path in pick(paths, r"\.rgi\.txt$"):
        records = read_table(path)
        cut = defaultdict(int)
        family = defaultdict(float)
        for row in records:
            cut[(row.get("Cut_Off") or "").strip()] += 1
            family[row.get("AMR Gene Family") or "unknown"] += 1
        classes = {c.strip() for row in records for c in (row.get("Drug Class") or "").split(";") if c.strip()}
        name = sample_name(path)
        rows[name] = {"hits": len(records), "perfect": cut["Perfect"], "strict": cut["Strict"], "loose": cut["Loose"], "drug_classes": len(classes)}
        families[name] = dict(family)
    headers = {
        "hits": fmt("Hits", "RGI main hits", "Reds"),
        "perfect": fmt("Perfect", "Perfect hits", "Reds"),
        "strict": fmt("Strict", "Strict hits", "Oranges"),
        "loose": fmt("Loose", "Loose hits (only with --include_loose)", "Greys"),
        "drug_classes": fmt("Drug classes", "Distinct drug classes", "Purples"),
    }
    return [
        table("microcosm_rgi_amr", "Bacterial genomes: AMR genes (RGI)", "RGI main on each bacterial genome. "
              "Files: Bacterial_Analysis/amr/rgi/.", headers, rows, "Genome"),
        family_plot("microcosm_rgi_amr_families", "Bacterial genomes: AMR gene families (RGI)", f"RGI hits per AMR gene family "
                    f"(top {top}).", families, top, "Hits"),
    ]


# Kleborate v3 column names depend on the preset (klebsiella_pneumo_complex__, klebsiella_oxytoca_complex__,
# escherichia__): each report column takes the first non-empty value among the matching Kleborate columns.
KLEBORATE_COLUMNS = {
    "species": (r"^enterobacterales__species__species$",),
    "species_match": (r"^enterobacterales__species__species_match$",),
    "st": (r"_complex__mlst__ST$", r"^escherichia__mlst_achtman__ST$"),
    "virulence_score": (r"__virulence_score__virulence_score$",),
    "resistance_score": (r"__resistance_score__resistance_score$",),
    "resistance_classes": (r"__num_resistance_classes$",),
    "resistance_genes": (r"__num_resistance_genes$",),
    "k_locus": (r"__kaptive__K_locus$",),
    "k_confidence": (r"__kaptive__K_locus_confidence$",),
    "o_type": (r"__kaptive__O_type$", r"^escherichia__ectyper__Serotype$"),
    "clermont": (r"^escherichia__ezclermont__Clermont_type$",),
    "pathotype": (r"^escherichia__pathovar__Pathotype$",),
    "yersiniabactin": (r"__ybst__Yersiniabactin$",),
    "colibactin": (r"__cbst__Colibactin$",),
    "aerobactin": (r"__abst__Aerobactin$",),
    "salmochelin": (r"__smst__Salmochelin$",),
    "rmpadc": (r"__rmst__RmpADC$",),
    "rmpa2": (r"__rmpa2__rmpA2$",),
}
KLEBORATE_NUMERIC = {"virulence_score", "resistance_score", "resistance_classes", "resistance_genes"}


def kleborate_value(row: dict, patterns: tuple) -> str:
    for pattern in patterns:
        for column, value in row.items():
            value = (value or "").strip()
            if re.search(pattern, column or "") and value not in ("", "-", "NA", "Not Tested"):
                return value
    return ""


def kleborate(paths: list[Path], _top: int) -> list:
    rows = {}
    for path in pick(paths, r"\.kleborate\.txt$"):
        for row in read_table(path):
            values = {key: kleborate_value(row, patterns) for key, patterns in KLEBORATE_COLUMNS.items()}
            # Scores are 0 for species whose modules were "Not Tested": keep them only when a typing module ran.
            tested = bool(values["st"] or values["k_locus"] or values["o_type"])
            values = {key: number(v) if key in KLEBORATE_NUMERIC else v for key, v in values.items()
                      if v and (tested or key not in KLEBORATE_NUMERIC)}
            if values:
                rows[(row.get("strain") or "").strip() or sample_name(path)] = values
    text = lambda title, description: {"title": title, "description": description, "scale": False}
    headers = {
        "species": text("Species", "Kleborate species"),
        "species_match": text("Match", "Species match strength"),
        "st": text("ST", "Kleborate MLST (KpSC / KoSC scheme, or Achtman for Escherichia)"),
        "virulence_score": fmt("Virulence", "Virulence score (0-5, KpSC)", "Reds", max=5, min=0),
        "resistance_score": fmt("Resistance", "Resistance score (0-3, KpSC)", "Reds", max=3, min=0),
        "resistance_classes": fmt("AMR classes", "Drug classes with acquired resistance", "Oranges"),
        "resistance_genes": fmt("AMR genes", "Acquired resistance genes", "Oranges"),
        "k_locus": text("K locus", "Capsule (K) locus from Kaptive"),
        "k_confidence": text("K confidence", "Kaptive K locus confidence"),
        "o_type": text("O type", "O antigen type (Kaptive) or serotype (ECTyper for Escherichia)"),
        "clermont": text("Clermont", "EzClermont phylogroup (Escherichia)"),
        "pathotype": text("Pathotype", "Escherichia pathotype"),
        "yersiniabactin": text("Yersiniabactin", "ybt lineage / ICEKp"),
        "colibactin": text("Colibactin", "clb lineage"),
        "aerobactin": text("Aerobactin", "iuc lineage"),
        "salmochelin": text("Salmochelin", "iro lineage"),
        "rmpadc": text("RmpADC", "Hypermucoidy locus rmpADC"),
        "rmpa2": text("rmpA2", "Hypermucoidy gene rmpA2"),
    }
    return [table("microcosm_kleborate", "Bacterial genomes: Kleborate", "Kleborate v3 on Klebsiella, Raoultella and "
                  "Escherichia genomes (preset from the lineage). Files: Bacterial_Analysis/klebsiella/kleborate/.",
                  headers, rows, "Genome")]


def salmonella(paths: list[Path], _top: int) -> list:
    """SISTR, SeqSero2, routing, GenoTyphi and Paratype tables, recognised by their headers."""
    rows: dict[str, dict] = defaultdict(dict)
    clean = lambda value: "" if (value or "").strip() in ("", "-", "NA", "missing") else value.strip()
    for path in pick(paths, r"\.(tab|tsv)$"):
        for row in read_table(path):
            if "serovar" in row and "cgmlst_ST" in row:  # SISTR
                antigens = ":".join(clean(row.get(k)) or "-" for k in ("o_antigen", "h1", "h2"))
                rows[clean(row.get("genome")) or sample_name(path)].update(
                    sistr_serovar=clean(row["serovar"]), serogroup=clean(row.get("serogroup")), antigens=antigens,
                    cgmlst_st=clean(row["cgmlst_ST"]), sistr_qc=clean(row.get("qc_status")))
            elif "Predicted serotype" in row:  # SeqSero2
                rows[clean(row.get("Sample name")) or sample_name(path)].update(
                    seqsero2_serotype=clean(row["Predicted serotype"]), seqsero2_profile=clean(row.get("Predicted antigenic profile")))
            elif "routing_reason" in row:  # microcosm Salmonella routing
                route = "GenoTyphi" if row.get("route_genotyphi") == "true" else "Paratype" if row.get("route_paratype") == "true" else "none"
                rows[clean(row.get("genome_id"))].update(route=route, conflict="yes" if row.get("tool_conflict") == "true" else "no")
            elif "final genotype" in row:  # GenoTyphi (Mykrobe)
                rows[clean(row.get("genome")).removesuffix(".mykrobe") or sample_name(path)].update(
                    genotype=clean(row["final genotype"]), genotype_confidence=clean(row.get("confidence")))
            elif "Genotype" in row and "Primary_clade" in row:  # Paratype
                rows[clean(row.get("Strain")) or sample_name(path)].update(
                    genotype=clean(row["Genotype"]), genotype_confidence=f"support {clean(row.get('Support'))}" if clean(row.get("Support")) else "")
    rows = {genome: {k: v for k, v in row.items() if v} for genome, row in rows.items() if genome}
    text = lambda title, description: {"title": title, "description": description, "scale": False}
    headers = {
        "sistr_serovar": text("SISTR serovar", "SISTR serovar (routes the follow-up tools)"),
        "serogroup": text("Serogroup", "SISTR serogroup"),
        "antigens": text("O:H1:H2", "SISTR antigenic formula"),
        "cgmlst_st": text("cgMLST ST", "SISTR cgMLST330 sequence type"),
        "sistr_qc": text("SISTR QC", "SISTR QC status"),
        "seqsero2_serotype": text("SeqSero2", "SeqSero2 predicted serotype"),
        "seqsero2_profile": text("SeqSero2 profile", "SeqSero2 antigenic profile"),
        "conflict": text("Conflict", "SISTR and SeqSero2 disagree"),
        "route": text("Follow-up", "Typhi: Mykrobe + GenoTyphi; Paratyphi A: Paratype"),
        "genotype": text("Genotype", "GenoTyphi (Typhi) or Paratype (Paratyphi A) genotype"),
        "genotype_confidence": text("Genotype support", "GenoTyphi confidence or Paratype support"),
    }
    return [table("microcosm_salmonella", "Bacterial genomes: Salmonella typing", "SISTR and SeqSero2 serovars, the "
                  "follow-up chosen from the SISTR serovar and the GenoTyphi / Paratype genotype. "
                  "Files: Bacterial_Analysis/salmonella/.", headers, rows, "Genome")]


TOOLS = {f.__name__: f for f in (contig_summary, pyrodigal, genomad, tiara_contigs, skani, tiara_bins, gunc, assembly_stats, rgi_bwt, rgi_amr, sylph, kleborate, salmonella)}


def expand(paths: list[Path]) -> list[Path]:
    files = []
    for path in paths:
        # os.walk, not rglob: Nextflow stages folders as symlinks, which rglob does not enter
        files.extend(sorted(Path(root) / name for root, _dirs, names in os.walk(path, followlinks=True) for name in names)
                     if path.is_dir() else [path])
    return files


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--tool", required=True, choices=sorted(TOOLS))
    parser.add_argument("inputs", type=Path, nargs="*")
    parser.add_argument("--top", type=int, default=10, help="categories kept in bar graphs (rest = Other)")
    parser.add_argument("--outdir", type=Path, default=Path("."))
    parser.add_argument("--version", action="version", version=VERSION)
    args = parser.parse_args(argv)

    args.outdir.mkdir(parents=True, exist_ok=True)
    written = []
    for content in TOOLS[args.tool](expand(args.inputs), args.top):
        if content:
            name = f"{content['id'].removeprefix('microcosm_')}_mqc.json"
            (args.outdir / name).write_text(json.dumps(content, indent=1) + "\n")
            written.append(name)
    print(f"{args.tool}: {written or 'no rows'}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
