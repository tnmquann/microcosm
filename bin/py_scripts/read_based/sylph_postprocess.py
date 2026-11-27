#!/usr/bin/env python3
"""Merge sylph-tax profiles (*.sylphmpa) and build microcosm-style phyloseq tables.

Inputs: a manifest TSV (sample_id, sylphmpa path; no header) and the
abundance column to use (--datatype taxonomic_abundance -> relative_abundance,
sequence_abundance -> sequence_abundance, i.e. read counts when sylph ran with
--estimate-read-counts).

Outputs (in --outdir):
  <prefix>_sylphtax_<datatype>.tsv  clade_name x sample matrix, every rank
  phyloseq/ (with --phyloseq)       tax_table.csv, otu_table.csv, sample_data.csv
                                    in the layout of post_processing.py: taxa are
                                    genomes (t__ rows), lineage d__..s__.
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

VERSION = "1.0.0"
DATATYPE_COLUMNS = {
    "taxonomic_abundance": "relative_abundance",
    "sequence_abundance": "sequence_abundance",
}
RANKS = ["domain", "phylum", "class", "order", "family", "genus", "species"]


def read_sylphmpa(path: Path, column: str) -> dict[str, str]:
    values: dict[str, str] = {}
    header: list[str] | None = None
    with path.open(newline="") as handle:
        for row in csv.reader(handle, delimiter="\t"):
            if not row or row[0].startswith("#"):
                continue
            if header is None:
                header = row
                if column not in header:
                    raise SystemExit(f"{path}: no '{column}' column (found {header})")
                continue
            values[row[0]] = row[header.index(column)]
    return values


def number(text: str) -> str:
    value = float(text)
    return str(int(value)) if value.is_integer() else format(value, ".15g")


def run(manifest: Path, datatype: str, outdir: Path, prefix: str, phyloseq: bool) -> int:
    column = DATATYPE_COLUMNS[datatype]
    samples: list[tuple[str, Path]] = []
    with manifest.open(newline="") as handle:
        for row in csv.reader(handle, delimiter="\t"):
            if row:
                samples.append((row[0], Path(row[1])))
    samples.sort()
    profiles = {sample: read_sylphmpa(path, column) for sample, path in samples}
    clades = sorted({clade for values in profiles.values() for clade in values})

    outdir.mkdir(parents=True, exist_ok=True)
    with (outdir / f"{prefix}_sylphtax_{datatype}.tsv").open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["clade_name"] + [sample for sample, _ in samples])
        for clade in clades:
            writer.writerow([clade] + [number(profiles[s].get(clade, "0")) for s, _ in samples])

    if phyloseq:
        genomes = [clade for clade in clades if clade.rsplit("|", 1)[-1].startswith("t__")]
        target = outdir / "phyloseq"
        target.mkdir(exist_ok=True)
        with (target / "tax_table.csv").open("w", newline="") as handle:
            writer = csv.writer(handle, lineterminator="\n")
            writer.writerow(["match_name"] + RANKS)
            for clade in genomes:
                parts = clade.split("|")
                ranks = (parts[:-1] + [""] * len(RANKS))[: len(RANKS)]
                writer.writerow([parts[-1][3:]] + ranks)
        with (target / "otu_table.csv").open("w", newline="") as handle:
            writer = csv.writer(handle, lineterminator="\n")
            writer.writerow(["match_name"] + [sample for sample, _ in samples])
            for clade in genomes:
                writer.writerow([clade.rsplit("|", 1)[-1][3:]] + [number(profiles[s].get(clade, "0")) for s, _ in samples])
        with (target / "sample_data.csv").open("w", newline="") as handle:
            writer = csv.writer(handle, lineterminator="\n")
            writer.writerow(["sample_name", "source_file"])
            for sample, path in samples:
                writer.writerow([sample, path.name])
    return len(clades)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--datatype", choices=sorted(DATATYPE_COLUMNS), default="sequence_abundance")
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--prefix", required=True)
    parser.add_argument("--phyloseq", action="store_true")
    parser.add_argument("--version", action="version", version=VERSION)
    args = parser.parse_args(argv)
    count = run(args.manifest, args.datatype, args.outdir, args.prefix, args.phyloseq)
    print(f"{count} clades written to {args.outdir}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
