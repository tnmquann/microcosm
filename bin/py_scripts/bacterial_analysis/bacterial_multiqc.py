#!/usr/bin/env python3
"""MultiQC table: one row per bacterial genome of the downstream layer.

Columns: sample, GTDB species, MLST (scheme:ST), Salmonella serovar (SISTR),
and the number of distinct AMR genes per tool from the hAMRonization summary.
Only --routing is required; every other table is optional.
"""

from __future__ import annotations

import argparse
import csv
import sys
from collections import defaultdict
from pathlib import Path

VERSION = "1.0.0"


def read_tsv(path: Path | None) -> list[dict]:
    if not path or not path.exists() or path.stat().st_size == 0:
        return []
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def build_rows(routing: list[dict], mlst: list[dict], salmonella: list[dict], amr: list[dict]) -> tuple[list[str], dict[str, dict]]:
    tools = sorted({row["analysis_software_name"] for row in amr if row.get("analysis_software_name")})
    genes: dict[str, dict[str, set]] = defaultdict(lambda: defaultdict(set))
    for row in amr:
        genes[row.get("input_file_name", "")][row.get("analysis_software_name", "")].add(row.get("gene_symbol", ""))
    st = {row["genome_id"]: f"{row.get('scheme', '-')}:{row.get('st', '-')}" for row in mlst if row.get("scheme", "-") != "-"}
    serovar = {row["genome_id"]: row.get("predicted_serovar_or_lineage", "") for row in salmonella if row.get("tool") == "sistr"}
    columns = ["sample_id", "species", "mlst", "serovar"] + [f"amr_{tool}" for tool in tools]
    rows = {}
    for row in routing:
        genome = row["genome_id"]
        values = {
            "sample_id": row.get("sample_id", ""),
            "species": (row.get("species") or "").removeprefix("s__"),
            "mlst": st.get(genome, ""),
            "serovar": serovar.get(genome, ""),
        }
        for tool in tools:
            values[f"amr_{tool}"] = len(genes[genome][tool] - {""})
        rows[genome] = values
    return columns, rows


def write_table(path: Path, columns: list[str], rows: dict[str, dict], tools_label: str) -> None:
    header = [
        "# id: 'microcosm_bacterial_overview'",
        "# section_name: 'Bacterial genomes'",
        f"# description: 'Genomes routed to the bacterial downstream analysis: GTDB species, MLST, SISTR serovar and distinct AMR genes per tool ({tools_label}). Full tables: Bacterial_Analysis/*/summary.'",
        "# plot_type: 'table'",
        "# pconfig:",
        "#     id: 'microcosm_bacterial_overview_table'",
        "#     namespace: 'Bacterial genomes'",
        "#     col1_header: 'Genome'",
    ]
    with path.open("w", newline="") as handle:
        handle.write("\n".join(header) + "\n")
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["Genome"] + columns)
        for genome in sorted(rows):
            writer.writerow([genome] + [rows[genome][column] for column in columns])


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--routing", type=Path, required=True, help="routing/bacterial_routing_manifest.tsv")
    parser.add_argument("--mlst", type=Path, help="general_typing summary <prefix>_mlst.tsv")
    parser.add_argument("--salmonella", type=Path, help="salmonella summary <prefix>_salmonella_typing.tsv")
    parser.add_argument("--amr", type=Path, help="hamronization_summary/amr_summary.tsv")
    parser.add_argument("--output", type=Path, default=Path("bacterial_overview_mqc.tsv"))
    parser.add_argument("--version", action="version", version=VERSION)
    args = parser.parse_args(argv)
    routing = read_tsv(args.routing)
    columns, rows = build_rows(routing, read_tsv(args.mlst), read_tsv(args.salmonella), read_tsv(args.amr))
    tools = [column.removeprefix("amr_") for column in columns if column.startswith("amr_")]
    write_table(args.output, columns, rows, ", ".join(tools) or "none run")
    print(f"{len(rows)} genomes -> {args.output}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
