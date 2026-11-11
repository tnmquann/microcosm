#!/usr/bin/env python3
"""Write a compact per-genome manifest for general bacterial typing."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


VERSION = "1.1.0"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--genome-id", required=True)
    parser.add_argument("--source-route", required=True)
    parser.add_argument("--virulencefinder-tsv", required=True)
    parser.add_argument("--plasmidfinder-tsv", required=True)
    parser.add_argument("--mlst-tsv", required=True)
    parser.add_argument("--blastn-output", required=True)
    parser.add_argument("--blastn-filtered-output", required=True)
    parser.add_argument("--blast-reference-type", required=True)
    parser.add_argument("--blast-reference-label", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--version", action="version", version=VERSION)
    return parser.parse_args()


def nonempty_rows(path: Path) -> list[str]:
    return [
        line.strip()
        for line in path.read_text(encoding="utf-8", errors="replace").splitlines()
        if line.strip()
    ]


def tabular_state(path: Path) -> str:
    rows = nonempty_rows(path)
    return "hits_reported" if len(rows) > 1 else "no_hits"


def mlst_state(path: Path) -> str:
    rows = nonempty_rows(path)
    if not rows:
        return "no_typing_match"
    fields = rows[0].split("\t")
    if len(fields) >= 3 and fields[1] not in {"", "-"} and fields[2] not in {"", "-"}:
        return "typing_result_reported"
    return "no_typing_match"


def blast_state(path: Path) -> str:
    return "matches_reported" if len(nonempty_rows(path)) > 1 else "no_hits"


def is_skipped_output(path: Path) -> bool:
    return ".skipped." in path.name


def build_rows(args: argparse.Namespace) -> list[dict[str, str]]:
    virulence = Path(args.virulencefinder_tsv)
    plasmid = Path(args.plasmidfinder_tsv)
    mlst = Path(args.mlst_tsv)
    blast = Path(args.blastn_output)
    blast_filtered = Path(args.blastn_filtered_output)
    common = {
        "sample_id": args.sample_id,
        "genome_id": args.genome_id,
        "source_route": args.source_route,
        "execution_status": "completed",
    }
    return [
        common
        | {
            "tool": "virulencefinder",
            "execution_status": (
                "skipped" if is_skipped_output(virulence) else "completed"
            ),
            "result_state": (
                "not_configured"
                if is_skipped_output(virulence)
                else tabular_state(virulence)
            ),
            "database_contract": "explicit_read_only_database",
            "raw_output": virulence.name,
            "filtered_output": "",
            "interpretation": (
                "VirulenceFinder skipped because no database was supplied."
                if is_skipped_output(virulence)
                else (
                    "Raw VirulenceFinder calls; database taxon coverage is limited and "
                    "no-hit is not evidence that virulence determinants are absent."
                )
            ),
        },
        common
        | {
            "tool": "plasmidfinder",
            "result_state": tabular_state(plasmid),
            "database_contract": "tool_package_bundled_database",
            "raw_output": plasmid.name,
            "filtered_output": "",
            "interpretation": "Raw PlasmidFinder result; no cross-tool conclusion assigned.",
        },
        common
        | {
            "tool": "mlst",
            "result_state": mlst_state(mlst),
            "database_contract": "tool_package_bundled_schemes",
            "raw_output": mlst.name,
            "filtered_output": "",
            "interpretation": "Raw MLST scheme/type result; no cross-tool conclusion assigned.",
        },
        (
            common | {"execution_status": "skipped"}
            if args.blast_reference_type == "not_configured"
            else common
        )
        | {
            "tool": "blastn",
            "result_state": (
                "not_configured"
                if args.blast_reference_type == "not_configured"
                else blast_state(blast_filtered)
            ),
            "database_contract": args.blast_reference_type,
            "raw_output": blast.name,
            "filtered_output": (
                "" if args.blast_reference_type == "not_configured" else blast_filtered.name
            ),
            "interpretation": (
                "BLASTN skipped because no reference FASTA or prebuilt database was supplied."
                if args.blast_reference_type == "not_configured"
                else (
                    f"Sequence similarity to user-supplied reference '{args.blast_reference_label}' "
                    "only; no cross-tool or biological conclusion assigned."
                )
            ),
        },
    ]


def write_rows(path: Path, rows: list[dict[str, str]]) -> None:
    fieldnames = [
        "sample_id",
        "genome_id",
        "source_route",
        "tool",
        "execution_status",
        "result_state",
        "database_contract",
        "raw_output",
        "filtered_output",
        "interpretation",
    ]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    args = parse_args()
    write_rows(Path(args.output), build_rows(args))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
