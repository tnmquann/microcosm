#!/usr/bin/env python3
"""Write source-preserving general-typing evidence and MLST summaries."""

from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path
from typing import Iterable


VERSION = "1.0.0"
EVIDENCE_FIELDS = (
    "sample_id",
    "genome_id",
    "source_tool",
    "observation_type",
    "normalized_identifier",
    "raw_identifier",
    "source_output",
)
MLST_FIELDS = (
    "sample_id",
    "genome_id",
    "scheme",
    "sequence_type",
    "allele_profile",
    "clonal_complex",
    "clonal_complex_lookup_status",
    "source_output",
)


def normalize_identifier(value: str) -> str:
    """Normalize presentation only, without collapsing allele/version syntax."""
    return re.sub(r"\s+", "_", value.strip().casefold())


def tabular_identifiers(path: Path, candidate_fields: tuple[str, ...]) -> list[str]:
    """Read identifiers from the first supported non-empty TSV field."""
    if not path.is_file() or ".skipped." in path.name:
        return []
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames:
            return []
        field = next((name for name in candidate_fields if name in reader.fieldnames), None)
        if field is None:
            return []
        return [
            row[field].strip()
            for row in reader
            if row.get(field) and row[field].strip() not in {"-", "NA", "Unknown"}
        ]


def build_evidence_rows(
    sample_id: str,
    genome_id: str,
    virulencefinder: Path,
    plasmidfinder: Path,
    blastn_filtered: Path,
) -> list[dict[str, str]]:
    """Build a deterministic long table without cross-tool consensus calls."""
    sources = (
        (
            "virulencefinder",
            "gene",
            virulencefinder,
            ("Virulence factor", "Gene", "gene", "Gene ID"),
        ),
        (
            "plasmidfinder",
            "replicon",
            plasmidfinder,
            ("Plasmid", "plasmid", "Database"),
        ),
        ("blastn", "reference_hit", blastn_filtered, ("sseqid",)),
    )
    rows: list[dict[str, str]] = []
    for tool, observation_type, path, fields in sources:
        for raw_identifier in tabular_identifiers(path, fields):
            rows.append(
                {
                    "sample_id": sample_id,
                    "genome_id": genome_id,
                    "source_tool": tool,
                    "observation_type": observation_type,
                    "normalized_identifier": normalize_identifier(raw_identifier),
                    "raw_identifier": raw_identifier,
                    "source_output": path.name,
                }
            )
    unique = {tuple(row[field] for field in EVIDENCE_FIELDS): row for row in rows}
    return [unique[key] for key in sorted(unique)]


def value_or_unknown(value: str) -> str:
    value = value.strip()
    return "Unknown" if value in {"", "-", "NA", "N/A"} else value


def build_mlst_rows(sample_id: str, genome_id: str, mlst: Path) -> list[dict[str, str]]:
    """Preserve MLST scheme, ST, and alleles; do not infer a clonal complex."""
    parsed: list[dict[str, str]] = []
    if mlst.is_file():
        with mlst.open(encoding="utf-8", newline="") as handle:
            for row_number, row in enumerate(csv.reader(handle, delimiter="\t"), start=1):
                if not row:
                    continue
                if len(row) < 3:
                    raise ValueError(
                        f"MLST row {row_number} has {len(row)} columns; at least 3 are required"
                    )
                parsed.append(
                    {
                        "sample_id": sample_id,
                        "genome_id": genome_id,
                        "scheme": value_or_unknown(row[1]),
                        "sequence_type": value_or_unknown(row[2]),
                        "allele_profile": (
                            ";".join(value_or_unknown(value) for value in row[3:])
                            if row[3:]
                            else "Unknown"
                        ),
                        "clonal_complex": "NA",
                        "clonal_complex_lookup_status": "not_configured",
                        "source_output": mlst.name,
                    }
                )
    if parsed:
        return sorted(
            parsed,
            key=lambda row: (
                row["scheme"],
                row["sequence_type"],
                row["allele_profile"],
            ),
        )
    return [
        {
            "sample_id": sample_id,
            "genome_id": genome_id,
            "scheme": "Unknown",
            "sequence_type": "Unknown",
            "allele_profile": "Unknown",
            "clonal_complex": "NA",
            "clonal_complex_lookup_status": "not_configured",
            "source_output": mlst.name,
        }
    ]


def write_tsv(path: Path, fields: Iterable[str], rows: Iterable[dict[str, str]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(fields), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--genome-id", required=True)
    parser.add_argument("--virulencefinder", required=True, type=Path)
    parser.add_argument("--plasmidfinder", required=True, type=Path)
    parser.add_argument("--mlst", required=True, type=Path)
    parser.add_argument("--blastn-filtered", required=True, type=Path)
    parser.add_argument("--output", required=True)
    parser.add_argument("--version", action="version", version=VERSION)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    write_tsv(
        Path(f"{args.output}.general_typing_evidence.tsv"),
        EVIDENCE_FIELDS,
        build_evidence_rows(
            args.sample_id,
            args.genome_id,
            args.virulencefinder,
            args.plasmidfinder,
            args.blastn_filtered,
        ),
    )
    write_tsv(
        Path(f"{args.output}.mlst_summary.tsv"),
        MLST_FIELDS,
        build_mlst_rows(args.sample_id, args.genome_id, args.mlst),
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
