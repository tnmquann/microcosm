#!/usr/bin/env python3
"""Normalize SeqSero2/SISTR Salmonella evidence and route follow-up tools.

Routing follows playground/modular/4_sal.sh: only the SISTR ``serovar`` column
decides follow-up, by exact value: ``Typhi`` -> Mykrobe + GenoTyphi,
``Paratyphi A`` -> Paratype, anything else -> none.  SeqSero2 is evidence only.
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path
from typing import Any


VERSION = "2.0.0"
EVIDENCE_FIELDS = (
    "sample_id",
    "genome_id",
    "tool",
    "predicted_serovar_or_lineage",
    "normalized_prediction",
    "confidence_or_status",
    "raw_result_path",
)
ROUTING_FIELDS = (
    "routing_key",
    "sample_id",
    "genome_id",
    "sistr_serovar",
    "seqsero2_prediction",
    "route_genotyphi",
    "route_paratype",
    "tool_conflict",
    "routing_reason",
)

# Exact normalizations only.  In particular, no substring matching is used.
SEROVAR_NORMALIZATION = {
    "typhi": "Typhi",
    "salmonella typhi": "Typhi",
    "salmonella enterica serovar typhi": "Typhi",
    "salmonella enterica subsp. enterica serovar typhi": "Typhi",
    "paratyphi a": "Paratyphi A",
    "salmonella paratyphi a": "Paratyphi A",
    "salmonella enterica serovar paratyphi a": "Paratyphi A",
    "salmonella enterica subsp. enterica serovar paratyphi a": "Paratyphi A",
    "typhimurium": "Typhimurium",
    "salmonella typhimurium": "Typhimurium",
    "salmonella enterica serovar typhimurium": "Typhimurium",
}
SEQSERO2_SEROVAR_FIELDS = (
    "Predicted serotype",
    "Predicted serovar",
    "predicted_serotype",
    "predicted_serovar",
)
SEQSERO2_STATUS_FIELDS = ("Note", "note", "Status", "status")
SISTR_SEROVAR_FIELDS = ("serovar", "Serovar", "predicted_serovar", "Predicted serovar")
SISTR_STATUS_FIELDS = ("serogroup", "Serogroup", "cgmlst_subspecies", "cgmlst_ST")

def parse_args() -> argparse.Namespace:
    """Parse the fixed, staged-input command-line contract."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--genome-id", required=True)
    parser.add_argument("--seqsero2-tsv", type=Path, help="omit when SeqSero2 was skipped")
    parser.add_argument("--sistr-tsv", required=True, type=Path)
    parser.add_argument("--evidence-output", required=True, type=Path)
    parser.add_argument("--routing-output", required=True, type=Path)
    parser.add_argument("--version", action="version", version=VERSION)
    return parser.parse_args()


def canonical_serovar(value: Any) -> str:
    """Return a recognized serovar only when ``value`` is an exact known form."""
    if not isinstance(value, str):
        return ""
    return SEROVAR_NORMALIZATION.get(" ".join(value.split()).casefold(), "")


def first_named_value(row: dict[str, str], names: tuple[str, ...]) -> str:
    """Read the first non-empty value from a fixed list of named columns."""
    for name in names:
        value = (row.get(name) or "").strip()
        if value:
            return value
    return ""


def read_single_tsv(path: Path, label: str) -> dict[str, str] | None:
    """Read at most one TSV row without silently collapsing duplicate columns."""
    with path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            return None
        if len(reader.fieldnames) != len(set(reader.fieldnames)):
            raise ValueError(f"{label} result contains duplicate column names: {path}")
        rows = list(reader)
    if not rows:
        return None
    if len(rows) != 1:
        raise ValueError(f"{label} result must contain one data row: {path}")
    if None in rows[0]:
        raise ValueError(
            f"{label} result row contains more values than the header defines: {path}"
        )
    return rows[0]


def parse_seqsero2(path: Path | None) -> tuple[str, str]:
    """Read SeqSero2's named prediction and status columns from one TSV row."""
    if path is None:
        return "", "not_run"
    row = read_single_tsv(path, "SeqSero2")
    if row is None:
        return "", "empty_result"
    prediction = first_named_value(row, SEQSERO2_SEROVAR_FIELDS)
    if not prediction:
        expected = ", ".join(SEQSERO2_SEROVAR_FIELDS)
        raise ValueError(
            f"SeqSero2 result lacks a supported predicted-serovar column ({expected}): {path}"
        )
    return prediction, first_named_value(row, SEQSERO2_STATUS_FIELDS) or "reported"


def parse_sistr(path: Path) -> tuple[str, str]:
    """Read SISTR's named serovar field; it is the only routing input."""
    row = read_single_tsv(path, "SISTR")
    if row is None:
        return "", "empty_result"
    prediction = first_named_value(row, SISTR_SEROVAR_FIELDS)
    status = first_named_value(row, SISTR_STATUS_FIELDS)
    return prediction, status or ("reported" if prediction else "prediction_not_reported")


def write_tsv(path: Path, fields: tuple[str, ...], rows: list[dict[str, str]]) -> None:
    """Write a deterministic TSV report."""
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def build_reports(
    sample_id: str,
    genome_id: str,
    seqsero2_tsv: Path | None,
    sistr_tsv: Path,
) -> tuple[list[dict[str, str]], dict[str, str]]:
    """Build normalized evidence and the SISTR-only route decision."""
    seqsero2_prediction, seqsero2_status = parse_seqsero2(seqsero2_tsv)
    sistr_prediction, sistr_status = parse_sistr(sistr_tsv)
    seqsero2_normalized = canonical_serovar(seqsero2_prediction)
    sistr_normalized = canonical_serovar(sistr_prediction)
    conflict = bool(
        seqsero2_normalized and sistr_normalized and seqsero2_normalized != sistr_normalized
    )
    # Exact match on the raw SISTR value, as `case "$serovar"` in 4_sal.sh.
    route_genotyphi = sistr_prediction == "Typhi"
    route_paratype = sistr_prediction == "Paratyphi A"

    evidence = [
        {
            "sample_id": sample_id,
            "genome_id": genome_id,
            "tool": tool,
            "predicted_serovar_or_lineage": prediction,
            "normalized_prediction": normalized,
            "confidence_or_status": status,
            "raw_result_path": path.name if path else "",
        }
        for tool, prediction, normalized, status, path in (
            ("seqsero2", seqsero2_prediction, seqsero2_normalized, seqsero2_status, seqsero2_tsv),
            ("sistr", sistr_prediction, sistr_normalized, sistr_status, sistr_tsv),
        )
    ]
    if route_genotyphi:
        reasons = ["sistr_serovar_Typhi"]
    elif route_paratype:
        reasons = ["sistr_serovar_Paratyphi_A"]
    elif sistr_prediction:
        reasons = ["sistr_serovar_has_no_follow_up"]
    else:
        reasons = ["sistr_serovar_missing"]
    if conflict:
        reasons.append("seqsero2_sistr_conflict_preserved")
    routing = {
        "routing_key": f"{sample_id}::{genome_id}",
        "sample_id": sample_id,
        "genome_id": genome_id,
        "sistr_serovar": sistr_prediction,
        "seqsero2_prediction": seqsero2_normalized or seqsero2_prediction,
        "route_genotyphi": str(route_genotyphi).lower(),
        "route_paratype": str(route_paratype).lower(),
        "tool_conflict": str(conflict).lower(),
        "routing_reason": ";".join(reasons),
    }
    return evidence, routing


def main() -> int:
    """Run the parser and write both reports."""
    args = parse_args()
    try:
        evidence, routing = build_reports(
            args.sample_id,
            args.genome_id,
            args.seqsero2_tsv,
            args.sistr_tsv,
        )
        write_tsv(args.evidence_output, EVIDENCE_FIELDS, evidence)
        write_tsv(args.routing_output, ROUTING_FIELDS, [routing])
    except (OSError, csv.Error, ValueError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
