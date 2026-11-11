#!/usr/bin/env python3
"""Normalize explicit Salmonella typing evidence for safe follow-up routing.

The script deliberately reads only documented, named result fields.  It does
not scan arbitrary text for the word ``Typhi``: this prevents a Typhimurium
call from activating the Typhi-only GenoTyphi workflow.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from pathlib import Path
from typing import Any


VERSION = "1.0.0"
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
    "route_genotyphi",
    "route_paratype",
    "tool_conflict",
    "seqsero2_prediction",
    "mykrobe_prediction",
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
MYKROBE_TYPHI_SPECIES_KEYS = {
    "typhi",
    "salmonella typhi",
    "salmonella_typhi",
    "salmonella enterica serovar typhi",
    "salmonella_enterica_serovar_typhi",
}


def parse_args() -> argparse.Namespace:
    """Parse the fixed, staged-input command-line contract."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--genome-id", required=True)
    parser.add_argument("--seqsero2-tsv", required=True, type=Path)
    parser.add_argument("--mykrobe-json", required=True, type=Path)
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


def parse_seqsero2(path: Path) -> tuple[str, str]:
    """Read SeqSero2's named prediction and status columns from one TSV row."""
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
    """Preserve SISTR's named serovar field without using it for follow-up gates."""
    row = read_single_tsv(path, "SISTR")
    if row is None:
        return "", "empty_result"
    prediction = first_named_value(row, SISTR_SEROVAR_FIELDS)
    status = first_named_value(row, SISTR_STATUS_FIELDS)
    return prediction, status or ("reported" if prediction else "prediction_not_reported")


def sample_payload(payload: Any) -> dict[str, Any]:
    """Return Mykrobe's result object without guessing from arbitrary nesting."""
    if not isinstance(payload, dict):
        raise ValueError("Mykrobe JSON root must be an object")
    if "species" in payload or "phylogenetics" in payload:
        return payload
    object_values = [value for value in payload.values() if isinstance(value, dict)]
    if len(payload) == 1 and len(object_values) == 1:
        return object_values[0]
    raise ValueError(
        "Mykrobe JSON must contain a result object with 'species' or 'phylogenetics'"
    )


def numeric_coverage(value: Any) -> float | None:
    """Coerce a documented coverage value while rejecting malformed values."""
    try:
        numeric = float(value)
    except (TypeError, ValueError):
        return None
    return numeric if math.isfinite(numeric) else None


def unique_json_object(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    """Build a JSON object while rejecting duplicate keys at every level."""
    result: dict[str, Any] = {}
    for key, value in pairs:
        if key in result:
            raise ValueError(f"Mykrobe JSON contains duplicate key {key!r}")
        result[key] = value
    return result


def parse_mykrobe(path: Path) -> tuple[str, str]:
    """Extract Typhi evidence from the structured Mykrobe species object.

    Mykrobe is run with its Typhi panel, so a panel name alone is not enough
    evidence.  An exact Typhi species key must additionally report positive
    ``percent_coverage`` (or ``median_depth`` in older output variants).
    """
    with path.open(encoding="utf-8") as handle:
        result = sample_payload(json.load(handle, object_pairs_hook=unique_json_object))
    species = result.get("species")
    if not isinstance(species, dict):
        return "", "species_not_reported"

    matching_entries = []
    for key, metrics in species.items():
        canonical_key = " ".join(str(key).replace("_", " ").split()).casefold()
        if canonical_key in MYKROBE_TYPHI_SPECIES_KEYS:
            matching_entries.append((key, metrics))
    if not matching_entries:
        return "", "typhi_species_not_reported"
    if len(matching_entries) > 1:
        keys = ", ".join(repr(key) for key, _metrics in matching_entries)
        raise ValueError(f"Mykrobe JSON contains multiple matching species keys: {keys}")

    _key, metrics = matching_entries[0]
    if not isinstance(metrics, dict):
        return "", "typhi_species_entry_malformed"
    coverage = numeric_coverage(metrics.get("percent_coverage"))
    metric_name = "percent_coverage"
    if coverage is None:
        coverage = numeric_coverage(metrics.get("median_depth"))
        metric_name = "median_depth"
    if coverage is not None and coverage > 0:
        return "Typhi", f"{metric_name}={coverage:g}"
    return "", f"{metric_name}_not_positive"


def write_tsv(path: Path, fields: tuple[str, ...], rows: list[dict[str, str]]) -> None:
    """Write a deterministic TSV report."""
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def build_reports(
    sample_id: str,
    genome_id: str,
    seqsero2_tsv: Path,
    mykrobe_json: Path,
    sistr_tsv: Path,
) -> tuple[list[dict[str, str]], dict[str, str]]:
    """Build normalized evidence and route decisions from exact parsed fields."""
    seqsero2_prediction, seqsero2_status = parse_seqsero2(seqsero2_tsv)
    mykrobe_prediction, mykrobe_status = parse_mykrobe(mykrobe_json)
    sistr_prediction, sistr_status = parse_sistr(sistr_tsv)
    seqsero2_normalized = canonical_serovar(seqsero2_prediction)
    mykrobe_normalized = canonical_serovar(mykrobe_prediction)
    sistr_normalized = canonical_serovar(sistr_prediction)
    normalized_predictions = [
        value
        for value in (seqsero2_normalized, sistr_normalized, mykrobe_normalized)
        if value
    ]
    conflict = len(set(normalized_predictions)) > 1
    typhi_evidence = seqsero2_normalized == "Typhi" or mykrobe_normalized == "Typhi"
    paratyphi_a_evidence = seqsero2_normalized == "Paratyphi A"
    routing_key = f"{sample_id}::{genome_id}"

    evidence = [
        {
            "sample_id": sample_id,
            "genome_id": genome_id,
            "tool": "seqsero2",
            "predicted_serovar_or_lineage": seqsero2_prediction,
            "normalized_prediction": seqsero2_normalized,
            "confidence_or_status": seqsero2_status,
            "raw_result_path": seqsero2_tsv.name,
        },
        {
            "sample_id": sample_id,
            "genome_id": genome_id,
            "tool": "sistr",
            "predicted_serovar_or_lineage": sistr_prediction,
            "normalized_prediction": sistr_normalized,
            "confidence_or_status": sistr_status,
            "raw_result_path": sistr_tsv.name,
        },
        {
            "sample_id": sample_id,
            "genome_id": genome_id,
            "tool": "mykrobe",
            "predicted_serovar_or_lineage": mykrobe_prediction,
            "normalized_prediction": mykrobe_normalized,
            "confidence_or_status": mykrobe_status,
            "raw_result_path": mykrobe_json.name,
        },
    ]
    reasons = []
    if typhi_evidence:
        reasons.append("explicit_Typhi_evidence")
    if paratyphi_a_evidence:
        reasons.append("explicit_Paratyphi_A_from_SeqSero2")
    if conflict:
        reasons.append("tool_conflict_preserved")
    if not reasons:
        reasons.append("no_follow_up_evidence")
    routing = {
        "routing_key": routing_key,
        "sample_id": sample_id,
        "genome_id": genome_id,
        "route_genotyphi": str(typhi_evidence).lower(),
        "route_paratype": str(paratyphi_a_evidence).lower(),
        "tool_conflict": str(conflict).lower(),
        "seqsero2_prediction": seqsero2_normalized or seqsero2_prediction,
        "mykrobe_prediction": mykrobe_normalized or mykrobe_prediction,
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
            args.mykrobe_json,
            args.sistr_tsv,
        )
        write_tsv(args.evidence_output, EVIDENCE_FIELDS, evidence)
        write_tsv(args.routing_output, ROUTING_FIELDS, [routing])
    except (OSError, csv.Error, json.JSONDecodeError, ValueError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
