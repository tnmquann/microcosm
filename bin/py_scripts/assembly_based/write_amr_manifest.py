#!/usr/bin/env python3
"""Write one deterministic per-genome AMR execution manifest."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

VERSION = "0.3.0"
HEADERS = (
    "sample_id",
    "genome_id",
    "source_route",
    "tool",
    "request_status",
    "execution_status",
    "result_status",
    "raw_output",
    "harmonization_status",
    "normalization_status",
    "database_path",
    "database_status",
    "database_details",
    "notes",
)


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--genome-id", required=True)
    parser.add_argument("--source-route", required=True)
    parser.add_argument("--amrfinder-report", required=True, type=Path)
    parser.add_argument("--abricate-report", required=True, type=Path)
    parser.add_argument("--rgi-report", required=True, type=Path)
    parser.add_argument("--resfinder-report", required=True, type=Path)
    parser.add_argument("--deeparg-report", required=True, type=Path)
    parser.add_argument("--pointfinder-species", default="")
    parser.add_argument("--abricate-normalizable", action="store_true")
    parser.add_argument("--database-validation", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--version", action="version", version=VERSION)
    return parser.parse_args()


def build_rows(args: argparse.Namespace) -> list[dict[str, str]]:
    """Build explicit execution/harmonization rows without biological inference."""
    common = {
        "sample_id": args.sample_id,
        "genome_id": args.genome_id,
        "source_route": args.source_route,
    }

    database_rows = load_database_validation(
        getattr(args, "database_validation", None)
    )

    def row(
        tool: str,
        output: str,
        harmonization: str,
        normalization: str,
        notes: str,
        execution: str = "completed",
        result: str = "completed",
    ) -> dict[str, str]:
        database = database_rows.get(tool, {})
        return {
            **common,
            "tool": tool,
            "request_status": "requested",
            "execution_status": execution,
            "result_status": result,
            "raw_output": output,
            "harmonization_status": harmonization,
            "normalization_status": normalization,
            "database_path": database.get("path", "unknown"),
            "database_status": database.get("status", "unknown"),
            "database_details": database.get("details", "unknown"),
            "notes": notes,
        }

    def skipped(path: Path) -> bool:
        return ".skipped." in path.name

    def skipped_notes(path: Path, fallback: str) -> str:
        if not skipped(path):
            return fallback
        text = path.read_text(encoding="utf-8", errors="replace").strip()
        return text or fallback

    rows = [
        row(
            "AMRFinderPlus",
            args.amrfinder_report.name,
            "not_run" if skipped(args.amrfinder_report) else "supported",
            "not_run" if skipped(args.amrfinder_report) else "supported",
            skipped_notes(
                args.amrfinder_report,
                "Raw, hAMRonization, and ARGnorm outputs are retained.",
            ),
            execution="skipped" if skipped(args.amrfinder_report) else "completed",
            result=result_status(args.amrfinder_report, "tabular"),
        ),
        row(
            "ABRicate",
            args.abricate_report.name,
            "supported",
            "supported" if args.abricate_normalizable else "not_supported",
            (
                "ARGnorm supports this ABRicate database identifier."
                if args.abricate_normalizable
                else "Raw and hAMRonization outputs are retained; ARGnorm does not "
                "support the configured ABRicate database identifier."
            ),
            result=result_status(args.abricate_report, "tabular"),
        ),
        row(
            "RGI",
            args.rgi_report.name,
            "not_run" if skipped(args.rgi_report) else "supported",
            "not_supported",
            skipped_notes(args.rgi_report, "Raw and hAMRonization outputs are retained."),
            execution="skipped" if skipped(args.rgi_report) else "completed",
            result=result_status(args.rgi_report, "tabular"),
        ),
        row(
            "ResFinder",
            args.resfinder_report.name,
            "not_supported",
            "not_supported",
            skipped_notes(
                args.resfinder_report,
                "Current requested integration preserves ResFinder raw outputs.",
            ),
            execution="skipped" if skipped(args.resfinder_report) else "completed",
            result=result_status(args.resfinder_report, "json"),
        ),
        row(
            "DeepARG",
            args.deeparg_report.name,
            "not_run" if skipped(args.deeparg_report) else "supported",
            "not_run" if skipped(args.deeparg_report) else "supported",
            skipped_notes(
                args.deeparg_report,
                "DeepARG runs on predicted nucleotide ORFs with the LS model.",
            ),
            execution="skipped" if skipped(args.deeparg_report) else "completed",
            result=result_status(args.deeparg_report, "tabular"),
        ),
    ]
    if skipped(args.resfinder_report):
        rows.append(
            row(
                "PointFinder",
                "",
                "not_run",
                "not_run",
                "PointFinder skipped because ResFinder/PointFinder databases were not configured.",
                execution="skipped",
                result="skipped",
            )
        )
    elif args.pointfinder_species:
        rows.append(
            row(
                "PointFinder",
                "PointFinder outputs from RESFINDER_AMR",
                "not_supported",
                "not_supported",
                f"Enabled with parsed taxonomy '{args.pointfinder_species}'.",
                result="completed",
            )
        )
    else:
        rows.append(
            row(
                "PointFinder",
                "",
                "not_supported",
                "not_supported",
                "No parsed species/genus matched a configured PointFinder database.",
                execution="skipped",
                result="skipped",
            )
        )
    return rows


def load_database_validation(path: Path | None) -> dict[str, dict[str, str]]:
    """Load validator-produced database provenance keyed by tool name."""
    if path is None or not path.is_file():
        return {}
    with path.open(encoding="utf-8", newline="") as handle:
        return {
            row["tool"]: row
            for row in csv.DictReader(handle, delimiter="\t")
            if row.get("tool")
        }


def result_status(path: Path, output_format: str) -> str:
    """Classify only explicit skipped/empty/non-empty result states."""
    if ".skipped." in path.name:
        return "skipped"
    if not path.is_file() or path.stat().st_size == 0:
        return "no_hits"
    if output_format == "json":
        try:
            value = json.loads(path.read_text(encoding="utf-8"))
        except (OSError, UnicodeDecodeError, json.JSONDecodeError):
            return "completed"
        return "no_hits" if value in ({}, []) else "completed"

    lines = [
        line.strip()
        for line in path.read_text(encoding="utf-8", errors="replace").splitlines()
        if line.strip() and not line.lstrip().startswith("#")
    ]
    if not lines:
        return "no_hits"
    first_field = lines[0].split("\t", maxsplit=1)[0].strip().lower()
    known_headers = {
        "protein identifier",
        "orf_id",
        "read_id",
        "query",
        "sequence",
    }
    if first_field in known_headers:
        return "no_hits" if len(lines) == 1 else "completed"
    return "completed"


def main() -> int:
    """Write the manifest TSV."""
    args = parse_args()
    with args.output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=HEADERS, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(build_rows(args))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
