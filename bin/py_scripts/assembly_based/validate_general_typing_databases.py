#!/usr/bin/env python3
"""Validate explicit databases used by bacterial general typing."""

from __future__ import annotations

import argparse
import csv
import math
import os
import re
from pathlib import Path
from typing import Iterable


VERSION = "1.1.0"
FASTA_SEQUENCE = re.compile(r"^[A-Za-z*.-]+$")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--virulencefinder-db")
    parser.add_argument("--blast-reference", default="")
    parser.add_argument("--blast-db", default="")
    parser.add_argument("--blast-db-label")
    parser.add_argument("--blast-min-identity", type=float)
    parser.add_argument("--blast-min-coverage", type=float)
    parser.add_argument("--validation-output")
    parser.add_argument("--version", action="version", version=VERSION)
    return parser.parse_args()


def require_readable(path: Path, kind: str) -> None:
    if not path.exists():
        raise ValueError(f"{kind} does not exist: {path}")
    if not os.access(path, os.R_OK):
        raise ValueError(f"{kind} is not readable: {path}")


def nonempty_lines(path: Path) -> Iterable[str]:
    with path.open(encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if line:
                yield line


def validate_virulencefinder_db(path_text: str) -> dict[str, str]:
    if not (path_text or "").strip():
        return {
            "component": "virulencefinder",
            "status": "skipped",
            "path": "",
            "details": (
                "no VirulenceFinder database supplied and no automatic default "
                "database setup is configured"
            ),
        }
    path = Path(path_text).expanduser().resolve()
    require_readable(path, "VirulenceFinder database directory")
    if not path.is_dir():
        raise ValueError(f"VirulenceFinder database is not a directory: {path}")

    config = path / "config"
    notes = path / "notes.txt"
    require_readable(config, "VirulenceFinder config")
    require_readable(notes, "VirulenceFinder notes.txt")

    extensions: list[str] = []
    database_rows: list[tuple[str, str, str]] = []
    for line in nonempty_lines(config):
        if line.startswith("#"):
            if "extensions:" in line:
                extensions = [
                    value.strip().lstrip(".")
                    for value in line.split("extensions:", 1)[1].split(",")
                    if value.strip().lstrip(".")
                ]
            continue
        fields = line.split("\t")
        if len(fields) != 3 or not all(field.strip() for field in fields):
            raise ValueError(
                "Every VirulenceFinder config data row must contain three "
                f"non-empty tab-separated fields: {line!r}"
            )
        database_rows.append(tuple(field.strip() for field in fields))

    if not extensions:
        raise ValueError(
            "VirulenceFinder config must declare at least one '# ... extensions:' value"
        )
    if not database_rows:
        raise ValueError("VirulenceFinder config contains no database rows")

    for prefix, _taxon, _description in database_rows:
        for extension in extensions:
            database_file = path / f"{prefix}.{extension}"
            require_readable(database_file, "VirulenceFinder database file")
            if not database_file.is_file() or database_file.stat().st_size == 0:
                raise ValueError(
                    f"VirulenceFinder database file is empty or not a file: {database_file}"
                )

    note_count = 0
    for line in nonempty_lines(notes):
        note_count += 1
        key, separator, _value = line.partition(":")
        if not separator or not key.strip():
            raise ValueError(
                "Every non-empty VirulenceFinder notes.txt row must contain "
                f"a non-empty key followed by ':': {line!r}"
            )
    if note_count == 0:
        raise ValueError("VirulenceFinder notes.txt contains no entries")

    return {
        "component": "virulencefinder",
        "status": "validated",
        "path": str(path),
        "details": (
            f"{len(database_rows)} database prefix(es); extensions="
            f"{','.join(extensions)}; {note_count} note row(s)"
        ),
    }


def validate_reference_fasta(path_text: str) -> dict[str, str]:
    path = Path(path_text).expanduser().resolve()
    require_readable(path, "Bacterial BLAST reference FASTA")
    if not path.is_file():
        raise ValueError(f"Bacterial BLAST reference is not a file: {path}")

    record_count = 0
    sequence_seen = False
    with path.open(encoding="utf-8") as handle:
        numbered_lines = enumerate(handle, start=1)
        lines = (
            (line_number, raw_line.strip())
            for line_number, raw_line in numbered_lines
            if raw_line.strip()
        )
        for line_number, line in lines:
            if line.startswith(">"):
                if record_count and not sequence_seen:
                    raise ValueError(
                        f"FASTA record before line {line_number} has no sequence in {path}"
                    )
                if len(line) == 1:
                    raise ValueError(
                        f"Empty FASTA identifier in bacterial BLAST reference: {path}"
                    )
                record_count += 1
                sequence_seen = False
                continue
            if record_count == 0:
                raise ValueError(
                    f"Sequence precedes the first FASTA header in {path} (line {line_number})"
                )
            if not FASTA_SEQUENCE.fullmatch(line):
                raise ValueError(
                    f"Invalid FASTA sequence characters in {path} (line {line_number})"
                )
            sequence_seen = True

    if record_count == 0 or not sequence_seen:
        raise ValueError(f"Bacterial BLAST reference contains no complete FASTA record: {path}")

    return {
        "component": "blastn_reference",
        "status": "validated",
        "path": str(path),
        "details": f"user-supplied reference FASTA; {record_count} record(s)",
    }


def validate_prebuilt_blast_db(path_text: str) -> dict[str, str]:
    path = Path(path_text).expanduser().resolve()
    require_readable(path, "Prebuilt bacterial BLAST database directory")
    if not path.is_dir():
        raise ValueError(f"Prebuilt bacterial BLAST database is not a directory: {path}")

    aliases = sorted(path.rglob("*.nal"))
    indexes = sorted(path.rglob("*.nin"))
    if len(aliases) > 1:
        raise ValueError(
            "Prebuilt bacterial BLAST database must expose one logical database; "
            f"found {len(aliases)} .nal aliases"
        )
    if not aliases and len(indexes) != 1:
        raise ValueError(
            "Prebuilt bacterial BLAST database must contain exactly one .nin index "
            "when no .nal alias is present"
        )

    if aliases:
        require_readable(aliases[0], "BLAST nucleotide alias")
        if aliases[0].stat().st_size == 0:
            raise ValueError(f"BLAST nucleotide alias is empty: {aliases[0]}")
        logical_db = aliases[0]
    else:
        logical_db = indexes[0]
        stem = logical_db.with_suffix("")
        for extension in (".nhr", ".nsq"):
            companion = Path(f"{stem}{extension}")
            require_readable(companion, "BLAST nucleotide index companion")
            if companion.stat().st_size == 0:
                raise ValueError(f"BLAST nucleotide index companion is empty: {companion}")

    return {
        "component": "blastn_reference",
        "status": "validated",
        "path": str(path),
        "details": f"user-supplied prebuilt nucleotide database; logical index={logical_db.name}",
    }


def validate_contract(args: argparse.Namespace) -> list[dict[str, str]]:
    reference = (args.blast_reference or "").strip()
    database = (args.blast_db or "").strip()
    label = (args.blast_db_label or "").strip()
    minimum_identity = args.blast_min_identity
    minimum_coverage = args.blast_min_coverage

    if reference and database:
        raise ValueError(
            "Provide at most one of --bacterial_blast_reference or --bacterial_blast_db"
        )
    if (reference or database) and not label:
        raise ValueError(
            "--bacterial_blast_db_label is required to identify the user-supplied BLAST reference"
        )
    if reference or database:
        for value, option in (
            (minimum_identity, "--bacterial_blast_min_identity"),
            (minimum_coverage, "--bacterial_blast_min_coverage"),
        ):
            if value is None:
                raise ValueError(
                    f"{option} is required when a bacterial BLAST reference is configured"
                )
            if not math.isfinite(value) or not 0 <= value <= 100:
                raise ValueError(f"{option} must be a finite number between 0 and 100")

    rows = [validate_virulencefinder_db(args.virulencefinder_db)]
    if reference or database:
        rows.append(
            validate_reference_fasta(reference)
            if reference
            else validate_prebuilt_blast_db(database)
        )
        rows.append(
            {
                "component": "blastn_interpretation",
                "status": "validated",
                "path": label,
                "details": (
                    "sequence similarity to the explicitly user-supplied reference only; "
                    f"minimum identity={minimum_identity:g}%; minimum subject "
                    f"coverage={minimum_coverage:g}%; no cross-tool or biological "
                    "conclusion is assigned"
                ),
            }
        )
    else:
        rows.append(
            {
                "component": "blastn_reference",
                "status": "skipped",
                "path": "",
                "details": "no BLAST reference FASTA or prebuilt database supplied",
            }
        )
        rows.append(
            {
                "component": "blastn_interpretation",
                "status": "skipped",
                "path": "",
                "details": "BLASTN skipped; no sequence-similarity interpretation assigned",
            }
        )
    return rows


def write_rows(path: Path, rows: list[dict[str, str]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["component", "status", "path", "details"],
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    args = parse_args()
    rows = validate_contract(args)
    write_rows(Path(args.validation_output), rows)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
