#!/usr/bin/env python3
"""Validate read-only database inputs for bacterial AMR analysis."""

from __future__ import annotations

import argparse
import csv
import os
import re
import sys
import tarfile
from pathlib import Path
from typing import Iterable

VERSION = "0.2.0"
VALIDATION_HEADERS = ("tool", "path", "status", "details")


class DatabaseValidationError(ValueError):
    """Raised when a configured AMR database cannot satisfy a tool contract."""


def _path_text(value: object) -> str:
    return str(value or "").strip()


def _readable_directory(path: Path, label: str) -> Path:
    resolved = path.expanduser().resolve()
    if not resolved.is_dir():
        raise DatabaseValidationError(f"{label} is not a directory: {resolved}")
    if not os.access(resolved, os.R_OK | os.X_OK):
        raise DatabaseValidationError(f"{label} is not readable: {resolved}")
    return resolved


def _required_file(path: Path, label: str) -> None:
    if not path.is_file() or not os.access(path, os.R_OK):
        raise DatabaseValidationError(f"{label} is missing or unreadable: {path}")


def _amrfinder_members(path: Path) -> tuple[set[str], str, str]:
    """Return member basenames, database version, and format version."""
    if path.is_dir():
        members = {item.name for item in path.iterdir() if item.is_file()}
        version_path = path / "version.txt"
        format_path = path / "database_format_version.txt"
        _required_file(version_path, "AMRFinderPlus version file")
        _required_file(format_path, "AMRFinderPlus database format file")
        version = version_path.read_text(encoding="utf-8").strip()
        database_format = format_path.read_text(encoding="utf-8").strip()
        return members, version, database_format

    if not path.is_file() or not tarfile.is_tarfile(path):
        raise DatabaseValidationError(
            "AMRFinderPlus database must be a readable directory or tar archive: "
            f"{path}"
        )
    with tarfile.open(path, "r:*") as archive:
        files = [member for member in archive.getmembers() if member.isfile()]
        members = {Path(member.name).name for member in files}
        text_members: dict[str, str] = {}
        for required_name in ("version.txt", "database_format_version.txt"):
            matching = [
                member
                for member in files
                if Path(member.name).name == required_name
            ]
            if not matching:
                raise DatabaseValidationError(
                    f"AMRFinderPlus archive has no {required_name}: {path}"
                )
            extracted = archive.extractfile(matching[0])
            if extracted is None:
                raise DatabaseValidationError(
                    f"AMRFinderPlus archive {required_name} is unreadable: {path}"
                )
            text_members[required_name] = (
                extracted.read().decode("utf-8").strip()
            )
    return (
        members,
        text_members["version.txt"],
        text_members["database_format_version.txt"],
    )


def validate_amrfinder(path: str) -> tuple[Path, str, str, list[str]]:
    """Validate AMRFinderPlus content and list mutation-capable taxgroups."""
    path_text = _path_text(path)
    if not path_text:
        return Path("not_configured"), "skipped", "skipped", []
    resolved = Path(path_text).expanduser().resolve()
    members, version, database_format = _amrfinder_members(resolved)
    for required in ("AMR.LIB", "version.txt", "database_format_version.txt"):
        if required not in members:
            raise DatabaseValidationError(
                f"AMRFinderPlus database is missing required file '{required}': {resolved}"
            )
    if not database_format.startswith("3.12."):
        raise DatabaseValidationError(
            "The local AMRFinderPlus 3.12.8 wrapper requires "
            "database_format_version 3.12.x, but found "
            f"'{database_format or 'unknown'}' in {resolved}"
        )
    organisms = sorted(
        {
            match.group(1)
            for name in members
            if (match := re.fullmatch(r"AMR_DNA-(.+?)(?:\.[A-Za-z0-9]+)?", name))
        }
    )
    return (
        resolved,
        version or "unknown",
        database_format or "unknown",
        organisms,
    )


def validate_abricate(path_value: str, database_id: str) -> tuple[str, str]:
    """Validate an optional custom ABRicate data directory."""
    if not path_value:
        return "bundled", f"container database id={database_id}"
    path = _readable_directory(Path(path_value), "ABRicate database directory")
    matching_dirs = [
        candidate
        for candidate in (path / database_id, path / "db" / database_id)
        if candidate.is_dir()
    ]
    if not matching_dirs:
        raise DatabaseValidationError(
            f"ABRicate database id '{database_id}' was not found below {path}"
        )
    return str(path), f"custom database id={database_id}"


def validate_rgi(path: str, include_wildcard: bool) -> tuple[Path, str]:
    """Validate the processed CARD directory expected by rgi/main."""
    path_text = _path_text(path)
    if not path_text:
        return Path("not_configured"), "skipped"
    resolved = _readable_directory(Path(path_text), "RGI CARD database")
    _required_file(resolved / "card.json", "RGI CARD card.json")
    all_models = sorted(resolved.glob("card_database_v*_all.fasta"))
    if not all_models:
        raise DatabaseValidationError(
            f"RGI CARD database has no card_database_v*_all.fasta: {resolved}"
        )
    match = re.search(r"card_database_v(.+?)_all\.fasta$", all_models[0].name)
    version = match.group(1) if match else "unknown"
    if include_wildcard:
        for required in (
            resolved / f"wildcard_database_v{version}.fasta",
            resolved / f"wildcard_database_v{version}_all.fasta",
            resolved / "wildcard" / "index-for-model-sequences.txt",
            resolved / "wildcard" / "all_amr_61mers.txt",
            resolved / "wildcard" / "61_kmer_db.json",
        ):
            _required_file(required, "RGI WildCARD database file")
    return resolved, version


def validate_resfinder(path: str) -> Path:
    """Validate the acquired-gene ResFinder database."""
    path_text = _path_text(path)
    if not path_text:
        return Path("not_configured")
    resolved = _readable_directory(Path(path_text), "ResFinder database")
    _required_file(resolved / "config", "ResFinder config")
    if not any(resolved.glob("*.fsa")):
        raise DatabaseValidationError(
            f"ResFinder database has no top-level *.fsa files: {resolved}"
        )
    return resolved


def validate_pointfinder(path: str) -> tuple[Path, list[str]]:
    """Validate PointFinder and return supported database directory names."""
    path_text = _path_text(path)
    if not path_text:
        return Path("not_configured"), []
    resolved = _readable_directory(Path(path_text), "PointFinder database")
    _required_file(resolved / "config", "PointFinder config")
    supported = sorted(
        directory.name
        for directory in resolved.iterdir()
        if directory.is_dir()
        and (directory / "genes.txt").is_file()
        and any(directory.glob("*.fsa"))
    )
    if not supported:
        raise DatabaseValidationError(
            f"PointFinder database has no usable species directories: {resolved}"
        )
    return resolved, supported


def validate_deeparg(path: str) -> Path:
    """Validate the DeepARG model/database bundle used by deeparg/predict."""
    path_text = _path_text(path)
    if not path_text:
        return Path("not_configured")
    resolved = _readable_directory(Path(path_text), "DeepARG database")
    database_dirs = [
        directory
        for directory in resolved.rglob("database")
        if directory.is_dir()
    ]
    if not database_dirs:
        raise DatabaseValidationError(
            f"DeepARG bundle contains no directory named 'database': {resolved}"
        )
    if not any((directory.parent / "model").is_dir() for directory in database_dirs):
        raise DatabaseValidationError(
            f"DeepARG bundle contains no model directory beside its database: {resolved}"
        )
    return resolved


def write_rows(path: Path, headers: Iterable[str], rows: Iterable[dict[str, str]]) -> None:
    """Write a deterministic TSV file."""
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(headers), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def validate_all(args: argparse.Namespace) -> None:
    """Validate all configured database contracts and write routing metadata."""
    (
        amrfinder_path,
        amrfinder_version,
        amrfinder_format,
        amrfinder_organisms,
    ) = validate_amrfinder(args.amrfinderplus_db)
    abricate_path, abricate_details = validate_abricate(
        args.abricate_db, args.abricate_db_id
    )
    rgi_path, rgi_version = validate_rgi(args.rgi_db, args.include_rgi_wildcard)
    resfinder_path = validate_resfinder(args.resfinder_db)
    pointfinder_path, pointfinder_species = validate_pointfinder(args.pointfinder_db)
    deeparg_path = validate_deeparg(args.deeparg_db)

    rows = [
        {
            "tool": "AMRFinderPlus",
            "path": str(amrfinder_path),
            "status": (
                "skipped"
                if str(amrfinder_path) == "not_configured"
                else "validated"
            ),
            "details": (
                "no database supplied; AMRFinderPlus skipped because this wrapper "
                "does not automatically set up a default database"
                if str(amrfinder_path) == "not_configured"
                else (
                    f"database_version={amrfinder_version};"
                    f"database_format={amrfinder_format};"
                    "tool_version=3.12.8"
                )
            ),
        },
        {
            "tool": "ABRicate",
            "path": abricate_path,
            "status": "validated",
            "details": abricate_details,
        },
        {
            "tool": "RGI",
            "path": str(rgi_path),
            "status": "skipped" if str(rgi_path) == "not_configured" else "validated",
            "details": (
                "no CARD database supplied; RGI skipped because this workflow "
                "does not automatically prepare CARD"
                if str(rgi_path) == "not_configured"
                else (
                    f"database_version={rgi_version};wildcard="
                    f"{str(args.include_rgi_wildcard).lower()}"
                )
            ),
        },
        {
            "tool": "ResFinder",
            "path": str(resfinder_path),
            "status": (
                "skipped" if str(resfinder_path) == "not_configured" else "validated"
            ),
            "details": (
                "no database supplied; ResFinder skipped because this wrapper "
                "does not automatically set up a default database"
                if str(resfinder_path) == "not_configured"
                else "acquired-gene database"
            ),
        },
        {
            "tool": "PointFinder",
            "path": str(pointfinder_path),
            "status": (
                "skipped" if str(pointfinder_path) == "not_configured" else "validated"
            ),
            "details": (
                "no PointFinder database supplied; mutation routing disabled"
                if str(pointfinder_path) == "not_configured"
                else f"supported_taxa={len(pointfinder_species)}"
            ),
        },
        {
            "tool": "DeepARG",
            "path": str(deeparg_path),
            "status": (
                "skipped" if str(deeparg_path) == "not_configured" else "validated"
            ),
            "details": (
                "no database supplied; DeepARG skipped because this workflow "
                "does not download databases automatically"
                if str(deeparg_path) == "not_configured"
                else f"database_version={args.deeparg_db_version}"
            ),
        },
    ]
    write_rows(args.validation_output, VALIDATION_HEADERS, rows)
    write_rows(
        args.amrfinder_organisms_output,
        ("organism",),
        ({"organism": organism} for organism in amrfinder_organisms),
    )
    write_rows(
        args.pointfinder_species_output,
        ("database_name",),
        ({"database_name": species} for species in pointfinder_species),
    )


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--amrfinderplus-db", default="")
    parser.add_argument("--abricate-db", default="")
    parser.add_argument("--abricate-db-id", required=True)
    parser.add_argument("--rgi-db", default="")
    parser.add_argument("--include-rgi-wildcard", action="store_true")
    parser.add_argument("--resfinder-db", default="")
    parser.add_argument("--pointfinder-db", default="")
    parser.add_argument("--deeparg-db", default="")
    parser.add_argument("--deeparg-db-version", required=True)
    parser.add_argument("--validation-output", required=True, type=Path)
    parser.add_argument("--amrfinder-organisms-output", required=True, type=Path)
    parser.add_argument("--pointfinder-species-output", required=True, type=Path)
    parser.add_argument("--version", action="version", version=VERSION)
    return parser.parse_args()


def main() -> int:
    """Run the database validation program."""
    args = parse_args()
    try:
        validate_all(args)
    except (OSError, tarfile.TarError, DatabaseValidationError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
