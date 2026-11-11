"""Targeted tests for explicit bacterial AMR database contracts."""

from __future__ import annotations

import argparse
import csv
import importlib.util
import tempfile
import unittest
from pathlib import Path

SCRIPT_PATH = Path(__file__).resolve().parents[1] / "validate_amr_databases.py"
SPEC = importlib.util.spec_from_file_location("validate_amr_databases", SCRIPT_PATH)
if SPEC is None or SPEC.loader is None:
    raise RuntimeError("Cannot load validate_amr_databases.py")
validator = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(validator)


class AmrDatabaseValidationTests(unittest.TestCase):
    """Validate accepted structures and actionable failures."""

    def setUp(self) -> None:
        self.temp_dir = tempfile.TemporaryDirectory()
        self.root = Path(self.temp_dir.name)

        self.amrfinder = self.root / "amrfinder"
        self.amrfinder.mkdir()
        (self.amrfinder / "AMR.LIB").touch()
        (self.amrfinder / "version.txt").write_text("2026-01-01.1\n")
        (self.amrfinder / "database_format_version.txt").write_text("3.12.1\n")
        (self.amrfinder / "AMR_DNA-Escherichia").touch()
        (self.amrfinder / "AMR_DNA-Salmonella.tab").touch()

        self.abricate = self.root / "abricate"
        (self.abricate / "ncbi").mkdir(parents=True)

        self.rgi = self.root / "rgi"
        (self.rgi / "wildcard").mkdir(parents=True)
        (self.rgi / "card.json").touch()
        (self.rgi / "card_database_v4.0.1_all.fasta").touch()
        (self.rgi / "wildcard_database_v4.0.1.fasta").touch()
        (self.rgi / "wildcard_database_v4.0.1_all.fasta").touch()
        for name in (
            "index-for-model-sequences.txt",
            "all_amr_61mers.txt",
            "61_kmer_db.json",
        ):
            (self.rgi / "wildcard" / name).touch()

        self.resfinder = self.root / "resfinder"
        self.resfinder.mkdir()
        (self.resfinder / "config").touch()
        (self.resfinder / "all.fsa").touch()

        self.pointfinder = self.root / "pointfinder"
        (self.pointfinder / "escherichia_coli").mkdir(parents=True)
        (self.pointfinder / "config").touch()
        (self.pointfinder / "escherichia_coli" / "genes.txt").touch()
        (self.pointfinder / "escherichia_coli" / "gyrA.fsa").touch()

        self.deeparg = self.root / "deeparg"
        (self.deeparg / "database").mkdir(parents=True)
        (self.deeparg / "model").mkdir()

    def tearDown(self) -> None:
        self.temp_dir.cleanup()

    def args(self) -> argparse.Namespace:
        """Build a complete validation argument set."""
        return argparse.Namespace(
            amrfinderplus_db=self.amrfinder,
            abricate_db=str(self.abricate),
            abricate_db_id="ncbi",
            rgi_db=self.rgi,
            include_rgi_wildcard=True,
            resfinder_db=self.resfinder,
            pointfinder_db=self.pointfinder,
            deeparg_db=self.deeparg,
            deeparg_db_version="1",
            validation_output=self.root / "validation.tsv",
            amrfinder_organisms_output=self.root / "amrfinder_organisms.tsv",
            pointfinder_species_output=self.root / "pointfinder_species.tsv",
        )

    def test_complete_database_contract(self) -> None:
        args = self.args()
        validator.validate_all(args)
        with args.validation_output.open(encoding="utf-8") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual([row["tool"] for row in rows], [
            "AMRFinderPlus",
            "ABRicate",
            "RGI",
            "ResFinder",
            "PointFinder",
            "DeepARG",
        ])
        self.assertIn("database_version=2026-01-01.1", rows[0]["details"])
        self.assertIn("database_format=3.12.1", rows[0]["details"])
        self.assertIn("tool_version=3.12.8", rows[0]["details"])
        self.assertIn(
            "Escherichia",
            args.amrfinder_organisms_output.read_text(encoding="utf-8"),
        )
        self.assertIn(
            "escherichia_coli",
            args.pointfinder_species_output.read_text(encoding="utf-8"),
        )

    def test_bundled_abricate_database_is_explicit(self) -> None:
        path, details = validator.validate_abricate("", "ncbi")
        self.assertEqual(path, "bundled")
        self.assertIn("database id=ncbi", details)

    def test_missing_optional_databases_are_recorded_without_failing(self) -> None:
        args = self.args()
        args.amrfinderplus_db = ""
        args.rgi_db = ""
        args.resfinder_db = ""
        args.pointfinder_db = ""
        args.deeparg_db = ""
        validator.validate_all(args)
        with args.validation_output.open(encoding="utf-8") as handle:
            rows = {row["tool"]: row for row in csv.DictReader(handle, delimiter="\t")}
        self.assertEqual(rows["AMRFinderPlus"]["status"], "skipped")
        self.assertEqual(rows["RGI"]["status"], "skipped")
        self.assertEqual(rows["ResFinder"]["status"], "skipped")
        self.assertEqual(rows["PointFinder"]["status"], "skipped")
        self.assertEqual(rows["DeepARG"]["status"], "skipped")
        self.assertIn("does not download", rows["DeepARG"]["details"])

    def test_missing_card_file_is_actionable(self) -> None:
        (self.rgi / "card.json").unlink()
        with self.assertRaisesRegex(
            validator.DatabaseValidationError, "RGI CARD card.json"
        ):
            validator.validate_all(self.args())

    def test_incompatible_amrfinder_format_is_actionable(self) -> None:
        (self.amrfinder / "database_format_version.txt").write_text("4.2.0\n")
        with self.assertRaisesRegex(
            validator.DatabaseValidationError,
            "requires database_format_version 3.12.x",
        ):
            validator.validate_all(self.args())

    def test_missing_pointfinder_taxa_is_actionable(self) -> None:
        (self.pointfinder / "escherichia_coli" / "genes.txt").unlink()
        with self.assertRaisesRegex(
            validator.DatabaseValidationError, "no usable species directories"
        ):
            validator.validate_all(self.args())


if __name__ == "__main__":
    unittest.main()
