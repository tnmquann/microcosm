"""Tests for SISTR-only Salmonella follow-up routing (playground/modular/4_sal.sh)."""

from __future__ import annotations

import csv
import importlib.util
import tempfile
import unittest
from pathlib import Path


SCRIPT = Path(__file__).resolve().parents[1] / "parse_salmonella_typing.py"
SPEC = importlib.util.spec_from_file_location("parse_salmonella_typing", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
parser = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(parser)


class SalmonellaTypingTests(unittest.TestCase):
    """Verify exact-field parsing and follow-up decisions."""

    def setUp(self) -> None:
        self.temporary_directory = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary_directory.name)

    def tearDown(self) -> None:
        self.temporary_directory.cleanup()

    def reports(self, seqsero2_serovar: str, sistr_serovar: str) -> dict[str, str]:
        seqsero2 = self.root / "seqsero2.tsv"
        seqsero2.write_text(f"Predicted serotype\tNote\n{seqsero2_serovar}\ttest\n", encoding="utf-8")
        sistr = self.root / "sistr.tsv"
        sistr.write_text(f"serovar\tserogroup\n{sistr_serovar}\ttest\n", encoding="utf-8")
        _evidence, routing = parser.build_reports("sample", "genome", seqsero2, sistr)
        return routing

    def test_sistr_typhi_routes_genotyphi_despite_seqsero2_conflict(self) -> None:
        routing = self.reports("Typhimurium", "Typhi")
        self.assertEqual(routing["route_genotyphi"], "true")
        self.assertEqual(routing["route_paratype"], "false")
        self.assertEqual(routing["tool_conflict"], "true")

    def test_sistr_paratyphi_a_routes_paratype_despite_seqsero2_conflict(self) -> None:
        routing = self.reports("Typhi", "Paratyphi A")
        self.assertEqual(routing["route_genotyphi"], "false")
        self.assertEqual(routing["route_paratype"], "true")

    def test_seqsero2_alone_never_routes(self) -> None:
        for seqsero2 in ("Typhi", "Paratyphi A"):
            with self.subTest(seqsero2=seqsero2):
                routing = self.reports(seqsero2, "Enteritidis")
                self.assertEqual(routing["route_genotyphi"], "false")
                self.assertEqual(routing["route_paratype"], "false")

    def test_other_or_missing_sistr_serovar_routes_nothing(self) -> None:
        for sistr in ("Typhimurium", "Paratyphi B", "typhi", ""):
            with self.subTest(sistr=sistr):
                routing = self.reports("Typhi", sistr)
                self.assertEqual(routing["route_genotyphi"], "false")
                self.assertEqual(routing["route_paratype"], "false")
        self.assertEqual(self.reports("Typhi", "")["routing_reason"].split(";")[0], "sistr_serovar_missing")

    def test_evidence_rows_are_seqsero2_then_sistr(self) -> None:
        seqsero2 = self.root / "seqsero2.tsv"
        seqsero2.write_text("Predicted serotype\tNote\nTyphimurium\ttest\n", encoding="utf-8")
        sistr = self.root / "sistr.tsv"
        sistr.write_text("serovar\tserogroup\nTyphi\tD1\n", encoding="utf-8")
        evidence, _routing = parser.build_reports("sample", "genome", seqsero2, sistr)
        self.assertEqual([row["tool"] for row in evidence], ["seqsero2", "sistr"])
        self.assertEqual(evidence[1]["normalized_prediction"], "Typhi")

    def test_duplicate_columns_are_rejected(self) -> None:
        duplicate_tsv = self.root / "duplicate.tsv"
        duplicate_tsv.write_text(
            "Predicted serotype\tPredicted serotype\nTyphi\tTyphimurium\n",
            encoding="utf-8",
        )
        with self.assertRaisesRegex(ValueError, "duplicate column"):
            parser.parse_seqsero2(duplicate_tsv)

    def test_bundle_validator_requires_every_upstream_resource(self) -> None:
        validator_spec = importlib.util.spec_from_file_location(
            "validate_paratype_bundle",
            Path(__file__).resolve().parents[1] / "validate_paratype_bundle.py",
        )
        assert validator_spec is not None and validator_spec.loader is not None
        validator = importlib.util.module_from_spec(validator_spec)
        validator_spec.loader.exec_module(validator)
        bundle = self.root / "bundle"
        bundle.mkdir()
        with self.assertRaisesRegex(ValueError, "paratype.py"):
            validator.validate_bundle(bundle)
        for filename in validator.REQUIRED_FILES:
            (bundle / filename).write_text("test\n", encoding="utf-8")
        self.assertEqual(validator.validate_bundle(bundle), bundle.resolve())


if __name__ == "__main__":
    unittest.main()
