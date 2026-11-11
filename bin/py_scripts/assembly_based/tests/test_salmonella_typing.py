"""Focused tests for exact Salmonella Typhi/Paratyphi routing evidence."""

from __future__ import annotations

import csv
import importlib.util
import json
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

    def reports(self, serovar: str, mykrobe_species: dict[str, object], sistr_serovar: str = "") -> dict[str, str]:
        seqsero2 = self.root / "seqsero2.tsv"
        with seqsero2.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(
                handle,
                fieldnames=["Predicted serotype", "Note"],
                delimiter="\t",
            )
            writer.writeheader()
            writer.writerow({"Predicted serotype": serovar, "Note": "test"})
        mykrobe = self.root / "mykrobe.json"
        mykrobe.write_text(json.dumps({"species": mykrobe_species}), encoding="utf-8")
        sistr = self.root / "sistr.tsv"
        with sistr.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=["serovar", "serogroup"], delimiter="\t")
            writer.writeheader()
            writer.writerow({"serovar": sistr_serovar, "serogroup": "test"})
        _evidence, routing = parser.build_reports("sample", "genome", seqsero2, mykrobe, sistr)
        return routing

    def test_seqsero2_typhi_routes_genotyphi(self) -> None:
        routing = self.reports("Typhi", {})
        self.assertEqual(routing["route_genotyphi"], "true")
        self.assertEqual(routing["route_paratype"], "false")

    def test_mykrobe_typhi_coverage_routes_genotyphi(self) -> None:
        routing = self.reports("Typhimurium", {"salmonella_typhi": {"percent_coverage": 99.1}})
        self.assertEqual(routing["route_genotyphi"], "true")
        self.assertEqual(routing["tool_conflict"], "true")

    def test_typhimurium_never_routes_genotyphi(self) -> None:
        routing = self.reports("Typhimurium", {})
        self.assertEqual(routing["route_genotyphi"], "false")
        self.assertEqual(routing["route_paratype"], "false")

    def test_only_paratyphi_a_routes_paratype(self) -> None:
        paratyphi_a = self.reports("Paratyphi A", {})
        paratyphi_b = self.reports("Paratyphi B", {})
        paratyphi_c = self.reports("Paratyphi C", {})
        self.assertEqual(paratyphi_a["route_paratype"], "true")
        self.assertEqual(paratyphi_b["route_paratype"], "false")
        self.assertEqual(paratyphi_c["route_paratype"], "false")

    def test_sistr_evidence_is_preserved_but_does_not_route_follow_up(self) -> None:
        seqsero2 = self.root / "seqsero2.tsv"
        seqsero2.write_text("Predicted serotype\tNote\nTyphimurium\ttest\n", encoding="utf-8")
        mykrobe = self.root / "mykrobe.json"
        mykrobe.write_text(json.dumps({"species": {}}), encoding="utf-8")
        sistr = self.root / "sistr.tsv"
        sistr.write_text("serovar\tserogroup\nTyphi\tD1\n", encoding="utf-8")
        evidence, routing = parser.build_reports("sample", "genome", seqsero2, mykrobe, sistr)
        self.assertEqual([row["tool"] for row in evidence], ["seqsero2", "sistr", "mykrobe"])
        self.assertEqual(evidence[1]["normalized_prediction"], "Typhi")
        self.assertEqual(routing["route_genotyphi"], "false")
        self.assertEqual(routing["tool_conflict"], "true")
        self.assertIn("tool_conflict_preserved", routing["routing_reason"])

    def test_structured_inputs_reject_duplicate_keys_and_columns(self) -> None:
        duplicate_json = self.root / "duplicate.json"
        duplicate_json.write_text(
            '{"species": {"typhi": {"percent_coverage": 1}, '
            '"typhi": {"percent_coverage": 2}}}',
            encoding="utf-8",
        )
        with self.assertRaisesRegex(ValueError, "duplicate key"):
            parser.parse_mykrobe(duplicate_json)

        duplicate_tsv = self.root / "duplicate.tsv"
        duplicate_tsv.write_text(
            "Predicted serotype\tPredicted serotype\nTyphi\tTyphimurium\n",
            encoding="utf-8",
        )
        with self.assertRaisesRegex(ValueError, "duplicate column"):
            parser.parse_seqsero2(duplicate_tsv)

    def test_nonfinite_coverage_is_not_positive(self) -> None:
        for value in ("NaN", "Infinity", "-Infinity"):
            with self.subTest(value=value):
                self.assertIsNone(parser.numeric_coverage(value))

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
