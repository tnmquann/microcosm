"""Targeted tests for the per-genome AMR execution manifest."""

from __future__ import annotations

import argparse
import importlib.util
import tempfile
import unittest
from pathlib import Path

SCRIPT_PATH = Path(__file__).resolve().parents[1] / "write_amr_manifest.py"
SPEC = importlib.util.spec_from_file_location("write_amr_manifest", SCRIPT_PATH)
if SPEC is None or SPEC.loader is None:
    raise RuntimeError("Cannot load write_amr_manifest.py")
manifest = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(manifest)


class AmrManifestTests(unittest.TestCase):
    """Check supported, unsupported, and skipped status reporting."""

    def args(self, pointfinder_species: str, normalizable: bool) -> argparse.Namespace:
        """Build manifest arguments."""
        return argparse.Namespace(
            sample_id="sample",
            genome_id="genome",
            source_route="standalone",
            amrfinder_report=Path("amrfinder.tsv"),
            abricate_report=Path("abricate.txt"),
            rgi_report=Path("rgi.txt"),
            resfinder_report=Path("resfinder.json"),
            deeparg_report=Path("deeparg.mapping.ARG"),
            pointfinder_species=pointfinder_species,
            abricate_normalizable=normalizable,
        )

    def test_pointfinder_skip_and_unsupported_normalization(self) -> None:
        rows = manifest.build_rows(self.args("", False))
        by_tool = {row["tool"]: row for row in rows}
        self.assertEqual(by_tool["PointFinder"]["execution_status"], "skipped")
        self.assertEqual(by_tool["ABRicate"]["normalization_status"], "not_supported")
        self.assertEqual(by_tool["ResFinder"]["harmonization_status"], "not_supported")
        self.assertNotIn("skipped", by_tool["AMRFinderPlus"]["notes"].lower())

    def test_pointfinder_activation_and_supported_abricate(self) -> None:
        rows = manifest.build_rows(self.args("Escherichia coli", True))
        by_tool = {row["tool"]: row for row in rows}
        self.assertEqual(by_tool["PointFinder"]["execution_status"], "completed")
        self.assertIn("Escherichia coli", by_tool["PointFinder"]["notes"])
        self.assertEqual(by_tool["ABRicate"]["normalization_status"], "supported")

    def test_missing_deeparg_database_is_reported_as_skipped(self) -> None:
        with tempfile.TemporaryDirectory() as tempdir:
            args = self.args("", False)
            args.deeparg_report = Path(tempdir) / "genome.deeparg.skipped.ARG"
            args.deeparg_report.write_text("DeepARG skipped: no database\n")
            rows = manifest.build_rows(args)
            deeparg = {row["tool"]: row for row in rows}["DeepARG"]
            self.assertEqual(deeparg["execution_status"], "skipped")
            self.assertEqual(deeparg["harmonization_status"], "not_run")
            self.assertEqual(deeparg["normalization_status"], "not_run")

    def test_no_hits_and_database_provenance_are_explicit(self) -> None:
        with tempfile.TemporaryDirectory() as tempdir:
            root = Path(tempdir)
            args = self.args("Escherichia coli", True)
            for attribute, filename in (
                ("amrfinder_report", "amrfinder.tsv"),
                ("abricate_report", "abricate.txt"),
                ("rgi_report", "rgi.txt"),
                ("deeparg_report", "deeparg.mapping.ARG"),
            ):
                report = root / filename
                report.write_text("", encoding="utf-8")
                setattr(args, attribute, report)
            args.resfinder_report = root / "resfinder.json"
            args.resfinder_report.write_text("{}\n", encoding="utf-8")
            args.database_validation = root / "amr_database_validation.tsv"
            args.database_validation.write_text(
                "tool\tpath\tstatus\tdetails\n"
                "AMRFinderPlus\t/db/amrfinder\tvalidated\tdatabase_version=1\n",
                encoding="utf-8",
            )

            rows = manifest.build_rows(args)
            by_tool = {row["tool"]: row for row in rows}
            self.assertEqual(by_tool["AMRFinderPlus"]["request_status"], "requested")
            self.assertEqual(by_tool["AMRFinderPlus"]["result_status"], "no_hits")
            self.assertEqual(
                by_tool["AMRFinderPlus"]["database_path"], "/db/amrfinder"
            )
            self.assertEqual(by_tool["ResFinder"]["result_status"], "no_hits")


if __name__ == "__main__":
    unittest.main()
