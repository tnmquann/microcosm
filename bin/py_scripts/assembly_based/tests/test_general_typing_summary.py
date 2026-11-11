"""Tests for source-preserving general-typing summaries."""

from __future__ import annotations

import importlib.util
import tempfile
import unittest
from pathlib import Path


SCRIPT = Path(__file__).resolve().parents[1] / "summarize_gene_hits.py"
SPEC = importlib.util.spec_from_file_location("summarize_gene_hits", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
summary = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(summary)


class GeneralTypingSummaryTests(unittest.TestCase):
    def test_evidence_retains_raw_identifiers_and_sources(self) -> None:
        with tempfile.TemporaryDirectory() as tempdir:
            root = Path(tempdir)
            virulence = root / "virulence.tsv"
            plasmid = root / "plasmid.tsv"
            blast = root / "blast.tsv"
            virulence.write_text(
                "Database\tVirulence factor\nVFDB\tGene A\n", encoding="utf-8"
            )
            plasmid.write_text("Database\tPlasmid\nDB\tIncFIB\n", encoding="utf-8")
            blast.write_text("qseqid\tsseqid\ncontig\tblaTEM-1\n", encoding="utf-8")
            rows = summary.build_evidence_rows(
                "sample", "genome", virulence, plasmid, blast
            )
            self.assertEqual(
                {(row["source_tool"], row["normalized_identifier"]) for row in rows},
                {
                    ("virulencefinder", "gene_a"),
                    ("plasmidfinder", "incfib"),
                    ("blastn", "blatem-1"),
                },
            )
            self.assertIn("blaTEM-1", {row["raw_identifier"] for row in rows})

    def test_mlst_preserves_profile_and_uses_explicit_na_for_cc(self) -> None:
        with tempfile.TemporaryDirectory() as tempdir:
            mlst = Path(tempdir) / "sample.mlst.tsv"
            mlst.write_text(
                "sample.fa\tecoli_achtman_4\t131\tadk(53)\tfumC(40)\n",
                encoding="utf-8",
            )
            row = summary.build_mlst_rows("sample", "genome", mlst)[0]
            self.assertEqual(row["scheme"], "ecoli_achtman_4")
            self.assertEqual(row["sequence_type"], "131")
            self.assertEqual(row["allele_profile"], "adk(53);fumC(40)")
            self.assertEqual(row["clonal_complex"], "NA")
            self.assertEqual(row["clonal_complex_lookup_status"], "not_configured")

    def test_no_mlst_match_is_unknown(self) -> None:
        with tempfile.TemporaryDirectory() as tempdir:
            mlst = Path(tempdir) / "sample.mlst.tsv"
            mlst.write_text("sample.fa\t-\t-\n", encoding="utf-8")
            row = summary.build_mlst_rows("sample", "genome", mlst)[0]
            self.assertEqual(row["scheme"], "Unknown")
            self.assertEqual(row["sequence_type"], "Unknown")
            self.assertEqual(row["allele_profile"], "Unknown")


if __name__ == "__main__":
    unittest.main()
