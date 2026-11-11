"""Regression tests for deterministic fixed-schema BLASTN hit filtering."""

from __future__ import annotations

import csv
import importlib.util
import tempfile
import unittest
from pathlib import Path


SCRIPT = Path(__file__).resolve().parents[1] / "filter_blasthits.py"
SPEC = importlib.util.spec_from_file_location("filter_blasthits", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
filters = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(filters)


class FilterBlastHitsTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary_directory = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary_directory.name)

    def tearDown(self) -> None:
        self.temporary_directory.cleanup()

    @staticmethod
    def hit(
        query: str,
        subject: str,
        identity: str,
        start: str,
        end: str,
        score: str,
        length: str = "100",
        subject_length: str = "100",
    ) -> list[str]:
        return [
            query,
            subject,
            identity,
            length,
            "0",
            "0",
            start,
            end,
            "1",
            length,
            "1e-20",
            score,
            subject_length,
        ]

    def write_hits(self, name: str, rows: list[list[str]]) -> Path:
        path = self.root / name
        with path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerows(rows)
        return path

    def test_thresholds_and_overlap_selection_are_order_independent(self) -> None:
        lower = self.hit("query", "gene-b", "91", "100", "1", "40")
        higher = self.hit("query", "gene-b", "99", "1", "100", "60")
        unrelated = self.hit("query", "gene-a", "95", "200", "299", "50")
        low_identity = self.hit("query", "gene-c", "89.9", "400", "499", "70")
        low_coverage = self.hit(
            "query", "gene-d", "99", "500", "549", "80", length="50", subject_length="100"
        )

        first_rows, first_counts = filters.filter_results(
            self.write_hits(
                "first.tsv", [lower, unrelated, low_identity, higher, low_coverage]
            ),
            90,
            80,
            1,
        )
        second_rows, second_counts = filters.filter_results(
            self.write_hits(
                "second.tsv", [low_coverage, higher, low_identity, unrelated, lower]
            ),
            90,
            80,
            1,
        )
        self.assertEqual(first_rows, second_rows)
        self.assertEqual(first_counts, second_counts)
        self.assertEqual(first_counts, {"gene-a": 1, "gene-b": 1})
        self.assertEqual([row[2] for row in first_rows if row[1] == "gene-b"], ["99"])
        self.assertTrue(all(row[-1] == "100" for row in first_rows))

        output_prefix = self.root / "result"
        filters.write_outputs(first_rows, first_counts, output_prefix)
        with Path(f"{output_prefix}.blastn_filtered_hits.tsv").open(
            encoding="utf-8", newline=""
        ) as handle:
            written = list(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual([row["sseqid"] for row in written], ["gene-a", "gene-b"])
        summary = Path(f"{output_prefix}.blastn_summary.tsv").read_text(
            encoding="utf-8"
        )
        self.assertLess(summary.index("gene-a"), summary.index("gene-b"))

    def test_malformed_rows_and_invalid_thresholds_fail(self) -> None:
        malformed = self.write_hits("malformed.tsv", [["too", "short"]])
        with self.assertRaisesRegex(ValueError, "expected 13"):
            filters.filter_results(malformed, 90, 80)

        nonfinite = self.hit("query", "gene-a", "NaN", "1", "100", "50")
        with self.assertRaisesRegex(ValueError, "non-finite"):
            filters.filter_results(self.write_hits("nonfinite.tsv", [nonfinite]), 90, 80)

        empty = self.write_hits("empty.tsv", [])
        for identity, coverage, message in (
            (101, 80, "minimum identity"),
            (90, -1, "minimum coverage"),
        ):
            with self.subTest(identity=identity, coverage=coverage):
                with self.assertRaisesRegex(ValueError, message):
                    filters.filter_results(empty, identity, coverage)
        with self.assertRaisesRegex(ValueError, "minimum overlap"):
            filters.filter_results(empty, 90, 80, 0)

    def test_empty_input_writes_schema_only_outputs(self) -> None:
        rows, counts = filters.filter_results(self.write_hits("empty.tsv", []), 90, 80)
        prefix = self.root / "empty"
        filters.write_outputs(rows, counts, prefix)
        detailed = Path(f"{prefix}.blastn_filtered_hits.tsv").read_text(encoding="utf-8")
        summary = Path(f"{prefix}.blastn_summary.tsv").read_text(encoding="utf-8")
        self.assertEqual(detailed.splitlines()[0].split("\t"), list(filters.FILTERED_COLUMNS))
        self.assertEqual(summary, "gene\tcount\n")


if __name__ == "__main__":
    unittest.main()
