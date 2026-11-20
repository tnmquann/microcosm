"""Targeted tests for genome-input and classification normalization."""

from __future__ import annotations

import csv
import importlib.util
import tempfile
import unittest
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parents[1]


def load_module(name: str, filename: str):
    """Load a repository script without requiring a Python package."""
    spec = importlib.util.spec_from_file_location(name, SCRIPT_DIR / filename)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Cannot load {filename}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


validator = load_module("validate_genome_input", "validate_genome_input.py")
normalizer = load_module(
    "normalize_genome_classifications", "normalize_genome_classifications.py"
)

BACTERIAL_LINEAGE = (
    "d__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria;"
    "o__Enterobacterales;f__Enterobacteriaceae;g__Escherichia;"
    "s__Escherichia coli"
)


class GenomeInputValidationTests(unittest.TestCase):
    """Exercise every required standalone input-validation behavior."""

    def setUp(self) -> None:
        self.temp_dir = tempfile.TemporaryDirectory()
        self.root = Path(self.temp_dir.name)
        self.fasta_a = self.root / "a.fna"
        self.fasta_b = self.root / "b.fasta.gz"
        self.fasta_a.write_text(">a\nACGT\n", encoding="utf-8")
        self.fasta_b.write_bytes(b"not-compressed-but-readable-for-validation")

    def tearDown(self) -> None:
        self.temp_dir.cleanup()

    def write_csv(self, headers: list[str], rows: list[list[str]]) -> Path:
        """Write an input fixture."""
        path = self.root / "input.csv"
        with path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle)
            writer.writerow(headers)
            writer.writerows(rows)
        return path

    def test_supplied_lineage_preserves_metadata_and_parses_ranks(self) -> None:
        path = self.write_csv(
            ["sample_id", "run_id", "group", "assembler", "fasta", "lineage"],
            [["sample-A", "run-X", "group-Y", "custom", str(self.fasta_a), BACTERIAL_LINEAGE]],
        )
        rows = validator.validate_rows(path, classification_available=False)
        self.assertEqual(rows[0]["sample_id"], "sample-A")
        self.assertEqual(rows[0]["run_id"], "run-X")
        self.assertEqual(rows[0]["group"], "group-Y")
        self.assertEqual(rows[0]["assembler"], "custom")
        self.assertEqual(rows[0]["domain"], "d__Bacteria")
        self.assertEqual(rows[0]["phylum"], "p__Pseudomonadota")
        self.assertEqual(rows[0]["class"], "c__Gammaproteobacteria")
        self.assertEqual(rows[0]["order"], "o__Enterobacterales")
        self.assertEqual(rows[0]["family"], "f__Enterobacteriaceae")
        self.assertEqual(rows[0]["genus"], "g__Escherichia")
        self.assertEqual(rows[0]["species"], "s__Escherichia coli")
        self.assertEqual(rows[0]["lineage_source"], "provided")

    def test_missing_and_mixed_lineage_require_classification_contract(self) -> None:
        path = self.write_csv(
            ["sample_id", "fasta", "lineage"],
            [
                ["provided", str(self.fasta_a), BACTERIAL_LINEAGE],
                ["missing", str(self.fasta_b), ""],
            ],
        )
        with self.assertRaisesRegex(
            validator.ValidationError, "classification is required.*row 3 sample 'missing'"
        ):
            validator.validate_rows(path, classification_available=False)
        rows = validator.validate_rows(path, classification_available=True)
        self.assertEqual(rows[1]["lineage_source"], "classified")
        self.assertEqual(rows[1]["lineage"], "")

    def test_non_bacterial_lineage_is_valid_but_not_bacterial(self) -> None:
        path = self.write_csv(
            ["sample_id", "fasta", "lineage"],
            [["archaea", str(self.fasta_a), "d__Archaea;p__Thermoproteota"]],
        )
        rows = validator.validate_rows(path, classification_available=False)
        self.assertEqual(rows[0]["domain"], "d__Archaea")

    def test_malformed_lineage_is_row_specific(self) -> None:
        path = self.write_csv(
            ["sample_id", "fasta", "lineage"],
            [["bad", str(self.fasta_a), "Bacteria;g__Escherichia"]],
        )
        with self.assertRaisesRegex(
            validator.ValidationError, "row 2 sample 'bad'.*invalid lineage token"
        ):
            validator.validate_rows(path, classification_available=False)

    def test_missing_header_and_missing_fasta_are_actionable(self) -> None:
        path = self.write_csv(["sample_id", "lineage"], [["x", BACTERIAL_LINEAGE]])
        with self.assertRaisesRegex(validator.ValidationError, "missing required header.*fasta"):
            validator.validate_rows(path, classification_available=False)

        path = self.write_csv(
            ["sample_id", "fasta"], [["x", str(self.root / "missing.fa")]]
        )
        with self.assertRaisesRegex(
            validator.ValidationError, "row 2 sample 'x': FASTA does not exist"
        ):
            validator.validate_rows(path, classification_available=True)

    def test_duplicate_sample_and_classification_key_collisions_fail(self) -> None:
        path = self.write_csv(
            ["sample_id", "fasta", "lineage"],
            [
                ["duplicate", str(self.fasta_a), BACTERIAL_LINEAGE],
                ["duplicate", str(self.fasta_b), BACTERIAL_LINEAGE],
            ],
        )
        with self.assertRaisesRegex(validator.ValidationError, "duplicate sample_id"):
            validator.validate_rows(path, classification_available=False)

        other_dir = self.root / "other"
        other_dir.mkdir()
        colliding = other_dir / self.fasta_a.name
        colliding.write_text(">b\nACGT\n", encoding="utf-8")
        path = self.write_csv(
            ["sample_id", "fasta"],
            [["first", str(self.fasta_a)], ["second", str(colliding)]],
        )
        with self.assertRaisesRegex(validator.ValidationError, "FASTA basename 'a' collides"):
            validator.validate_rows(path, classification_available=True)

    def test_extra_csv_values_are_rejected(self) -> None:
        path = self.root / "input.csv"
        path.write_text(
            f"sample_id,fasta\nsample,{self.fasta_a},unexpected\n",
            encoding="utf-8",
        )
        with self.assertRaisesRegex(validator.ValidationError, "more values than the header"):
            validator.validate_rows(path, classification_available=True)


class ClassificationNormalizationTests(unittest.TestCase):
    """Validate exact lineage parsing for existing BIN_CLASSIFICATION outputs."""

    def setUp(self) -> None:
        self.temp_dir = tempfile.TemporaryDirectory()
        self.root = Path(self.temp_dir.name)

    def tearDown(self) -> None:
        self.temp_dir.cleanup()

    def write_classification(self, rows: list[list[str]]) -> Path:
        """Write a Sourmash-like classification CSV."""
        path = self.root / "classification.csv"
        with path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle)
            writer.writerow(["query_name", "status", "lineage"])
            writer.writerows(rows)
        return path

    def test_resolved_and_unresolved_rows(self) -> None:
        path = self.write_classification(
            [
                ["genome-a", "match", BACTERIAL_LINEAGE],
                ["genome-b", "no_match", ""],
            ]
        )
        resolved, unresolved = normalizer.normalize([path])
        self.assertEqual(resolved[0]["domain"], "d__Bacteria")
        self.assertEqual(resolved[0]["phylum"], "p__Pseudomonadota")
        self.assertEqual(resolved[0]["family"], "f__Enterobacteriaceae")
        self.assertEqual(resolved[0]["genus"], "g__Escherichia")
        self.assertEqual(resolved[0]["species"], "s__Escherichia coli")
        self.assertEqual(unresolved[0]["query_name"], "genome-b")
        self.assertEqual(unresolved[0]["detail"], "no_match")

    def test_unresolved_detail_prefers_gtdbtk_warnings(self) -> None:
        path = Path(self.temp_dir.name) / "gtdbtk.csv"
        path.write_text(
            "query_name,status,lineage,note,warnings\n"
            "rand_a,nomatch,,N/A,No bacterial or archaeal marker\n"
            "rand_b,nomatch,,  ,N/A\n",
            encoding="utf-8",
        )
        _, unresolved = normalizer.normalize([path])
        self.assertEqual([row["detail"] for row in unresolved], ["No bacterial or archaeal marker", "nomatch"])

    def test_malformed_and_duplicate_classifications_fail(self) -> None:
        path = self.write_classification(
            [["genome-a", "match", "g__Escherichia;d__Bacteria"]]
        )
        with self.assertRaisesRegex(
            normalizer.ClassificationError, "duplicated or out of canonical order"
        ):
            normalizer.normalize([path])

        path = self.write_classification(
            [
                ["genome-a", "match", BACTERIAL_LINEAGE],
                ["genome-a", "match", BACTERIAL_LINEAGE],
            ]
        )
        with self.assertRaisesRegex(
            normalizer.ClassificationError, "duplicate classification"
        ):
            normalizer.normalize([path])

    def test_duplicate_headers_and_extra_values_fail(self) -> None:
        path = self.root / "classification.csv"
        path.write_text(
            "query_name,lineage,lineage\ngenome-a,d__Bacteria,d__Bacteria\n",
            encoding="utf-8",
        )
        with self.assertRaisesRegex(normalizer.ClassificationError, "duplicate header"):
            normalizer.normalize([path])

        path.write_text(
            "query_name,lineage\ngenome-a,d__Bacteria,unexpected\n",
            encoding="utf-8",
        )
        with self.assertRaisesRegex(normalizer.ClassificationError, "more values"):
            normalizer.normalize([path])


if __name__ == "__main__":
    unittest.main()
