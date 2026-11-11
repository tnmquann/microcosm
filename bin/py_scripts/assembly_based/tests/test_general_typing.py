from __future__ import annotations

import csv
import importlib.util
import tempfile
import unittest
from argparse import Namespace
from pathlib import Path


SCRIPT_ROOT = Path(__file__).resolve().parents[1]


def load_module(name: str, filename: str):
    spec = importlib.util.spec_from_file_location(name, SCRIPT_ROOT / filename)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


validator = load_module(
    "validate_general_typing_databases", "validate_general_typing_databases.py"
)
manifest = load_module(
    "write_general_typing_manifest", "write_general_typing_manifest.py"
)


class GeneralTypingDatabaseValidationTests(unittest.TestCase):
    def setUp(self) -> None:
        self.tempdir = tempfile.TemporaryDirectory()
        self.root = Path(self.tempdir.name)
        self.virulence = self.root / "virulence"
        self.virulence.mkdir()
        (self.virulence / "config").write_text(
            "# extensions: fsa\nvirulence_test\tTest taxon\tTest database\n",
            encoding="utf-8",
        )
        (self.virulence / "notes.txt").write_text(
            "test_gene: test function\n", encoding="utf-8"
        )
        (self.virulence / "virulence_test.fsa").write_text(
            ">test_gene:1:TEST001\nACGT\n", encoding="utf-8"
        )
        self.reference = self.root / "reference.fa"
        self.reference.write_text(">reference\nACGTACGT\n", encoding="utf-8")

    def tearDown(self) -> None:
        self.tempdir.cleanup()

    def args(self, reference: str = "", database: str = "") -> Namespace:
        return Namespace(
            virulencefinder_db=str(self.virulence),
            blast_reference=reference,
            blast_db=database,
            blast_db_label="unit_test_reference",
            blast_min_identity=90.0,
            blast_min_coverage=80.0,
        )

    def test_reference_contract_validates(self) -> None:
        rows = validator.validate_contract(self.args(reference=str(self.reference)))
        self.assertEqual([row["component"] for row in rows], [
            "virulencefinder",
            "blastn_reference",
            "blastn_interpretation",
        ])
        self.assertIn("reference FASTA", rows[1]["details"])

    def test_reference_rejects_an_incomplete_nonfinal_record(self) -> None:
        self.reference.write_text(
            ">missing_sequence\n>complete\nACGT\n", encoding="utf-8"
        )
        with self.assertRaisesRegex(ValueError, "has no sequence"):
            validator.validate_contract(self.args(reference=str(self.reference)))

    def test_prebuilt_database_contract_validates(self) -> None:
        database = self.root / "blast_db"
        database.mkdir()
        for extension in ("nin", "nhr", "nsq"):
            (database / f"reference.{extension}").write_text(
                "index\n", encoding="utf-8"
            )
        rows = validator.validate_contract(self.args(database=str(database)))
        self.assertIn("prebuilt", rows[1]["details"])

    def test_blast_source_is_optional_but_not_ambiguous(self) -> None:
        rows = validator.validate_contract(self.args())
        self.assertEqual(rows[1]["status"], "skipped")
        self.assertEqual(rows[2]["status"], "skipped")
        with self.assertRaisesRegex(ValueError, "at most one"):
            validator.validate_contract(
                self.args(reference=str(self.reference), database=str(self.root))
            )

    def test_blast_thresholds_are_required_and_bounded(self) -> None:
        missing = self.args(reference=str(self.reference))
        missing.blast_min_identity = None
        with self.assertRaisesRegex(ValueError, "bacterial_blast_min_identity"):
            validator.validate_contract(missing)

        invalid = self.args(reference=str(self.reference))
        invalid.blast_min_coverage = 101
        with self.assertRaisesRegex(ValueError, "between 0 and 100"):
            validator.validate_contract(invalid)

    def test_virulence_database_requires_all_declared_files(self) -> None:
        (self.virulence / "virulence_test.fsa").unlink()
        with self.assertRaisesRegex(ValueError, "does not exist"):
            validator.validate_contract(self.args(reference=str(self.reference)))

    def test_prebuilt_database_rejects_ambiguous_indexes(self) -> None:
        database = self.root / "blast_db"
        database.mkdir()
        for prefix in ("first", "second"):
            for extension in ("nin", "nhr", "nsq"):
                (database / f"{prefix}.{extension}").write_text(
                    "index\n", encoding="utf-8"
                )
        with self.assertRaisesRegex(ValueError, "exactly one .nin"):
            validator.validate_contract(self.args(database=str(database)))


class GeneralTypingManifestTests(unittest.TestCase):
    def test_no_hit_outputs_are_completed_not_failed(self) -> None:
        with tempfile.TemporaryDirectory() as tempdir:
            root = Path(tempdir)
            virulence = root / "virulence.tsv"
            plasmid = root / "plasmid.tsv"
            mlst = root / "mlst.tsv"
            blast = root / "blast.txt"
            virulence.write_text("Database\tVirulence factor\n", encoding="utf-8")
            plasmid.write_text("Database\tPlasmid\n", encoding="utf-8")
            mlst.write_text("sample.fa\t-\t-\n", encoding="utf-8")
            blast.write_text("", encoding="utf-8")
            args = Namespace(
                sample_id="sample",
                genome_id="genome",
                source_route="standalone",
                virulencefinder_tsv=str(virulence),
                plasmidfinder_tsv=str(plasmid),
                mlst_tsv=str(mlst),
                blastn_output=str(blast),
                blastn_filtered_output=str(blast),
                blast_reference_type="user_supplied_reference_fasta",
                blast_reference_label="unit_test_reference",
            )
            rows = manifest.build_rows(args)
            self.assertEqual(len(rows), 4)
            self.assertTrue(all(row["execution_status"] == "completed" for row in rows))
            self.assertEqual(
                [row["result_state"] for row in rows],
                ["no_hits", "no_hits", "no_typing_match", "no_hits"],
            )

            output = root / "manifest.tsv"
            manifest.write_rows(output, rows)
            with output.open(encoding="utf-8") as handle:
                written = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual([row["tool"] for row in written], [
                "virulencefinder",
                "plasmidfinder",
                "mlst",
                "blastn",
            ])
            self.assertIn("only", written[-1]["interpretation"])
            self.assertEqual(written[-1]["filtered_output"], blast.name)

    def test_skipped_general_outputs_are_reported_as_skipped(self) -> None:
        with tempfile.TemporaryDirectory() as tempdir:
            root = Path(tempdir)
            virulence = root / "genome.virulencefinder.skipped.tsv"
            plasmid = root / "plasmid.tsv"
            mlst = root / "mlst.tsv"
            blast = root / "genome.blastn.skipped.txt"
            virulence.write_text("status\treason\nskipped\tno database\n", encoding="utf-8")
            plasmid.write_text("Database\tPlasmid\n", encoding="utf-8")
            mlst.write_text("sample.fa\t-\t-\n", encoding="utf-8")
            blast.write_text("BLASTN skipped\n", encoding="utf-8")
            args = Namespace(
                sample_id="sample",
                genome_id="genome",
                source_route="standalone",
                virulencefinder_tsv=str(virulence),
                plasmidfinder_tsv=str(plasmid),
                mlst_tsv=str(mlst),
                blastn_output=str(blast),
                blastn_filtered_output=str(blast),
                blast_reference_type="not_configured",
                blast_reference_label="",
            )
            rows = manifest.build_rows(args)
            self.assertEqual(rows[0]["execution_status"], "skipped")
            self.assertEqual(rows[0]["result_state"], "not_configured")
            self.assertEqual(rows[-1]["execution_status"], "skipped")
            self.assertEqual(rows[-1]["result_state"], "not_configured")


if __name__ == "__main__":
    unittest.main()
