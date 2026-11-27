"""f_taxonomic_at_rank / match_genome_bp and --datatype in post_processing.py (raw_microcosm input).

Needs pandas + numpy (the POSTPROCESS_READBASED container); skipped where they are missing.
"""
from __future__ import annotations

import csv
import importlib.util
import sys
import tempfile
import unittest
from pathlib import Path

try:
    import pandas  # noqa: F401
except ImportError:  # pragma: no cover - host Python without pandas
    pandas = None

SCRIPT = Path(__file__).resolve().parents[1] / "post_processing.py"
HEADER = [
    "WGS_ID", "lineage", "organism_name", "query_md5", "query_filename", "ksize", "scaled",
    "f_unique_to_query", "unique_intersect_bp", "intersect_bp", "f_match", "query_bp",
    "total_weighted_hashes", "median_abund", "average_abund", "relative_abundance", "file_name",
    "p_vals", "actual_confidence_with_coverage", "alt_confidence_mut_rate_with_coverage",
]
LINEAGE = "d__Bacteria;p__P;c__C;o__O;f__F;g__G;s__G {}"


def genome_row(sample: str, name: str, genome_bp: int, unique_bp: int, *, f_match: float | None = None) -> dict:
    # scaled 1000, median_abund 1: sequence share = (unique_bp / 1000) / total_weighted_hashes
    return {
        "WGS_ID": sample,
        "lineage": LINEAGE.format(name),
        "organism_name": name,
        "query_md5": f"md5{sample}",
        "query_filename": f"{sample}.fq.gz",
        "ksize": 31,
        "scaled": 1000,
        "f_unique_to_query": 0.1,
        "unique_intersect_bp": unique_bp,
        "intersect_bp": unique_bp + 7000,  # overlap with an earlier match must not change the genome size
        "f_match": unique_bp / genome_bp if f_match is None else f_match,
        "query_bp": 10_000_000,
        "total_weighted_hashes": 10_000,
        "median_abund": 1,
        "average_abund": 1,
        "relative_abundance": 0.5,
        "file_name": "min_coverage0.05",
        "p_vals": 0,
        "actual_confidence_with_coverage": 0.95,
        "alt_confidence_mut_rate_with_coverage": 0.95,
    }


@unittest.skipIf(pandas is None, "pandas is not installed")
class TaxonomicAbundanceTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        spec = importlib.util.spec_from_file_location("post_processing", SCRIPT)
        cls.tool = importlib.util.module_from_spec(spec)
        sys.modules[spec.name] = cls.tool  # dataclasses look the module up while the script loads
        spec.loader.exec_module(cls.tool)

    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp())

    def write_input(self, rows: list[dict]) -> Path:
        path = self.tmp / "merged_sourmash_yacht.csv"
        with path.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=HEADER)
            writer.writeheader()
            writer.writerows(rows)
        return path

    def run_tool(self, path: Path, output_format: str, *extra: str) -> Path:
        out = self.tmp / output_format
        out.mkdir(exist_ok=True)
        code = self.tool.main([
            "--input", str(path), "--input_format", "raw_microcosm", "--min_coverage", "0.05",
            "--output", str(out), "--output_format", output_format, *extra,
        ])
        self.assertEqual(code, 0)
        return out

    def read_csv(self, path: Path) -> list[dict]:
        with path.open(newline="") as handle:
            return list(csv.DictReader(handle))

    def two_genomes(self) -> Path:
        # equal sequence share (2,000 hashes each = 20%), genomes of 4.2 Mb and 2.1 Mb
        return self.write_input([
            genome_row("s1", "big", 4_200_000, 2_000_000),
            genome_row("s1", "small", 2_100_000, 2_000_000),
        ])

    def test_genome_size_normalised_share(self):
        out = self.run_tool(self.two_genomes(), "default")
        rows = {row["match_name"]: row for row in self.read_csv(next(out.glob("*.processed_microcosm.csv")))}
        self.assertAlmostEqual(float(rows["big"]["f_taxonomic_at_rank"]), 1 / 3)
        self.assertAlmostEqual(float(rows["small"]["f_taxonomic_at_rank"]), 2 / 3)
        self.assertAlmostEqual(float(rows["big"]["match_genome_bp"]), 4_200_000)
        self.assertAlmostEqual(float(rows["small"]["match_genome_bp"]), 2_100_000)
        # sequence share unchanged
        self.assertAlmostEqual(float(rows["big"]["f_weighted_at_rank"]), 0.2)
        self.assertAlmostEqual(float(rows["unclassified"]["f_weighted_at_rank"]), 0.6)

    def test_unclassified_is_empty_and_classified_sums_to_one(self):
        out = self.run_tool(self.two_genomes(), "default")
        rows = self.read_csv(next(out.glob("*.processed_microcosm.csv")))
        unclassified = [row for row in rows if row["lineage"] == "unclassified"]
        self.assertEqual(len(unclassified), 1)
        self.assertEqual(unclassified[0]["f_taxonomic_at_rank"], "")
        self.assertEqual(unclassified[0]["match_genome_bp"], "")
        total = sum(float(row["f_taxonomic_at_rank"]) for row in rows if row["lineage"] != "unclassified")
        self.assertAlmostEqual(total, 1.0)

    def test_phyloseq_otu_table_follows_datatype(self):
        path = self.two_genomes()
        seq = self.run_tool(path, "phyloseq")
        seq_otu = {row["match_name"]: row["s1"] for row in self.read_csv(next(seq.glob("*.phyloseq")) / "otu_table.csv")}
        self.assertEqual(seq_otu, {"big": "2000", "small": "2000", "unclassified": "6000"})

        tax_dir = self.tmp / "phyloseq_tax"
        tax_dir.mkdir()
        code = self.tool.main([
            "--input", str(path), "--input_format", "raw_microcosm", "--min_coverage", "0.05",
            "--output", str(tax_dir), "--output_format", "phyloseq", "--datatype", "taxonomic_abundance",
        ])
        self.assertEqual(code, 0)
        phyloseq = next(tax_dir.glob("*.phyloseq"))
        tax_otu = {row["match_name"]: float(row["s1"]) for row in self.read_csv(phyloseq / "otu_table.csv")}
        self.assertEqual(set(tax_otu), {"big", "small"})
        self.assertAlmostEqual(tax_otu["big"], 1 / 3)
        self.assertAlmostEqual(tax_otu["small"], 2 / 3)
        taxa = {row["match_name"] for row in self.read_csv(phyloseq / "tax_table.csv")}
        self.assertEqual(taxa, {"big", "small"})
        self.assertEqual([row["sample_name"] for row in self.read_csv(phyloseq / "sample_data.csv")], ["s1"])

    def test_zero_f_match_is_an_error(self):
        path = self.write_input([
            genome_row("s1", "big", 4_200_000, 2_000_000),
            genome_row("s1", "broken", 2_100_000, 2_000_000, f_match=0.0),
        ])
        out = self.tmp / "err"
        out.mkdir()
        with self.assertLogs(level="ERROR") as logs:
            code = self.tool.main([
                "--input", str(path), "--input_format", "raw_microcosm", "--min_coverage", "0.05",
                "--output", str(out), "--output_format", "default",
            ])
        self.assertEqual(code, 1)
        self.assertIn("s1 / broken", "\n".join(logs.output))
        self.assertIn("f_match must be > 0", "\n".join(logs.output))


if __name__ == "__main__":
    unittest.main()
