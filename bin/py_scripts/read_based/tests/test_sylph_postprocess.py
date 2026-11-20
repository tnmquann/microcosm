from __future__ import annotations

import csv
import importlib.util
import tempfile
import unittest
from pathlib import Path

SCRIPT = Path(__file__).resolve().parents[1] / "sylph_postprocess.py"
SPEC = importlib.util.spec_from_file_location("sylph_postprocess", SCRIPT)
tool = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(tool)

HEADER = "clade_name\trelative_abundance\tsequence_abundance\tANI (if strain-level)\tCoverage (if strain-level)\n"
LINEAGE = "d__Bacteria|p__Bacillota|c__Bacilli|o__Staphylococcales|f__Staphylococcaceae|g__Staphylococcus"


def sylphmpa(rows: list[tuple[str, str, str]]) -> str:
    body = "".join(f"{clade}\t{rel}\t{seq}\tNA\tNA\n" for clade, rel, seq in rows)
    return "#SampleID\treads_1.fastq.gz\tTaxonomies_used:['x']\n" + HEADER + body


class SylphPostprocessTest(unittest.TestCase):
    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp())
        a = self.tmp / "a.sylphmpa"
        b = self.tmp / "b.sylphmpa"
        a.write_text(sylphmpa([
            ("d__Bacteria", "100", "30.0"),
            (f"{LINEAGE}|s__Staphylococcus aureus", "100", "30.0"),
            (f"{LINEAGE}|s__Staphylococcus aureus|t__GCF_1.1", "100", "30.0"),
        ]))
        b.write_text(sylphmpa([
            ("d__Bacteria", "100", "12.0"),
            (f"{LINEAGE}|s__Staphylococcus epidermidis|t__GCF_2.1", "100", "12.0"),
        ]))
        self.manifest = self.tmp / "manifest.tsv"
        self.manifest.write_text(f"sB\t{b}\nsA\t{a}\n")

    def read(self, name: str, delimiter: str = ",") -> list[list[str]]:
        with (self.tmp / "out" / name).open() as handle:
            return list(csv.reader(handle, delimiter=delimiter))

    def test_sequence_abundance_matrix_and_phyloseq(self):
        tool.run(self.manifest, "sequence_abundance", self.tmp / "out", "run1", True)
        matrix = self.read("run1_sylphtax_sequence_abundance.tsv", "\t")
        self.assertEqual(matrix[0], ["clade_name", "sA", "sB"])
        self.assertIn(["d__Bacteria", "30", "12"], matrix)
        otu = self.read("phyloseq/otu_table.csv")
        self.assertEqual(otu, [["match_name", "sA", "sB"], ["GCF_1.1", "30", "0"], ["GCF_2.1", "0", "12"]])
        tax = self.read("phyloseq/tax_table.csv")
        self.assertEqual(tax[0], ["match_name"] + tool.RANKS)
        self.assertEqual(tax[1][0], "GCF_1.1")
        self.assertEqual(tax[1][-1], "s__Staphylococcus aureus")
        samples = self.read("phyloseq/sample_data.csv")
        self.assertEqual(samples[1], ["sA", "a.sylphmpa"])

    def test_taxonomic_abundance_without_phyloseq(self):
        tool.run(self.manifest, "taxonomic_abundance", self.tmp / "out", "run1", False)
        self.assertTrue((self.tmp / "out" / "run1_sylphtax_taxonomic_abundance.tsv").exists())
        self.assertFalse((self.tmp / "out" / "phyloseq").exists())


if __name__ == "__main__":
    unittest.main()
