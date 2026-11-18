from __future__ import annotations

import csv
import importlib.util
import tempfile
import unittest
from pathlib import Path

SCRIPT = Path(__file__).resolve().parents[1] / "microcosm_multiqc.py"
SPEC = importlib.util.spec_from_file_location("microcosm_multiqc", SCRIPT)
tool = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(tool)

ECOLI = "d__Bacteria;p__Pseudomonadota;c__G;o__E;f__E;g__Escherichia;s__Escherichia coli"
KPN = "d__Bacteria;p__Pseudomonadota;c__G;o__E;f__E;g__Klebsiella;s__Klebsiella pneumoniae"
SAUR = "d__Bacteria;p__Bacillota;c__B;o__S;f__S;g__Staphylococcus;s__Staphylococcus aureus"


def row(sample, lineage, value, level="min_coverage0.05"):
    return {"sample_name": sample, "lineage": lineage, "f_weighted_at_rank": str(value), "min_coverage": level}


class MicrocosmMultiqcTest(unittest.TestCase):
    def test_top_other_unclassified(self):
        rows = [row("s1", ECOLI, 0.5), row("s1", KPN, 0.1), row("s1", SAUR, 0.05), row("s1", "unclassified", 0.35),
                row("s2", ECOLI, 0.2), row("s2", "unclassified", 0.8)]
        tables = tool.build_tables(rows, top=1)
        species = tables["species"]
        self.assertEqual(list(species["s1"]), ["s__Escherichia coli", "Other", "unclassified"])
        self.assertAlmostEqual(species["s1"]["Other"], 15.0)
        self.assertAlmostEqual(species["s2"]["unclassified"], 80.0)
        self.assertAlmostEqual(tables["phylum"]["s1"]["p__Pseudomonadota"], 60.0)

    def test_levels_become_samples_and_file_written(self):
        rows = [row("s1", ECOLI, 0.4, "min_coverage0.05"), row("s1", ECOLI, 0.5, "min_coverage0.01"),
                row("s1", "unclassified", 0.5, "")]
        tables = tool.build_tables(rows, top=5)
        self.assertEqual(sorted(tables["genus"]), ["s1 (min_coverage0.01)", "s1 (min_coverage0.05)"])
        self.assertAlmostEqual(tables["genus"]["s1 (min_coverage0.05)"]["unclassified"], 50.0)
        out = Path(tempfile.mkdtemp()) / "x_mqc.tsv"
        tool.write_table(out, "genus", tables["genus"])
        lines = out.read_text().splitlines()
        self.assertTrue(lines[0].startswith("# id: 'microcosm_genus'"))
        data = list(csv.reader([l for l in lines if not l.startswith("#")], delimiter="\t"))
        self.assertEqual(data[0], ["Sample", "g__Escherichia", "Other", "unclassified"])


if __name__ == "__main__":
    unittest.main()
