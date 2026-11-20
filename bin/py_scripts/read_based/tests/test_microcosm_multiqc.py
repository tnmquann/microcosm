from __future__ import annotations

import importlib.util
import json
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


def row(sample, lineage, value, match="", level="min_coverage0.05"):
    return {"sample_name": sample, "lineage": lineage, "f_weighted_at_rank": str(value),
            "match_name": match, "min_coverage": level}


class MicrocosmMultiqcTest(unittest.TestCase):
    def test_rank_tabs_top_other_unclassified(self):
        rows = [row("s1", ECOLI, 0.5, "GCF_1"), row("s1", KPN, 0.1, "GCF_2"), row("s1", SAUR, 0.05, "GCF_3"),
                row("s1", "unclassified", 0.35, "unclassified", ""), row("s2", ECOLI, 0.2, "GCF_1"),
                row("s2", "unclassified", 0.8, "unclassified", "")]
        datasets, categories, labels = tool.build(rows, top=1)
        self.assertEqual(labels, ["Strain", "Species", "Genus", "Family", "Order", "Class", "Phylum", "Domain"])
        species = datasets[labels.index("Species")]
        self.assertEqual(list(species["s1"]), ["s__Escherichia coli", tool.OTHER, tool.UNCLASSIFIED])
        self.assertAlmostEqual(species["s1"][tool.OTHER], 15.0)
        self.assertAlmostEqual(species["s2"][tool.UNCLASSIFIED], 80.0)
        self.assertEqual(list(datasets[0]["s1"])[0], "GCF_1")
        self.assertEqual(categories[0][tool.OTHER]["name"], "Other")

    def test_levels_and_json(self):
        # post_processing.py layout: a genome lists every level it passes; unclassified rows carry no level
        rows = [row("s1", ECOLI, 0.4, "GCF_1", "min_coverage0.05 | min_coverage0.01"), row("s1", KPN, 0.1, "GCF_2", "min_coverage0.01"),
                row("s1", "unclassified", 0.6, "unclassified", ""), row("s1", "unclassified", 0.5, "unclassified", "")]
        datasets, categories, labels = tool.build(rows, top=5, min_coverage=0.05)
        self.assertEqual(datasets[0], {"s1": {"GCF_1": 40.0, tool.OTHER: 0.0, tool.UNCLASSIFIED: 60.0}})
        low = tool.build(rows, top=5, min_coverage=0.01)[0]
        self.assertEqual(low[1]["s1"][tool.UNCLASSIFIED], 50.0)
        # --remove_unclassified: no Unclassified category, taxa rescaled to 100 %
        rescaled, cats, _ = tool.build(rows, top=5, min_coverage=0.01, remove_unclassified=True)
        self.assertEqual(rescaled[1]["s1"], {"s__Escherichia coli": 80.0, "s__Klebsiella pneumoniae": 20.0, tool.OTHER: 0.0})
        self.assertNotIn(tool.UNCLASSIFIED, cats[1])
        self.assertEqual(tool.build(rows, top=5, min_coverage=0.5)[0], [])
        out = Path(tempfile.mkdtemp()) / "microcosm_mqc.json"
        tool.write_json(out, datasets, categories, labels, 5, 0.05, False)
        content = json.loads(out.read_text())
        self.assertEqual(content["plot_type"], "bargraph")
        self.assertEqual(content["pconfig"]["data_labels"], labels)
        self.assertEqual(len(content["data"]), len(content["categories"]))
        self.assertIn("min_coverage 0.05", content["description"])


if __name__ == "__main__":
    unittest.main()
