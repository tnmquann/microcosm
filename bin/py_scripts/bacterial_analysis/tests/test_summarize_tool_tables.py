from __future__ import annotations

import csv
import importlib.util
import tempfile
import unittest
from pathlib import Path

SCRIPT = Path(__file__).resolve().parents[1] / "summarize_tool_tables.py"
SPEC = importlib.util.spec_from_file_location("summarize_tool_tables", SCRIPT)
tool = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(tool)


def read(path: Path) -> list[dict]:
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


class SummarizeToolTablesTest(unittest.TestCase):
    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp())

    def write(self, name: str, text: str) -> Path:
        path = self.tmp / name
        path.write_text(text)
        return path

    def manifest(self, entries) -> Path:
        return self.write("manifest.tsv", "".join(f"{g}\t{p}\n" for g, p in entries))

    def test_header_union_sorted_and_hash_stripped(self):
        b = self.write("b.tsv", "#FILE\tGENE\nb.fa\tblaTEM\n")
        a = self.write("a.tsv", "#FILE\tGENE\tEXTRA\na.fa\tsul1\tx\na.fa\ttetA\ty\n")
        empty = self.write("c.tsv", "#FILE\tGENE\n")
        out = self.tmp / "out.tsv"
        count = tool.summarize(self.manifest([("g_b", b), ("g_a", a), ("g_c", empty)]), out, False, [])
        rows = read(out)
        self.assertEqual(count, 3)
        self.assertEqual(list(rows[0]), ["genome_id", "FILE", "GENE", "EXTRA"])
        self.assertEqual([r["genome_id"] for r in rows], ["g_a", "g_a", "g_b"])
        self.assertEqual(rows[2]["EXTRA"], "")

    def test_no_header_mlst(self):
        m = self.write("m.tsv", "x.fa\tecoli\t131\tadk(53)\tfumC(40)\n")
        out = self.tmp / "out.tsv"
        tool.summarize(self.manifest([("g1", m)]), out, True, ["file", "scheme", "st", "alleles"])
        rows = read(out)
        self.assertEqual(rows, [{"genome_id": "g1", "file": "x.fa", "scheme": "ecoli", "st": "131",
                                 "alleles": "adk(53) fumC(40)"}])

    def test_empty_inputs_write_header_only(self):
        e = self.write("e.tsv", "")
        out = self.tmp / "out.tsv"
        self.assertEqual(tool.summarize(self.manifest([("g1", e)]), out, False, []), 0)
        self.assertEqual(out.read_text(), "genome_id\n")


if __name__ == "__main__":
    unittest.main()
