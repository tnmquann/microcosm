import csv
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import gtdbtk_to_taxgenome as tool  # noqa: E402

COLUMNS = ["user_genome", "classification", "closest_genome_reference", "closest_genome_ani", "classification_method", "note", "warnings"]
BAC = [
    ["MEGAHIT-MetaBAT2-s1.1", "d__Bacteria;p__Pseudomonadota;c__G;o__E;f__E;g__Escherichia;s__Escherichia coli", "GCF_1", "98.9", "ani_screen", "", "N/A"],
    ["MEGAHIT-MetaBAT2-s1.2", "d__Bacteria;p__Bacillota_A;c__C;o__O;f__F;g__G;s__", "N/A", "N/A", "taxonomic classification defined by topology and ANI", "", "N/A"],
    ["MEGAHIT-MetaBAT2-s1.3", "Unclassified Bacteria", "N/A", "N/A", "N/A", "Insufficient number of amino acids in MSA", "N/A"],
    ["MEGAHIT-MetaBAT2-s1.4", "Unclassified", "N/A", "N/A", "N/A", "", "N/A"],
]
ARC = [["MEGAHIT-MetaBAT2-s1.0", "d__Archaea;p__Methanobacteriota;c__M;o__M;f__M;g__Methanobrevibacter;s__Methanobrevibacter smithii", "GCF_9", "99.1", "ani_screen", "", "N/A"]]


def write(path: Path, rows: list[list[str]]) -> Path:
    with path.open("w", newline="") as handle:
        csv.writer(handle, delimiter="\t").writerows([COLUMNS] + rows)
    return path


class GtdbtkToTaxgenomeTest(unittest.TestCase):
    def test_convert(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            out = tmp / "out.csv"
            self.assertEqual(tool.main(["--summary", str(write(tmp / "b.tsv", BAC)), str(write(tmp / "a.tsv", ARC)), "--output", str(out)]), 0)
            with out.open() as handle:
                rows = {row["query_name"]: row for row in csv.DictReader(handle)}
        self.assertEqual(list(rows), sorted(rows))
        self.assertEqual(rows["MEGAHIT-MetaBAT2-s1.1"]["rank"], "species")
        self.assertEqual(rows["MEGAHIT-MetaBAT2-s1.1"]["closest_genome_ani"], "98.9")
        # empty species rank dropped; GTDB suffix kept (routing strips it)
        self.assertEqual(rows["MEGAHIT-MetaBAT2-s1.2"]["lineage"], "d__Bacteria;p__Bacillota_A;c__C;o__O;f__F;g__G")
        self.assertEqual(rows["MEGAHIT-MetaBAT2-s1.2"]["rank"], "genus")
        self.assertEqual((rows["MEGAHIT-MetaBAT2-s1.3"]["lineage"], rows["MEGAHIT-MetaBAT2-s1.3"]["rank"]), ("d__Bacteria", "superkingdom"))
        self.assertEqual((rows["MEGAHIT-MetaBAT2-s1.4"]["status"], rows["MEGAHIT-MetaBAT2-s1.4"]["lineage"]), ("nomatch", ""))
        self.assertTrue(rows["MEGAHIT-MetaBAT2-s1.0"]["lineage"].startswith("d__Archaea;"))

    def test_routing_contract(self):
        # the output must pass the bacterial routing normalizer unchanged
        sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "bacterial_analysis"))
        import normalize_genome_classifications as normalizer

        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            out = tmp / "out.csv"
            tool.main(["--summary", str(write(tmp / "b.tsv", BAC)), "--output", str(out)])
            resolved, unresolved = normalizer.normalize([out])
        self.assertEqual([row["query_name"] for row in unresolved], ["MEGAHIT-MetaBAT2-s1.4"])
        species = {row["query_name"]: row["species"] for row in resolved}
        self.assertEqual(species["MEGAHIT-MetaBAT2-s1.1"], "s__Escherichia coli")
        self.assertEqual(species["MEGAHIT-MetaBAT2-s1.2"], "")

    def test_rejects_other_tables(self):
        with tempfile.TemporaryDirectory() as tmp:
            bad = Path(tmp) / "bad.tsv"
            bad.write_text("query_name\tlineage\nx\td__Bacteria\n")
            self.assertEqual(tool.main(["--summary", str(bad), "--output", str(Path(tmp) / "o.csv")]), 2)


if __name__ == "__main__":
    unittest.main()
