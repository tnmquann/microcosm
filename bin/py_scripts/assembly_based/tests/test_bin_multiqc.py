import csv
import json
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import bin_multiqc  # noqa: E402

TAX = """query_name,status,rank,fraction,lineage,query_md5,query_filename,f_weighted_at_rank,bp_match_at_rank,query_ani_at_rank
MEGAHIT-COMEBin-s1.1,match,species,0.7,d__Bacteria;p__P;c__C;o__O;f__F;g__Salmonella;s__Salmonella enterica,x,f,0.7,1,0.99
MEGAHIT-COMEBin-s1.2,match,genus,0.4,d__Bacteria;p__P;c__C;o__O;f__F;g__Klebsiella,x,f,0.4,1,0.95
MEGAHIT-MetaBAT2-s1.1,nomatch,,0,,x,f,0,0,0
"""


def write(path: Path, rows: list[list[str]]) -> Path:
    with path.open("w", newline="") as handle:
        csv.writer(handle, delimiter="\t").writerows(rows)
    return path


class BinMultiqcTest(unittest.TestCase):
    def test_outputs(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            (tmp / "tax.csv").write_text(TAX)
            checkm2 = write(tmp / "checkm2_summary.tsv", [["Name", "Completeness", "Contamination", "Completeness_Model_Used"], ["MEGAHIT-COMEBin-s1.1", "98.5", "1.2", "NN"]])
            busco = write(tmp / "busco_summary.tsv", [["Input_file", "Dataset", "Complete"], ["MEGAHIT-COMEBin-s1.1.fa", "bacteria_odb12", "90.1"]])
            quast = write(tmp / "quast_summary.tsv", [["Assembly", "# contigs (>= 0 bp)", "Total length (>= 0 bp)", "N50", "GC (%)"],
                                                       ["MEGAHIT-COMEBin-s1.1.fa.gz", "10", "5000000", "90000", "52.1"],
                                                       ["MEGAHIT-COMEBin-s1.2.fa.gz", "20", "4000000", "50000", "57.0"]])
            out = tmp / "out"
            self.assertEqual(bin_multiqc.main(["--taxgenome", str(tmp / "tax.csv"), "--qc", str(checkm2), str(busco), "--quast", str(quast), "--outdir", str(out)]), 0)

            plot = json.loads((out / "bin_taxonomy_mqc.json").read_text())
            labels = plot["pconfig"]["data_labels"]
            species = plot["data"][labels.index("Species")]
            self.assertEqual(species["MEGAHIT-COMEBin-s1"], {"s__Salmonella enterica": 1, "zzy_other": 0, "zzz_unclassified": 1})
            self.assertEqual(species["MEGAHIT-MetaBAT2-s1"]["zzz_unclassified"], 1)
            genus = plot["data"][labels.index("Genus")]
            self.assertEqual(genus["MEGAHIT-COMEBin-s1"]["g__Klebsiella"], 1)

            table = json.loads((out / "bin_summary_mqc.json").read_text())
            row = table["data"]["MEGAHIT-COMEBin-s1.1"]
            self.assertEqual((row["checkm2_completeness"], row["busco_complete"], row["n50"], row["taxon"]), (98.5, 90.1, 90000.0, "s__Salmonella enterica"))
            self.assertEqual(table["data"]["MEGAHIT-MetaBAT2-s1.1"]["taxon"], "unclassified")
            self.assertNotIn("checkm_completeness", table["headers"])

            report = list(csv.reader((out / "report_bins.tsv").open(), delimiter="\t"))
            self.assertEqual(report[0], ["Assembly", "MEGAHIT-COMEBin-s1.1.fa.gz", "MEGAHIT-COMEBin-s1.2.fa.gz"])
            self.assertIn(["N50", "90000", "50000"], report)

    def test_no_inputs_writes_nothing(self):
        with tempfile.TemporaryDirectory() as tmp:
            self.assertEqual(bin_multiqc.main(["--outdir", tmp]), 0)
            self.assertEqual(list(Path(tmp).iterdir()), [])


if __name__ == "__main__":
    unittest.main()
