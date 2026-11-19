import gzip
import json
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import contig_multiqc as tool  # noqa: E402


class ContigMultiqcTest(unittest.TestCase):
    def test_outputs(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            (tmp / "MEGAHIT-s1.seqkit.txt").write_text("#id\tlength\tGC\nk141_1\t1000\t50\nk141_2\t3000\t40\nk141_3\t6000\t60\n")
            (tmp / "MEGAHIT-s1.coverage.txt").write_text(
                "#rname\tstartpos\tendpos\tnumreads\tcovbases\tcoverage\tmeandepth\tmeanbaseq\tmeanmapq\n"
                "k141_1\t1\t1000\t10\t1000\t100\t2\t30\t60\nk141_2\t1\t3000\t50\t1500\t50\t10\t30\t60\n")
            (tmp / "MEGAHIT-s1.pyrodigal.faa").write_text(">g1\nMK\n>g2\nMK\n>g3\nMK\n")
            (tmp / "MEGAHIT-s1.contigs_plasmid_summary.tsv").write_text("seq_name\tlength\ttopology\tamr_genes\nk141_2\t3000\tDTR\tNA\n")
            (tmp / "MEGAHIT-s1.contigs_virus_summary.tsv").write_text(
                "seq_name\tlength\ttopology\nk141_3|provirus_1_200\t200\tProvirus\nk141_1\t1000\tNo terminal repeats\n")
            with gzip.open(tmp / "MEGAHIT-s1.tiara.txt.gz", "wt") as handle:
                handle.write("sequence_id\tclass_fst_stage\tclass_snd_stage\nk141_2\tbacteria\tn/a\nk141_3 flag=1\teukarya\tunknown\n")
            out = tmp / "out"
            self.assertEqual(tool.main([
                "--seqkit", str(tmp / "MEGAHIT-s1.seqkit.txt"), "--coverage", str(tmp / "MEGAHIT-s1.coverage.txt"),
                "--faa", str(tmp / "MEGAHIT-s1.pyrodigal.faa"), "--plasmid", str(tmp / "MEGAHIT-s1.contigs_plasmid_summary.tsv"),
                "--virus", str(tmp / "MEGAHIT-s1.contigs_virus_summary.tsv"), "--tiara", str(tmp / "MEGAHIT-s1.tiara.txt.gz"),
                "--outdir", str(out)]), 0)
            row = json.loads((out / "contig_summary_mqc.json").read_text())["data"]["MEGAHIT-s1"]
            self.assertEqual((row["contigs"], row["length"], row["mapped_reads"], row["genes"]), (3, 10000.0, 60.0, 3))
            self.assertAlmostEqual(row["mean_depth"], (2 * 1000 + 10 * 3000) / 4000)
            self.assertAlmostEqual(row["breadth"], 2500 / 4000 * 100)
            self.assertEqual((row["plasmids"], row["plasmid_amr"], row["viruses"], row["proviruses"], row["virus_bp"]), (1, 0, 2, 1, 1200.0))
            tiara = json.loads((out / "contig_tiara_mqc.json").read_text())
            self.assertEqual(tiara["pconfig"]["data_labels"], ["Contigs", "Length (bp)"])
            self.assertEqual(tiara["data"][0]["MEGAHIT-s1"]["eukarya"], 1)
            self.assertEqual(tiara["data"][1]["MEGAHIT-s1"]["bacteria"], 3000.0)

    def test_skani_best_hits(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            (tmp / "MEGAHIT-s1.seqkit.txt").write_text("#id\tlength\nk1\t5000\nk2\t3000\nk3\t2000\n")
            header = "Ref_file\tQuery_file\tANI\tAlign_fraction_ref\tAlign_fraction_query\tRef_name\tQuery_name\n"
            (tmp / "MEGAHIT-s1.skani.tsv").write_text(header +
                "kpn.fa\tq\t99.9\t0.8\t99\tNC_1 Klebsiella pneumoniae HS11286 chromosome, complete genome\tk1 flag=0 len=5000\n"
                "kae.fa\tq\t92.0\t0.7\t80\tNC_2 Klebsiella aerogenes KCTC 2190, complete sequence\tk1 flag=0 len=5000\n"
                "sty.fa\tq\t98.0\t0.5\t90\tNC_3 Salmonella enterica Typhi CT18\tk2 flag=0 len=3000\n")
            tool.main(["--seqkit", str(tmp / "MEGAHIT-s1.seqkit.txt"), "--skani", str(tmp / "MEGAHIT-s1.skani.tsv"), "--top", "1", "--outdir", str(tmp)])
            plot = json.loads((tmp / "contig_skani_mqc.json").read_text())
            self.assertEqual(plot["categories"], ["Klebsiella pneumoniae HS11286", tool.SKANI_OTHER, tool.SKANI_NONE])
            self.assertEqual(plot["data"][0]["MEGAHIT-s1"], {"Klebsiella pneumoniae HS11286": 5000.0, tool.SKANI_OTHER: 3000.0, tool.SKANI_NONE: 2000.0})
            self.assertEqual(plot["data"][1]["MEGAHIT-s1"][tool.SKANI_NONE], 1)
            row = json.loads((tmp / "contig_summary_mqc.json").read_text())["data"]["MEGAHIT-s1"]
            self.assertEqual((row["skani_contigs"], row["skani_bp"]), (2, 80.0))

    def test_skipped_tools_leave_no_columns(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            (tmp / "SPAdes-s2.seqkit.txt").write_text("#id\tlength\nNODE_1\t500\n")
            tool.main(["--seqkit", str(tmp / "SPAdes-s2.seqkit.txt"), "--outdir", str(tmp)])
            table = json.loads((tmp / "contig_summary_mqc.json").read_text())
            self.assertEqual(list(table["headers"]), ["contigs", "length"])
            self.assertFalse((tmp / "contig_tiara_mqc.json").exists())


if __name__ == "__main__":
    unittest.main()
