import gzip
import json
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import tool_multiqc as tool  # noqa: E402


class ToolMultiqcTest(unittest.TestCase):
    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.tmp = Path(self._tmp.name)
        self.out = self.tmp / "out"

    def tearDown(self):
        self._tmp.cleanup()

    def write(self, name: str, text: str) -> Path:
        path = self.tmp / name
        path.parent.mkdir(parents=True, exist_ok=True)
        if name.endswith(".gz"):
            with gzip.open(path, "wt") as handle:
                handle.write(text)
        else:
            path.write_text(text)
        return path

    def run_tool(self, name: str, *inputs: Path) -> dict:
        self.assertEqual(tool.main(["--tool", name, *map(str, inputs), "--outdir", str(self.out), "--top", "1"]), 0)
        return {path.name: json.loads(path.read_text()) for path in self.out.glob("*_mqc.json")}

    def seqkit(self) -> Path:
        return self.write("MEGAHIT-s1.seqkit.txt", "#id\tlength\tGC\nk1\t1000\t50\nk2\t3000\t40\nk3\t6000\t60\n")

    def test_contig_summary(self):
        cov = self.write("MEGAHIT-s1.coverage.txt", "#rname\tstartpos\tendpos\tnumreads\tcovbases\tcoverage\tmeandepth\n"
                         "k1\t1\t1000\t10\t1000\t100\t2\nk2\t1\t3000\t50\t1500\t50\t10\n")
        row = self.run_tool("contig_summary", self.seqkit(), cov)["contig_summary_mqc.json"]["data"]["MEGAHIT-s1"]
        self.assertEqual((row["contigs"], row["length"], row["mapped_reads"]), (3, 10000.0, 60.0))
        self.assertAlmostEqual(row["mean_depth"], (2 * 1000 + 10 * 3000) / 4000)
        self.assertAlmostEqual(row["breadth"], 2500 / 4000 * 100)

    def test_pyrodigal_and_genomad(self):
        faa = self.write("MEGAHIT-s1.pyrodigal.faa.gz", ">g1 # 1 # 9 # 1 # ID=1_1;partial=00\nMKV*\n>g2 # 1 # 9 # 1 # ID=1_2;partial=10\nMK\n")
        row = self.run_tool("pyrodigal", faa)["pyrodigal_mqc.json"]["data"]["MEGAHIT-s1"]
        self.assertEqual((row["genes"], row["complete"], row["partial"], row["mean_length"]), (2, 1, 1, 2.5))
        plasmid = self.write("g/MEGAHIT-s1.contigs_plasmid_summary.tsv", "seq_name\tlength\tconjugation_genes\tamr_genes\nk2\t3000\tNA\tblaTEM;sul1\n")
        virus = self.write("g/MEGAHIT-s1.contigs_virus_summary.tsv", "seq_name\tlength\ttopology\nk3|provirus_1_200\t200\tProvirus\nk1\t1000\tNo terminal repeats\n")
        row = self.run_tool("genomad", plasmid.parent)["genomad_mqc.json"]["data"]["MEGAHIT-s1"]
        self.assertEqual((row["plasmids"], row["plasmid_amr"], row["plasmid_conjugation"], row["viruses"], row["proviruses"], row["virus_bp"]),
                         (1, 2, 0, 2, 1, 1200.0))

    def test_tiara_contigs_weighs_length(self):
        tiara = self.write("MEGAHIT-s1.tiara.txt.gz", "sequence_id\tclass_fst_stage\tclass_snd_stage\nk2\tbacteria\tn/a\nk3 flag=1\teukarya\tunknown\n")
        plot = self.run_tool("tiara_contigs", tiara, self.seqkit())["contig_tiara_mqc.json"]
        self.assertEqual(plot["pconfig"]["data_labels"], ["Contigs", "Length (bp)"])
        self.assertEqual((plot["data"][0]["MEGAHIT-s1"]["eukarya"], plot["data"][1]["MEGAHIT-s1"]["bacteria"]), (1, 3000.0))

    def test_skani_best_hits(self):
        self.write("MEGAHIT-s1.seqkit.txt", "#id\tlength\nk1\t5000\nk2\t3000\nk3\t2000\n")
        header = "Ref_file\tQuery_file\tANI\tAlign_fraction_ref\tAlign_fraction_query\tRef_name\tQuery_name\n"
        skani = self.write("MEGAHIT-s1.skani.tsv", header +
                           "kpn.fa\tq\t99.9\t0.8\t99\tNC_1 Klebsiella pneumoniae HS11286 chromosome, complete genome\tk1 flag=0\n"
                           "kae.fa\tq\t92.0\t0.7\t80\tNC_2 Klebsiella aerogenes KCTC 2190, complete sequence\tk1 flag=0\n"
                           "sty.fa\tq\t98.0\t0.5\t90\tNC_3 Salmonella enterica Typhi CT18\tk2 flag=0\n")
        plot = self.run_tool("skani", skani, self.tmp / "MEGAHIT-s1.seqkit.txt")["contig_skani_mqc.json"]
        self.assertEqual(plot["categories"], ["Klebsiella pneumoniae HS11286", tool.OTHER, tool.SKANI_NONE])
        self.assertEqual(plot["data"][0]["MEGAHIT-s1"], {"Klebsiella pneumoniae HS11286": 5000.0, tool.OTHER: 3000.0, tool.SKANI_NONE: 2000.0})
        self.assertEqual(plot["data"][1]["MEGAHIT-s1"][tool.SKANI_NONE], 1)

    def test_bin_tools(self):
        tiara = self.write("MEGAHIT-MetaBAT2-rawbins-s1.binclassification.tsv",
                           "BinID\tprokarya\teukarya\tclassification\nb1\t0.9\t0.1\tbacteria\nb2\t0.1\t0.9\teukarya\nb3\t0.5\t0.5\tweird\n")
        bins = self.run_tool("tiara_bins", tiara)["bin_tiara_mqc.json"]["data"][0]["MEGAHIT-MetaBAT2-rawbins-s1"]
        self.assertEqual((bins["bacteria"], bins["eukarya"], bins["unknown"]), (1, 1, 1))
        # GUNC writes the same file name for every bin set: both copies are read
        header = "genome\tn_genes_called\tn_genes_mapped\tn_contigs\ttaxonomic_level\tclade_separation_score\tcontamination_portion\treference_representation_score\tpass.GUNC\n"
        self.write("in1/GUNC.progenomes_2.1.maxCSS_level.tsv", header + "bin_a\t4000\t3900\t200\tkingdom\t0.0\t0.0\t0.98\tTrue\n")
        self.write("in2/GUNC.progenomes_2.1.maxCSS_level.tsv", header + "bin_b\t5000\t4900\t900\tphylum\t0.99\t0.02\t0.98\tFalse\n")
        gunc = self.run_tool("gunc", self.tmp / "in1", self.tmp / "in2")["gunc_mqc.json"]["data"]
        self.assertEqual((gunc["bin_a"]["pass"], gunc["bin_b"]["clade_separation_score"]), ("True", 0.99))
        stats = self.write("x.assembly_stats.csv", "bin_name,assembly_length_bp,contig_count,contig_N50_bp,contig_N90_bp,GC_perc,gaps_sum_bp\nb1,1000,2,600,400,50.5,0\n")
        row = self.run_tool("assembly_stats", stats)["assembly_stats_mqc.json"]["data"]["b1"]
        self.assertEqual((row["contig_count"], row["GC_perc"]), (2.0, 50.5))

    def test_rgi_bwt_and_rgi_amr(self):
        self.write("s1/s1.overall_mapping_stats.txt", "Total reads:       1000\nMapped reads:      50\t(5%)\n")
        self.write("s1/s1.gene_mapping_data.txt", "ARO Term\tAll Mapped Reads\tAMR Gene Family\tDrug Class\n"
                   "TEM-1\t30\tTEM beta-lactamase\tpenam; cephalosporin\nsul1\t10\tsul\tsulfonamide\n")
        (self.tmp / "staged").mkdir()
        (self.tmp / "staged" / "s1").symlink_to(self.tmp / "s1")  # Nextflow stages output folders as symlinks
        written = self.run_tool("rgi_bwt", self.tmp / "staged")
        row = written["rgi_bwt_mqc.json"]["data"]["s1"]
        self.assertEqual((row["genes"], row["arg_reads"], row["drug_classes"], row["mapped_pct"]), (2, 40.0, 3, 5.0))
        self.assertEqual(written["rgi_bwt_families_mqc.json"]["data"][0]["s1"], {"TEM beta-lactamase": 30.0, tool.OTHER: 10.0})
        self.write("g1.rgi.txt", "ORF_ID\tCut_Off\tDrug Class\tAMR Gene Family\no1\tPerfect\tpenam\tTEM\no2\tStrict\tpenam; carbapenem\tKPC\n")
        row = self.run_tool("rgi_amr", self.tmp / "g1.rgi.txt")["rgi_amr_mqc.json"]["data"]["g1"]
        self.assertEqual((row["hits"], row["perfect"], row["strict"], row["drug_classes"]), (2, 1, 1, 2))

    def test_no_rows_writes_nothing(self):
        self.run_tool("gunc")
        self.assertEqual(list(self.out.glob("*")), [])


if __name__ == "__main__":
    unittest.main()
