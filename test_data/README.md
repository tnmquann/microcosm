# microcosm test data

Small inputs for local tests (each set under 20 MB). Large data is not stored here: simulated reads are rebuilt by a
script, and real databases are passed by path. Paths inside the `genome_input*.csv` files are relative to the CSV,
so the folder can be used from any checkout.

## Bacterial genomes (`bacterial_genomes/`)

Eight NCBI RefSeq genomes (gzipped FASTA, 11 MB):

| File | Organism | Assembly | Used for |
|---|---|---|---|
| `kpneumoniae_HS11286.fasta.gz` | *Klebsiella pneumoniae* HS11286 | GCF_000240185.1 | Kleborate, MLST ST11, sim reads |
| `koxytoca_KCTC1686.fasta.gz` | *Klebsiella oxytoca* KCTC 1686 | GCF_000247855.1 | Kleborate (KoSC) |
| `kaerogenes_KCTC2190.fasta.gz` | *Klebsiella aerogenes* KCTC 2190 | GCF_000215745.1 | supplied-lineage case |
| `ecoli_K12.fasta.gz` | *Escherichia coli* K-12 MG1655 | GCF_000005845.2 | MLST ST10, VirulenceFinder, DeepARG |
| `styphi_CT18.fasta.gz` | *Salmonella* Typhi CT18 | GCF_000195995.1 | SISTR -> GenoTyphi, sim reads |
| `sparatyphiA_ATCC9150.fasta.gz` | *Salmonella* Paratyphi A ATCC 9150 | GCF_000011885.1 | SISTR -> Paratype |
| `styphimurium_LT2.fasta.gz` | *Salmonella* Typhimurium LT2 | GCF_000006945.2 | SISTR, no serovar-specific tool |
| `msmithii_PS.fasta.gz` | *Methanobrevibacter smithii* (archaeon) | GCF_000016525.1 | excluded (non-bacterial) genome |

- `genome_input.csv`: Kpn, *E. coli* and the three *Salmonella* (empty lineage -> classified by the pipeline).
- `genome_input_all.csv`: all eight genomes.
- `decision_table_steps.tsv`: `--bacterial_decision_table` with y/n step columns only.
- `decision_table_params.tsv`: `--bacterial_decision_table` with per-lineage parameter columns (`mlst_args`,
  `amrfinderplus_args`, `amrfinderplus_organism`, `sistr_args`); the *Salmonella* values are test values.

```bash
nextflow run main.nf -profile singularity --genome_input test_data/bacterial_genomes/genome_input.csv \
  --sourmash_database <gtdb-rs226-k51 rocksdb> --sourmash_ksize 51 --sourmash_taxonomy_csv <gtdb-rs226 lineages> \
  --bacterial_decision_table test_data/bacterial_genomes/decision_table_params.tsv --outdir results_genomes
# GTDB-Tk instead of sourmash: --bin_classifier gtdbtk --gtdbtk_db <GTDB-Tk R226 data folder>
```

## skani database (`skani_db/`)

skani 0.2.2 sketches of the eight genomes above. Built with
`skani sketch <genomes>.fasta -o skani_db` (FASTA uncompressed, same file names without `.gz`); rebuild it with
the skani version of `modules/nf-core/skani/search` if that module is updated. Use: `--skani_db test_data/skani_db`
with `--skip_contig_analysis false --skip_skani_search false`.

## GTDB-Tk failure case (`gtdbtk_no_markers/`)

Two random 50 kb sequences (Python `random.seed(1)`) with no marker genes. `--genome_input
test_data/gtdbtk_no_markers/genome_input.csv --bin_classifier gtdbtk --gtdbtk_db <R226>` must finish and list both
genomes in `Bacterial_Analysis/routing/unresolved_genomes.tsv` with reason
`no_resolved_lineage_from_BIN_CLASSIFICATION (gtdbtk: No bacterial or archaeal marker)`.

## Several sourmash databases (`sourmash_multidb/`)

`mock_part1.zip` / `mock_part2.zip`: the five `GCA_/GCF_*_genomic.fna.gz` genomes at the top of this folder split
over two k=31 zip databases; `tax_mock_part{1,2}.csv`: their lineages. Read-based test with two databases and two
taxonomies (YACHT is skipped when several databases are given):

```bash
nextflow run main.nf -profile test,singularity \
  --sourmash_database test_data/sourmash_multidb/mock_part1.zip,test_data/sourmash_multidb/mock_part2.zip \
  --sourmash_taxonomy_csv test_data/sourmash_multidb/tax_mock_part1.csv,test_data/sourmash_multidb/tax_mock_part2.csv
```

## Simulated reads (`sim_reads/simulate.sh`)

`bash test_data/sim_reads/simulate.sh <outdir>` writes `simmix_{1,2}.fastq.gz` (K. pneumoniae HS11286 + S. Typhi CT18,
2x150 bp, ~20x each, fixed seeds; needs `wgsim`). The reads (~104 MB) are not kept here. They drive the
assembly -> bins -> classification -> bacterial downstream tests (14 bins with MEGAHIT + MetaBAT2/MaxBin2/SemiBin2/COMEBin).
Make a samplesheet pointing at the two files and pass it with `--input`.

## Older fixtures

- `minigut/`: two small paired-end read sets (read-based `-profile test`).
- `mockupdb_r226_k31*`, `signatures/`, `test_sketches.zip`, `test_batch.csv`: mock sourmash / YACHT databases for `conf/test.config`.
- `readbased_postprocess_*`: inputs of the read-based post-processing tests.
- `hostile/`: small host index for `conf/test.config`.
- `binette/`: bin folders for the Binette refinement tests.
- `genome_input/`: tiny genome-input validation cases (valid, mixed, missing lineage, archaea).
- `general_typing/`, `amr_databases/`: trimmed Paratype / VirulenceFinder / CARD fixtures.
