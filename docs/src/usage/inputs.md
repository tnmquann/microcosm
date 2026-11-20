# Inputs

## Input Specifications
**microcosm** supports two main read input formats via the `--input_format` parameter: `directory` and `csv`. Genome FASTA files for the bacterial genome analysis are given separately with `--genome_input` ([see below](#genome-input)).

### Samplesheet CSV Input

You can specify a CSV samplesheet containing paths to your FASTQ files and optional metadata. Run merging occurs after short read QC (adapter clipping, host removal, etc.), and before taxonomic profiling and assembly.

A valid samplesheet can contain the following columns:

`sample_id,run_id,group,short_reads_1,short_reads_2,long_reads`

Example available at `/test_data/sample_e2e.csv`.

**Requirements:**
- Minimum 5 comma-separated columns
- File extension: `.csv`
- Must contain header: `sample_id,group,short_reads_1,short_reads_2,long_reads` (`run_id` is optional)
- Run IDs must be unique within a multi-run sample; samples with multiple runs are automatically concatenated
- Recommended: FastQ files should be compressed (`.fastq.gz`, `.fq.gz`)
- Only single-end or paired-end short reads per samplesheet
- For single-end reads, specify `--single_end` and provide reads in `short_reads_1`

**Example command:**
```bash
nextflow run main.nf \
    --input /path/to/your/samples.csv \
    --input_format csv
```
### Directory Input

This format currently supports only paired-end reads. Supported file patterns: `*_{1,2}*.{fastq,fastq.gz,fq,fq.gz}`.

**Example command:**
```bash
nextflow run main.nf \
    --input /path/to/your/fastq_directory \
    --input_format directory
```

### Genome Input

`--genome_input` takes a CSV of genome FASTA files (assemblies or MAGs made elsewhere) and runs only the [bacterial genome analysis](cookbook.md#bacterial-genome-analysis): read preprocessing, assembly and binning are skipped, and `--input` is not needed.

| Column | Required | Description |
| :--- | :--- | :--- |
| `sample_id` | Yes | Genome name, used in the output file names. |
| `fasta` | Yes | Path to the FASTA file (`.fa`, `.fasta`, `.fna`, optionally `.gz`). A relative path is resolved against the folder of the CSV file. |
| `lineage` | No | GTDB lineage (`d__Bacteria;p__...;s__...`). When valid, the genome is not reclassified; when blank, it is classified with sourmash (`--sourmash_database`, `--sourmash_taxonomy_csv`) or GTDB-Tk (`--bin_classifier gtdbtk --gtdbtk_db`). |
| `run_id`, `group`, `assembler` | No | Kept as metadata. |

```text
sample_id,fasta,lineage
kpneumoniae_HS11286,genomes/kpneumoniae_HS11286.fasta.gz,
styphi_CT18,genomes/styphi_CT18.fasta.gz,d__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae;g__Salmonella;s__Salmonella enterica
```

Only genomes resolved as `d__Bacteria` are analysed; the others are listed in `Bacterial_Analysis/routing/excluded_genomes.tsv` or `unresolved_genomes.tsv` with the reason. A ready-to-use example is in `test_data/bacterial_genomes/genome_input.csv`.
