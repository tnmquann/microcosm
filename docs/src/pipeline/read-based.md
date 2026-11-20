# Read-based analysis

Read-based analysis answers *who is in the sample* and *which resistance genes are in the reads* without assembling.
Turn it on with `--enable_readbase`.

| What | Tools | Switch |
| :--- | :--- | :--- |
| Taxonomic profile | sourmash gather + YACHT (`microcosm`, default), sylph (+ sylph-tax) | `--readbased_taxotools microcosm\|sylph\|microcosm,sylph` |
| AMR genes in reads | RGI bwt against CARD | `--enable_rgi_bwt` |
| Ready-to-use tables | standard table, phyloseq, taxburst, RGI bwt tables | `--readbased_postprocess all\|phyloseq\|taxburst\|rgi_bwt` |

:::{tip}
sourmash and YACHT databases must be built from the same references with the same k-mer size and scaled value
(see [Databases](../getting-started/databases.md)). Several sourmash databases can be given comma-separated;
YACHT is then skipped.
:::

This step includes taxonomic assignment and detection of antimicrobial resistance genes (ARGs) from the preprocessed reads. Users must provide the `--enable_readbase` parameter to enable this step.

## Taxonomic classification of trimmed reads

### sourmash
[sourmash](https://sourmash.readthedocs.io/en/latest/) is a tool for computing and comparing MinHash sketches of DNA sequences. It is used in this pipeline for taxonomic classification of trimmed reads. This significantly expands the range of operations that can be done quickly and in low memory.

Currently, the pipeline does not automatically download the sourmash database. Users must obtain pre-built databases from [ctbrown’s farm](https://farm.cse.ucdavis.edu/~ctbrown/sourmash-db/) or build their own following the instructions [here](../getting-started/databases.md#sourmash), and provide the `--sourmash_database` and `--sourmash_ksize` parameters corresponding to the database path and its k-mer length.

:::{dropdown} Output files
- `Readbased_Analsysis/Sourmash_YACHT/`
  - `reads_collected/sourmash_csv_collected_mag.csv`: Path to the reads of `[sample_id]`, located in the Nextflow `/work/` cache folder. The CSV header is `name,genome_filename,protein_filename`.
  - `sketches/batch.batch.manysketch.zip`: Result of running `sourmash scripts manysketch`.
  - `sketches/single_sketches/manysketch_output/zip_files/[sample_id]_[ksize].sig.zip`: Result for `[sample_id]` and `[ksize]` extracted from `batch.batch.manysketch.zip`. This profile can be adjusted via the `--sourmash_profile` parameter.
  - `fastmultigather/batch.batch.manysketch_sourmash_gather.csv`: Result of running `sourmash fastmultigather` on `batch.batch.manysketch.zip` with the database provided in `--sourmash_database`. This file contains all raw information returned by `sourmash fastmultigather` for all `[sample_id]`.
  - `taxannotate/batch.batch.manysketch_sourmash_gather.with-lineages.csv.gz`: Result of running `sourmash taxannotate` on `batch.batch.manysketch_sourmash_gather.csv` with the database provided in `--sourmash_taxonomy_csv`. This file contains all raw information returned by `sourmash taxannotate` for all `[sample_id]`.
  - `taxmetagenome/batch.summarized.csv`: Result of running `sourmash taxmetagenome` on `batch.batch.manysketch_sourmash_gather.with-lineages.csv.gz`. This file contains aggregated taxonomic abundance information for all `[sample_id]`.
:::

### YACHT
[YACHT](https://github.com/KoslickiLab/YACHT) is a mathematically rigorous hypothesis test for the presence or absence of organisms in a metagenomic sample, based on average nucleotide identity (ANI).

Currently, the pipeline does not automatically download the YACHT database. Users must obtain pre-built databases or build their own following the instructions [here](../getting-started/databases.md#yacht), and provide the `--yacht_database` and `--yacht_ksize` parameters corresponding to the database path and its k-mer length. For this sub-workflow to function, users must ensure the YACHT and sourmash reference databases are identical (same k-mer length and reference genomes).

:::{dropdown} Output files
- `Readbased_Analsysis/Sourmash_YACHT/`
  - `yacht_results/[sample_id]_yacht.xlsx`: YACHT result for `[sample_id]`.
  - `yacht_results/temp_yacht_results/csv/[sample_id]_min_coverage[min_coveragevalue].csv`: Result from `[sample_id]_yacht.xlsx` split by `[min_coveragevalue]` defined in the `--min_coverage_list` parameter of the YACHT module.
  - `yacht_results/temp_yacht_results/yacht_mergedresults.{xlsx,csv}`: Merged YACHT results by `[sample_id]`.
  - `yacht_results/temp_yacht_results/yacht_merged_processed.csv`: Merged YACHT results by `[sample_id]` and `[min_coveragevalue]`.
  - `final_results/merged_sourmash_yacht.{xlsx,csv}`: Raw merged results from sourmash and YACHT.
  - `final_results/merged_sourmash_yacht_removeuncl_mincoverage[coveragevalue].{xlsx,csv}`: Merged results from sourmash and YACHT, processed with thresholds `--remove_unclassified --min_coverage 0.05` in the `PROCESS_READBASED_RESULTS` module.
:::

## Antimicrobial resistance genes (ARGs) read-based analysis

### RGI_BWT
[rgi_bwt](https://github.com/arpcard/rgi/blob/master/docs/rgi_bwt.rst) tool is a specialized module within the Resistance Gene Identifier (RGI) suite designed to analyze metagenomic reads in FASTQ format by aligning them against the Comprehensive Antibiotic Resistance Database (CARD). It supports popular aligners such as KMA (default), Bowtie2, and BWA to identify antibiotic resistance genes across various model types, including protein homologs, variants, and rRNA mutations. Users can further expand the search to include CARD’s "WildCARD" dataset to capture broader allelic diversity from environmental and clinical samples. The tool generates comprehensive tab-delimited reports at both the allele and gene levels, providing detailed metrics on read mapping, percentage coverage, depth, and SNPs to facilitate in-depth resistome profiling.

:::{dropdown} Output files
- `Readbased_Analsysis/rgi_bwt/[sample_id]/`
  - `[sample_id].allele_mapping_data.{json,txt}`: RGI bwt read mapping results at allele level.
  - `[sample_id].gene_mapping_data.txt`:RGI bwt read mapping results at gene level.
  - `[sample_id].artifacts_mapping_stats.txt`: Statistics for read mapping artifacts.
  - `[sample_id].overall_mapping_stats.txt`: Statistics for overall read mapping results.
  - `[sample_id].reference_mapping_stats.txt`: Statistics for reference matches.
:::

If `--rgi_preparecarddb_dir` is specified, the output directory will also contain the prepared database.

:::{dropdown} Output files
- `Databases/RGI/`
:::

## Read-based post-processing

With `--readbased_postprocess all|phyloseq|taxburst|rgi_bwt`, the sourmash + YACHT results are converted into ready-to-use tables (`bin/py_scripts/read_based/post_processing.py`). Filters are set with `--readbased_postprocess_options` (default `--min_coverage 0.05`). `[prefix]` is the last folder name of `--outdir`.

:::{dropdown} Output files
- `Readbased_Analysis/Post_processing/` (batch mode)
  - `[prefix]_default/[prefix].processed_microcosm.csv`: one standard table (taxonomy per rank, relative abundance per sample); written by `default` and `all`.
  - `[prefix]_phyloseq/{otu_table,tax_table,sample_data}.csv`: tables for [`phyloseq`](https://joey711.github.io/phyloseq/) (`otu_table()`, `tax_table()`, `sample_data()`).
  - `[prefix]_taxburst/*.taxburst.csv`, `*.taxburst.split.csv`: input tables for `taxburst` charts, merged and per sample.
  - `[prefix]_rgi_bwt/`: RGI bwt gene tables collected over the samples (needs `--enable_rgi_bwt`).
- `Readbased_Analysis/Post_processing/[sample_id]/{default,phyloseq,taxburst,rgi_bwt}/` (single-sketch mode, `--enable_singlesketch`), plus the all-sample aggregates.
:::

## sylph read-based profiling

[sylph](https://github.com/bluenote-1577/sylph) profiles the paired reads of each sample directly against a sylph database (`--readbased_taxotools sylph` or `microcosm,sylph`). With `--enable_sylphtax`, [sylph-tax](https://github.com/bluenote-1577/sylph-tax) adds the taxonomy.

:::{dropdown} Output files
- `Readbased_Analysis/sylph/`
  - `profile/[sample_id].tsv`: `sylph profile` output, one row per detected genome (taxonomic and sequence abundance, adjusted ANI, effective coverage).
  - `sylphtax/[sample_id].sylphmpa`: `sylph-tax taxprof` profile (MetaPhlAn-like format; `--enable_sylphtax`).
  - `[prefix]_sylphtax_[datatype].tsv`: clade × sample table merged over the samples (`--sylph_datatype`).
  - `phyloseq/{otu_table,tax_table,sample_data}.csv`: phyloseq tables (`--readbased_postprocess phyloseq|all` with `--enable_sylphtax`).
:::
