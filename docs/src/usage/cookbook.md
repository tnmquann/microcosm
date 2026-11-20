# Cookbook

Recipes by goal. Each one is a complete command you can copy and adapt; replace `/db/...` with your database paths.

:::{important}
**Database dependencies**
- microcosm relies on `sourmash` and `YACHT` for robust taxonomic classification. 
- If your workflow involves **read-based analysis** (requires both sourmash & YACHT) or **assembly-based analysis with taxonomic classification** (requires sourmash), you must prepare these databases prior to execution.
- For detailed setup instructions, please refer to the [Manual Database Setup](../getting-started/databases.md) section.
:::

:::{caution}
**Compatibility Warning** - To ensure seamless integration, your `sourmash` and `YACHT` databases must be generated using identical parameters: **database type** (e.g., GTDB), **k-mer size**, and **scaled value**. Mismatched parameters will cause pipeline failures.
:::


## Read-based analysis

Once your `sourmash` and `YACHT` databases are ready, you can launch the `microcosm` pipeline. Below are examples for different input formats.

### Option 1: Using a samplesheet (CSV)
This is the recommended method for complex datasets.

```bash
nextflow run main.nf \
    --input /path/to/your/samples.csv \
    --input_format csv \
    -profile conda \
    --outdir /path/to/output/directory \
    --sourmash_database /path/to/your/sourmash_database \
    --yacht_database /path/to/your/yacht_database.json \
    --enable_readbase
```

### Option 2: Using a Directory of FASTQ files
Useful for quick runs on a folder of raw data.

```bash
nextflow run main.nf \
    --input /path/to/your/fastq_directory \
    --input_format directory \
    -profile conda \
    --outdir /path/to/output/directory \
    --sourmash_database /path/to/your/sourmash_database \
    --yacht_database /path/to/your/yacht_database.json \
    --enable_readbase
```

:::{note}
**Configuration details**
- The examples above use the `conda` profile.
- Default settings assume a **k-size of 31** for both databases.
- The pipeline is configured to perform **only** read-based taxonomic classification.
- To include read-based Antimicrobial Resistance Gene (ARG) analysis, append the `--enable_rgi_bwt` flag.
:::

:::{tip}
**Database recommendation**
For comprehensive bacterial and archaeal analysis, we recommend using the **GTDB-R226** database.
- Use `k=31` for **species-level** resolution.
- Use `k=51` for **strain-level** resolution.
:::


## Read-based profiling with sylph

sylph can replace or complement sourmash + YACHT (`--readbased_taxotools sylph` or `microcosm,sylph`):

```bash
nextflow run tnmquann/microcosm -r 2.0.0 -profile singularity \
    --input samples.csv --outdir results_sylph --enable_readbase \
    --readbased_taxotools sylph \
    --sylph_database /path/to/gtdb-r226-c200-dbv2.syl2db \
    --enable_sylphtax --sylph_taxonomy /path/to/gtdb_r226_metadata.tsv.gz \
    --readbased_postprocess phyloseq
```

- Without `--enable_sylphtax`, you get the sylph profiles (`Readbased_Analysis/sylph/profile/`) and a MultiQC section per reference genome.
- `--sylph_datatype sequence_abundance` reports estimated read counts instead of taxonomic abundance.
- `--enable_rgi_bwt` also works with sylph alone.


## Read-based post-processing (phyloseq, taxburst)

```bash
nextflow run tnmquann/microcosm -r 2.0.0 -profile singularity \
    --input samples.csv --outdir project_A --enable_readbase \
    --sourmash_database /path/to/db.rocksdb --yacht_database /path/to/yacht.json \
    --readbased_postprocess all \
    --readbased_postprocess_options "--min_coverage 0.1 --remove_unclassified"
```

`all` writes `project_A.processed_microcosm.csv`, the phyloseq tables and the taxburst tables below `Readbased_Analysis/Post_processing/` (plus RGI bwt tables with `--enable_rgi_bwt`). In R:

```r
library(phyloseq)
otu <- read.csv("project_A_phyloseq/otu_table.csv", row.names = 1, check.names = FALSE)
tax <- read.csv("project_A_phyloseq/tax_table.csv", row.names = 1)
meta <- read.csv("project_A_phyloseq/sample_data.csv", row.names = 1)
ps <- phyloseq(otu_table(as.matrix(otu), taxa_are_rows = TRUE), tax_table(as.matrix(tax)), sample_data(meta))
```


## Assembly-based analysis

For assembly-based workflows, ensure your `sourmash` database is prepared.

```bash
nextflow run main.nf \
    --input /path/to/your/samples.csv \
    --input_format csv \
    -profile conda \
    --outdir /path/to/output/directory \
    --sourmash_database /path/to/your/sourmash_database
```

:::{note}
**Workflow overview**
Running the command above (defaulting to `k=31` for sourmash) triggers a comprehensive analysis pipeline:
1.  **Preprocessing:** Quality control of raw sequences.
2.  **Assembly:** Parallel assembly using **MEGAHIT** and **metaSPAdes**.
3.  **QC:** Assembly quality assessment via **QUAST**.
4.  **Contig analysis:** Basic analysis of assembled contigs (skani is excluded).
5.  **Binning:** Multi-tool binning with **MetaBAT2**, **MaxBin2**, **COMEBin**, and **SemiBin2**.
6.  **Refinement:** Bin consolidation and refinement using **Binette**.
7.  **Bin QC:** Quality assessment with **CheckM2**.
8.  **Annotation:** Bin annotation via **Prokka**.
9.  **Classification:** Taxonomic classification of bins using **sourmash**.

*Note: Apart from the user-provided `sourmash` database, all other required reference data will be automatically downloaded by the pipeline.*
:::

:::{tip}
**Bin classification with GTDB-Tk**
To classify bins with GTDB-Tk instead of sourmash, add `--bin_classifier gtdbtk --gtdbtk_db /path/to/gtdbtk/release226` (extracted GTDB-Tk reference data; GTDB-Tk needs about 64-128 GB of memory). To turn bin classification off, add `--skip_bin_classification`.
:::


## Bacterial genome analysis

### From your own genome FASTA files

```bash
nextflow run tnmquann/microcosm -r 2.0.0 -profile singularity \
    --genome_input genomes.csv --outdir results_genomes \
    --sourmash_database /path/to/gtdb-rs226-k51.rocksdb --sourmash_ksize 51 \
    --sourmash_taxonomy_csv /path/to/gtdb-rs226.lineages.csv
```

`genomes.csv` has the columns `sample_id,fasta[,lineage]` (see [Usage](inputs.md#genome-input)). Genomes with a GTDB lineage are not reclassified; the sourmash (or GTDB-Tk) options are only needed for rows with a blank lineage. A small example is in `test_data/bacterial_genomes/` of the repository.

### From the bins of an assembly run

Add `--enable_bacterial_downstream` to an assembly-based run: every bin classified as `d__Bacteria` is analysed.

### What runs, and how to change it

| Genome | Tools |
| :--- | :--- |
| Every bacterium | MLST, VirulenceFinder, PlasmidFinder, AMRFinderPlus (`--organism` from the lineage), ABRicate, RGI, ResFinder + PointFinder + DisinFinder, DeepARG, hAMRonization, ARGnorm |
| *Klebsiella*, *Raoultella*, *Escherichia* | + Kleborate (`kpsc`, `kosc` or `escherichia` preset) |
| *Salmonella* | + SISTR, SeqSero2; Typhi → Mykrobe + GenoTyphi; Paratyphi A → Paratype |

- Turn a tool off for every genome: `--skip_<tool>` (e.g. `--skip_deeparg --skip_abricate`).
- Change a tool's options: `--<tool>_args` (e.g. `--mlst_args "--scheme ecoli_achtman_4"`). The defaults (`--plus` for AMRFinderPlus, `--acquired --disinfectant` for ResFinder, ...) are part of the value, so repeat them if you still want them.
- Databases are downloaded when no path is given and saved to `${outdir}/Databases/<Tool>/`; pass that folder in later runs (e.g. `--rgi_preparecarddb_dir results/Databases/RGI`) to skip the download.

### Per-lineage decision table

`--bacterial_decision_table` chooses the steps (and options) per lineage. Rows are matched on whole GTDB taxa at any rank; the most specific match wins. Blank cells follow the default route, and `--skip_<tool>` always wins.

```text
lineage_pattern	label	amrfinderplus	abricate	rgi	resfinder	deeparg	mlst	virulencefinder	plasmidfinder	kleborate	mlst_args	amrfinderplus_organism
p__Pseudomonadota	Pseudomonadota	y	n	n	y	n	y	n	n				
g__Escherichia	Escherichia	y	n	n	y	n	y	y	y	n	--scheme ecoli_achtman_4	Escherichia
s__Klebsiella pneumoniae	K. pneumoniae	y	n	n	y	n	y	n	n	y	--scheme klebsiella	
```

The row and steps chosen for each genome are listed in `Bacterial_Analysis/routing/bacterial_routing_manifest.tsv`. A complete example is `assets/bacterial_decision_table_example.tsv`.

### Why is a genome missing from the results?

Check `Bacterial_Analysis/routing/`:
- `excluded_genomes.tsv`: the genome is not a bacterium (e.g. an archaeon).
- `unresolved_genomes.tsv`: no lineage could be assigned; the `reason` column gives the classifier's message (for GTDB-Tk e.g. `No bacterial or archaeal marker`, usually a very incomplete or contaminated genome).


## The MultiQC report

Open `MultiQC/<prefix>_multiqc_report.html`. Use `--multiqc_replace_names` / `--multiqc_sample_names` (TSV files) to rename samples, `--multiqc_extra_config` to change the report, and `--skip_multiqc` to turn it off.
