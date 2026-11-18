# metaflow

[![Nextflow](https://img.shields.io/badge/version-%E2%89%A526.04.0-green?style=flat&logo=nextflow&logoColor=white&color=%230DC09D&link=https%3A%2F%2Fnextflow.io)](https://www.nextflow.io/)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with apptainer](https://img.shields.io/badge/run%20with-apptainer-1d355c.svg?labelColor=000000)](https://apptainer.org/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)
[![run with slurm](https://img.shields.io/badge/run%20with-slurm-1AAEE8.svg?labelColor=000000)](https://www.schedmd.com)
---

**metaflow** is a robust, modular pipeline built on [Nextflow](https://www.nextflow.io/) for comprehensive analysis of Illumina **short-read** metagenomic data. Designed with practical experience in mind, this workflow streamlines both rapid, read-level classification and high-resolution genome assembly, making it accessible to users with minimal computational background. It delivers fast, actionable results for professionals in molecular epidemiology, enabling timely insights for outbreak investigation and surveillance.

This pipeline was developed based on the concepts and design principles of [`nf-core/mag`](https://nf-co.re/mag/latest/) and [`SOMA`](https://github.com/ukhsa-collaboration/SOMA). Many modifications have been made to the original pipelines, reflecting both practical experience and methodological preferences. Most importantly, **metaflow** is optimized to run efficiently on resource-limited devices, ensuring broad accessibility in both laboratory and field settings.


<p align="center">
  <img src="https://github.com/tnmquann/profiler_sourmash/blob/main/docs/images/metaflow_currentbeta_cleaned_metromap_light.png" width="70%" alt="metaflow workflow diagram"/>
</p>

## Overview

This pipeline performs taxonomic profiling of shotgun metagenomics data through the following stages:
- **Read-based Analysis**  
  - Taxonomic assignment for each sequencing read is performed using an optimized combination of [`sourmash`](https://github.com/sourmash-bio/sourmash) and [`YACHT`](https://github.com/KoslickiLab/YACHT). This approach ensures fast, scalable processing, even for large datasets. Reference databases can be customized to fit your specific research needs.
  - Detection of antimicrobial resistance genes (ARGs) is integrated via [`rgi_bwt`](https://github.com/arpcard/rgi/blob/master/docs/rgi_bwt.rst), providing actionable insights into resistance profiles.

- **Metagenome Assembly & Binning**  
  - The pipeline assembles contigs and performs binning to recover Metagenome Assembled Genomes (MAGs), enabling *in silico* characterization of microbial communities.

For step-by-step installation instructions, parameter explanations, and advanced usage tutorials, please refer to [wiki page](https://github.com/tnmquann/profiler_sourmash/wiki).

## Pipeline requirements and Installation overview

### Software environment
- Nextflow version 26.04.0 or later is required as the workflow engine for metaflow v2.0.0 strict syntax compatibility.
- Conda or Mamba is used for reproducible package and environment management, with Mamba preferred for faster dependency resolution.
- Container runtimes such as Docker or Apptainer are optional but improve portability and reproducibility across heterogeneous systems.

### Hardware and System requirements
- A POSIX-compatible operating system is required, with Linux or macOS preferred; Windows is supported via WSL for compatibility.
- Memory requirements scale with the assembler: a minimum of 32 GB RAM for MEGAHIT and up to 128 GB for metaSPAdes.
- Disk usage is substantial, requiring at least 256 GB to accommodate databases, environments, intermediate files, and per-sample outputs; GPU acceleration is optional for select assembly-based steps.

### Installation workflow
- The pipeline can be obtained via Git cloning, downloading a release archive, or pulling directly with Nextflow for version-controlled deployment.
- Nextflow installation involves validating Java, downloading the executable, optionally adding it to the system PATH, and verifying functionality.
- Conda or Mamba installation finalizes the environment setup, with optional container runtime installation to enhance execution consistency.
- Execution profiles include `conda`, `docker`, and `apptainer`; GPU-capable execution is opt-in with profile combinations such as `docker,gpu` or `apptainer,gpu` on hosts that provide compatible GPU runtime support.
> [!NOTE]
> metaflow relies on `sourmash` and `YACHT` for robust taxonomic classification. For detailed setup instructions, please refer to the [Manual Database Setup](https://github.com/tnmquann/profiler_sourmash/wiki/Installation#manual-database-setup-instructions) section.

## Input specifications
*   **Sample CSV format:** `sample_id` is required and must be non-empty for every row.
    ```csv
    sample_id,run_id,group,short_reads_1,short_reads_2,long_reads
    sample1,1,0,data/sample1_1.fq.gz,data/sample1_2.fq.gz,data/sample1.fastq.gz
    sample2,0,0,data/sample2_1.fastq.gz,data/sample2_2.fastq.gz,
    ```
*   **Directory:** A directory containing paired-end FASTQ files.

## Quick start

### Read-based analysis
Once your `sourmash` and `YACHT` databases are ready, you can launch the `metaflow` pipeline. Below are examples for different input formats.
```bash
# With csv file
nextflow run main.nf \
    --input /path/to/your/samples.csv \
    --input_format csv \
    -profile conda \
    --outdir /path/to/output/directory \
    --sourmash_database /path/to/your/sourmash_database \
    --yacht_database /path/to/your/yacht_database.json \
    --enable_readbase

# With directory of FASTQ files
nextflow run main.nf \
    --input /path/to/your/fastq_directory \
    --input_format directory \
    -profile conda \
    --outdir /path/to/output/directory \
    --sourmash_database /path/to/your/sourmash_database \
    --yacht_database /path/to/your/yacht_database.json \
    --enable_readbase
```

The examples above use `-profile conda`. Use `-profile docker` or `-profile apptainer` for containerized runs, and add `gpu` only when GPU runtime support is available, for example `-profile docker,gpu`.

### Optional read-based post-processing

Add `--readbased_postprocess` to either read-based execution mode to create downstream-ready outputs from the merged Sourmash-YACHT CSV:

```bash
nextflow run main.nf \
    --input /path/to/your/samples.csv \
    --input_format csv \
    --outdir /path/to/project_A \
    --enable_readbase \
    --readbased_postprocess all \
    --postprocess_options "--min_coverage 0.05"
```

Valid modes are `all`, `phyloseq`, `taxburst`, and `rgi_bwt`. The explicit `rgi_bwt` mode requires `--enable_rgi_bwt`. With `all`, the pipeline also writes the standard default processed CSV; RGI-BWT results are included only when RGI-BWT is enabled, otherwise the pipeline warns, skips that output group, and continues with default, phyloseq, and taxburst.

### Read-based profilers: microcosm (sourmash + YACHT) and sylph

`--readbased_taxotools` chooses the read-based profiler: `microcosm` (default,
sourmash + YACHT), `sylph`, or `microcosm,sylph` to run both. sylph profiles
the paired reads directly (no read merging), one task per sample, and writes to
`Readbased_Analysis/sylph/`; the microcosm outputs stay in `Sourmash-YACHT/` and
`Post_processing/`. `--enable_singlesketch` only changes the microcosm branch.

```bash
nextflow run main.nf --input samples.csv --outdir /path/to/project_A --enable_readbase \
    --readbased_taxotools microcosm,sylph \
    --sourmash_database /path/to/gtdb-rs226-reps.k31.rocksdb --yacht_database /path/to/yacht.json \
    --sylph_database /path/to/gtdb-r226-c200-dbv2.syl2db \
    --enable_sylphtax --sylph_taxonomy /path/to/gtdb_r226_metadata.tsv.gz \
    --readbased_postprocess phyloseq
```

| Parameter | Meaning |
|---|---|
| `--sylph_database` | sylph database(s) (`*.syldb`/`*.syl2db`), comma-separated for several |
| `--sylph_args` | extra `sylph profile` arguments |
| `--enable_sylphtax` | run `sylph-tax taxprof`; needs `--sylph_taxonomy` (sylph-tax metadata file(s) matching the databases, comma-separated) |
| `--sylph_datatype` | `taxonomic_abundance` (default, sylph-tax `relative_abundance`) or `sequence_abundance` (adds `--estimate-read-counts`, so values are estimated read counts) |

With `--enable_sylphtax`, `Readbased_Analysis/sylph/` holds `profile/<sample>.tsv`,
`sylphtax/<sample>.sylphmpa` and `<prefix>_sylphtax_<datatype>.tsv` (clade x
sample). `--readbased_postprocess phyloseq` (or `all`) also writes
`sylph/phyloseq/{tax_table,otu_table,sample_data}.csv` in the same layout as the
microcosm phyloseq tables (taxa = genomes, values = `--sylph_datatype`); this
needs `--enable_sylphtax`. `taxburst` and `rgi_bwt` post-processing use the
microcosm result only. RGI-BWT (`--enable_rgi_bwt`) also runs with `sylph` alone.

Several sourmash databases: `--sourmash_database "path/to/db1.zip,path/to/db2.zip"`
(give one taxonomy per source with `--sourmash_taxonomy_csv "tax1.csv,tax2.csv"`).
fastmultigather combines zip/sig databases; several RocksDB indexes cannot be
combined. YACHT takes one database, so YACHT (and the post-processing that needs
it) is skipped when several databases are given.

`--enable_copymergedseqs`, `--enable_copysketch`, `--enable_copysketch_singles`
and `--enable_copyintermediate` (YACHT tables, collected read CSV) turn the
publishing of those intermediate files on or off.

`--postprocess_options` accepts one shell-style option string and safely tokenizes it before invoking `post_processing.py`. For example, the general form `--postprocess_options "--params_1 abc def --params_2 sss ddd"` is passed as the corresponding argument sequence. Use only options supported by the current Python script; a current valid example is `--postprocess_options "--min_coverage 0.05 --use_average"`. If one individual value contains spaces, quote that value inside the outer option string.

A MultiQC report is written to `${outdir}/MultiQC/<prefix>_multiqc_report.html`
(turn off with `--skip_multiqc`), styled by `assets/multiqc_config.yml` with the
microcosm logo. Sections: microcosm top taxa (one bar graph with a tab per rank,
Strain to Domain, from `<prefix>.processed_metaflow.csv`, so it needs
`--readbased_postprocess all`), sylph-tax (`--enable_sylphtax`), a bacterial
genome overview (species, MLST, SISTR serovar, AMR genes per tool), FastQC/fastp,
bin QC (CheckM2, CheckM, BUSCO, QUAST), Prokka, the run summary (parameters that
differ from the defaults) and the software versions. `--multiqc_config` and
`--multiqc_extra_config` are applied after the default config; `--multiqc_logo`,
`--multiqc_replace_names` and `--multiqc_sample_names` are passed to MultiQC.

Software versions of every tool used in a run are collected in one file,
`${outdir}/pipeline_info/microcosm_software_versions.yml`, next to the
execution report.

### Assembly-based analysis
For assembly-based workflows, ensure your `sourmash` database is prepared.
`--binqc_tool` accepts one or more of `busco`, `checkm`, `checkm2`
(e.g. `--binqc_tool checkm2,checkm`); each tool writes its own
`Binning/QC/<tool>_summary.tsv`.
```bash
nextflow run main.nf \
    --input /path/to/your/samples.csv \
    --input_format csv \
    -profile conda \
    --outdir /path/to/output/directory \
    --sourmash_database /path/to/your/sourmash_database
```

## Output structure
### Read-based subworkflow
The following tree illustrates the organization of output files generated by the read-based subworkflow:

```text
outdir/
├── Databases/
│   └── RGI/
├── Preprocess/
│   ├── Merged_sequences/
│   └── QC/
│       ├── Raw_reads/
│       ├── Remove_HostGenome/
│       └── Trimming/
├── Remove_HostGenome/
├── Trimming/
├── Readbased_Analysis/
│   ├── rgi_bwt/
│   ├── Post_processing/
│   │   ├── <batch-prefix>_default/
│   │   ├── <batch-prefix>_phyloseq/
│   │   ├── <batch-prefix>_taxburst/
│   │   ├── <batch-prefix>_rgi_bwt/
│   │   └── <sample_id>/
│   │       ├── default/
│   │       ├── phyloseq/
│   │       ├── taxburst/
│   │       └── rgi_bwt/
│   ├── sylph/                  (--readbased_taxotools sylph)
│   │   ├── profile/  sylphtax/  phyloseq/
│   │   └── <batch-prefix>_sylphtax_<datatype>.tsv
│   └── Sourmash-YACHT/
│       ├── final_results/
│       ├── fastmultigather/
│       ├── reads_collected/
│       ├── sketches/
│       │   ├── single_sketches/
│       │   └── batch.batch.manysketch.zip
│       ├── taxannotate/
│       ├── taxmetagenome/
│       └── yacht_results/
```

Each directory is created automatically by the pipeline and contains outputs relevant to its analysis step. The results of the ARG read-based analysis will be located in the folder `./Readbased_Analysis/rgi_bwt/`, and the results of the read-based taxonomic classification will be located in `./Readbased_Analysis/Sourmash-YACHT/final_results/`.

Post-processing outputs are always published below `${outdir}/Readbased_Analysis/Post_processing`. Batch mode keeps the existing `<batch-prefix>_<mode>` layout, where the prefix is the basename of `outdir` (`/path/to/out/` becomes `out`). With `--readbased_postprocess all`, it also writes `<batch-prefix>_default/<batch-prefix>.processed_metaflow.csv`. Single-sketch mode uses the required CSV `sample_id`, or the normalized pairing key derived from directory input filenames (`SAMPLE_001_R1.fastq.gz` and `SAMPLE_001_R2.fastq.gz` become `SAMPLE_001`), and publishes all-mode outputs as `<sample_id>/<default|phyloseq|taxburst|rgi_bwt>`; its default CSV is `<sample_id>/default/<sample_id>.processed_metaflow.csv` and its RGI-BWT directory contains only direct `*.gene_mapping_data.txt` files. It also concatenates every sample's merged read-based table and publishes combined `<batch-prefix>_default`, `<batch-prefix>_phyloseq`, and `<batch-prefix>_taxburst` sets; its combined RGI-BWT group contains those same per-sample-named mapping files directly in `<batch-prefix>_rgi_bwt/`.

### Assembly-based subworkflow

#### Bacterial genome downstream analysis

Supply `--genome_input` with a CSV containing `sample_id,fasta` and optional
`run_id,group,assembler,lineage` to analyse standalone genome FASTAs. A valid
GTDB-style supplied lineage bypasses reclassification; a blank lineage is
classified through the existing bin-classification route. Only records resolved
as `d__Bacteria` enter the downstream tools. The ordinary assembly/bin route
uses the same layer when `--enable_bacterial_downstream` is set.

Each `d__Bacteria` genome then goes through:

- General typing: MLST, VirulenceFinder and PlasmidFinder.
- AMR: AMRFinderPlus (`--plus`, with `--organism` chosen from the lineage),
  ABRicate, RGI (contig mode), ResFinder with PointFinder and DisinFinder, and
  DeepARG (model `SS`); every report is converted with hAMRonization, summarised,
  and normalised with ARGnorm (RGI is not supported by ARGnorm).
- Kleborate v3 for Klebsiella, Raoultella and Escherichia, with the preset
  chosen from the lineage (`kpsc`, `kosc`, `escherichia`, or generic modules).
- SISTR and SeqSero2 for Salmonella. The SISTR serovar alone selects follow-up:
  `Typhi` runs Mykrobe and GenoTyphi, `Paratyphi A` runs Paratype.

Database parameters (`--amrfinderplus_db`, `--rgi_preparecarddb_dir`,
`--resfinder_db`, `--pointfinder_db`, `--disinfinder_db`, `--deeparg_db`,
`--virulencefinder_db`, `--paratype_bundle`) are optional. When one is not set,
the database is downloaded into the run work directory at a fixed version
(AMRFinderPlus uses the current release). PlasmidFinder and ABRicate use the
databases bundled with their containers. Nothing is downloaded when no
bacterial genome reaches the downstream layer.

Every step can be customised with parameters (see `nextflow_schema.json`):

| What | Parameters |
|---|---|
| Turn a tool off | `--skip_mlst`, `--skip_virulencefinder`, `--skip_plasmidfinder`, `--skip_kleborate`, `--skip_sistr`, `--skip_seqsero2`, `--skip_genotyphi`, `--skip_paratype`, `--skip_amrfinderplus`, `--skip_abricate`, `--skip_rgi`, `--skip_resfinder`, `--skip_deeparg`, `--skip_argnorm` |
| Tool arguments | `--<tool>_args` (e.g. `--mlst_args`, `--kleborate_args`, `--sistr_args`). Defaults taken from the reference scripts are part of the value and are replaced when you set it: `amrfinderplus_args='--plus'`, `resfinder_args='--acquired --disinfectant'`, `rgi_main_args='--clean -g PYRODIGAL'`, `seqsero2_args='-m k -t 4'` |
| Lineage rules | `--kleborate_preset_map` (default `assets/kleborate_presets.tsv`), `--amrfinderplus_organism_map` (default `assets/amrfinderplus_organisms.tsv`): TSV files; the first `lineage_pattern` found in the GTDB lineage wins. GTDB placeholder suffixes (`_A`, `_B`, ...) are ignored when matching, on both the lineage and the pattern (`s__Enterococcus_B faecium` matches `s__Enterococcus faecium`); published lineages keep the suffixes. |
| Database versions downloaded | `--virulencefinder_db_revision`, `--resfinder_db_revision`, `--paratype_revision` |
| Other | `--deeparg_model` (`SS` or `LS`), `--deeparg_db_version`, `--abricate_db_id`, `--virulencefinder_min_coverage`, `--virulencefinder_min_identity` |

SISTR alone routes Salmonella follow-up, so `--skip_sistr` also stops
GenoTyphi and Paratype. CPUs/memory per tool are set with the usual
`process.withName` selectors in a custom config.

Results are published below `${outdir}/Assembly/Bacterial_Analysis/`, one
folder per tool and one file per genome (`<genome_id>.<tool>.*`). Folders of
tools that did not run are absent. `<prefix>` is the last component of
`--outdir`, as for the read-based post-processing tables.

```text
Bacterial_Analysis/
├── routing/          bacterial_routing_manifest.tsv, excluded_genomes.tsv, unresolved_genomes.tsv
├── amr/
│   ├── amrfinderplus/ abricate/ rgi/ resfinder/ deeparg/
│   ├── hamronization/          <genome_id>.<tool>_hamronized.tsv
│   ├── argnorm/                <genome_id>.<tool>_hamronized_argnorm.tsv
│   ├── hamronization_summary/  amr_summary.tsv, amr_summary_argnorm.tsv
│   └── summary/                <prefix>_{amrfinderplus,abricate,rgi,resfinder,pointfinder,disinfinder,deeparg}.tsv
├── general_typing/
│   ├── mlst/ virulencefinder/ plasmidfinder/
│   └── summary/                <prefix>_{mlst,virulencefinder,plasmidfinder}.tsv
├── klebsiella/
│   ├── kleborate/
│   └── summary/                <prefix>_kleborate.tsv
└── salmonella/
    ├── sistr/ seqsero2/ genotyphi/ paratype/ routing/
    └── summary/                <prefix>_{sistr,seqsero2,genotyphi,paratype,salmonella_typing}.tsv
```

Each summary table is the tool's per-genome tables stacked, with `genome_id` as
the first column.

The following tree illustrates the organization of output files generated by the assembly-based subworkflow:

```text
outdir/
├── Databases/
│   ├── genomad/
│   ├── bakta/
│   ├── CheckM/
│   ├── CheckM2/
│   ├── GUNC/
│   └── BUSCO/
├── Preprocess/
│   ├── Merged_sequences/
│   └── QC/
│       ├── Raw_reads/
│       ├── Remove_HostGenome/
│       └── Trimming/
├── Remove_HostGenome/
├── Trimming/
├── Assembly/
│   ├── MEGAHIT/
│   │   └── QC/
│   └── SPAdes/
│       └── QC/
├── Binning/
│   ├── Annotation/
│   │   ├── Prokka/
│   │   └── Bakta/
│   ├── Binette/
│   ├── DASTool/
│   ├── COMEBin/
│   ├── MaxBin2/
│   ├── MetaBAT2/
│   ├── SemiBin2/
│   ├── CONCOCT/
│   ├── VAMB/
│   ├── depths/
│   ├── mapping/
│   └── QC/
│       ├── AssemblyStats/
│       ├── CheckM/
│       ├── CheckM2/
│       ├── QUAST/
│       ├── GUNC/
│       └── BUSCO/
├── ContigsAnalysis/
│   ├── annotation/
│   │   ├── Prokka/
│   │   └── pyrodigal/
│   ├── contig_stats/
│   ├── coverage/
│   ├── genomad/
│   ├── skani/
│   └── mapped_reads/
├── Taxonomy/
│   ├── sourmash/
│   └── Tiara/
```

## Error Handling

The pipeline implements an automatic error recovery strategy:
-   Maximum of 3 retries per task.
-   Automatic resource scaling on retry.
-   Configurable error strategy.

## Author
*   Ton Ngoc Minh Quan
*   Mai Thu Si Nguyen

## Citation

If you use this pipeline, please cite:
> Minh-Quan T.N., Nguyen M.T.S. (2025). metaflow: A Nextflow pipeline for comprehensive analysis of short-read metagenomic data using sourmash and YACHT.

## Acknowledgements

This pipeline uses code and infrastructure developed and maintained by the [nf-core](https://nf-co.re) community, reused here under the [MIT license](https://github.com/nf-core/tools/blob/master/LICENSE). We gratefully acknowledge the nf-core community for their tooling, framework, and support, which have been essential to the development of this pipeline.

> The nf-core framework for community-curated bioinformatics pipelines.  
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.  
> *Nat Biotechnol.* 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://doi.org/10.1038/s41587-020-0439-x).

> Empowering bioinformatics communities with Nextflow and nf-core.  
> Langer, B.E., Amaral, A., Baudement, M.O., Bonath, F. et al.  
> *Genome Biol.* 26, 228 (2025). doi: [10.1186/s13059-025-03673-9](https://doi.org/10.1186/s13059-025-03673-9)
