---
myst:
  html_meta:
    description: microcosm, a Nextflow pipeline for short-read metagenomics, from reads to typed bacterial genomes.
---

# microcosm

```{raw} html
<div class="mc-hero">
  <img class="mc-hero-logo mc-only-light" src="_static/microcosm_logo.png" alt="microcosm"/>
  <img class="mc-hero-logo mc-only-dark" src="_static/microcosm_logo_dark.png" alt="microcosm"/>
  <p class="mc-hero-tagline">Short-read metagenomics, from reads to typed bacterial genomes</p>
  <p class="mc-hero-text">A modular <a href="https://www.nextflow.io/">Nextflow</a> pipeline for Illumina shotgun
  metagenomes: fast read-based profiling, genome-resolved analysis, bacterial typing and AMR screening, and one
  MultiQC report for the whole run. Built for molecular epidemiology, tuned for modest hardware.</p>
  <div class="mc-hero-buttons">
    <a class="mc-btn mc-btn-primary" href="getting-started/quickstart/">Quick start</a>
    <a class="mc-btn" href="usage/builder/">Build my command</a>
    <a class="mc-btn" href="pipeline/overview/">How it works</a>
    <a class="mc-btn" href="https://github.com/tnmquann/microcosm">GitHub</a>
  </div>
  <div class="mc-badges">
    <span class="mc-badge">v2.0.0</span><span class="mc-badge">Nextflow ≥ 26.04</span>
    <span class="mc-badge">Singularity · Apptainer · Docker · Conda</span><span class="mc-badge">MIT license</span>
  </div>
</div>
```

:::{note}
**microcosm 2.0.0** is the first release under the new name (formerly *metaflow*). It adds sylph profiling,
read-based post-processing, a full bacterial typing and AMR layer, GTDB-Tk and one MultiQC report.
Coming from 1.x? Read [Upgrading from metaflow 1.x](reference/troubleshooting.md#upgrading-from-metaflow-1x).
:::

## What microcosm does

::::{grid} 1 2 2 2
:gutter: 3

:::{grid-item-card} {octicon}`telescope;1.2em;mc-icon` Read-based profiling
:link: pipeline/read-based
:link-type: doc
:class-card: mc-card

**sourmash + YACHT** (default) or **sylph** (+ sylph-tax), or both. AMR genes on reads with **RGI bwt**.
Ready-to-use **phyloseq** and **taxburst** tables.
:::

:::{grid-item-card} {octicon}`stack;1.2em;mc-icon` Genome-resolved metagenomics
:link: pipeline/assembly
:link-type: doc
:class-card: mc-card

**MEGAHIT / metaSPAdes**, contig analysis (geNomad, Tiara, skani), six binners with **Binette / DAS Tool**
refinement, **CheckM2, CheckM, BUSCO, GUNC**, sourmash or **GTDB-Tk** taxonomy.
:::

:::{grid-item-card} {octicon}`shield-check;1.2em;mc-icon` Bacterial genome analysis
:link: pipeline/bacterial
:link-type: doc
:class-card: mc-card

MLST, VirulenceFinder, PlasmidFinder; **AMRFinderPlus, ABRicate, RGI, ResFinder, DeepARG** harmonised with
hAMRonization + ARGnorm; **Kleborate**; **SISTR, SeqSero2, GenoTyphi, Paratype**. On MAGs or your own genomes.
:::

:::{grid-item-card} {octicon}`graph;1.2em;mc-icon` One report per run
:link: pipeline/report
:link-type: doc
:class-card: mc-card

A **MultiQC** report from read QC to typed genomes, with every tool version and a summary of the parameters you
changed.
:::
::::

```{image} ../images/microcosm_currentbeta_cleaned_metromap_light.png
:alt: microcosm workflow
:class: mc-metromap
:align: center
```

## Three routes, one pipeline

| Route | Use it to | Switch | Start from |
| :--- | :--- | :--- | :--- |
| **Read-based** | profile the community and find AMR genes in the reads | `--enable_readbase` | paired FASTQ |
| **Assembly-based** | recover and characterise genomes (MAGs), optionally type the bacterial ones | default (+ `--enable_bacterial_downstream`) | paired FASTQ or assemblies |
| **Genome input** | type and screen genome FASTA files (isolates, external MAGs) | `--genome_input` | FASTA |

## Where to start

::::{grid} 1 2 3 3
:gutter: 2

:::{grid-item-card} New to microcosm
:link: getting-started/quickstart
:link-type: doc
:class-card: mc-card-small

Six steps from installation to your first report.
:::

:::{grid-item-card} Setting up a server
:link: getting-started/installation
:link-type: doc
:class-card: mc-card-small

Nextflow, software profiles and [databases](getting-started/databases.md).
:::

:::{grid-item-card} A recipe for my goal
:link: usage/cookbook
:link-type: doc
:class-card: mc-card-small

Copy-paste commands for each analysis.
:::

:::{grid-item-card} Every parameter
:link: reference/parameters
:link-type: doc
:class-card: mc-card-small

Generated from the pipeline schema.
:::

:::{grid-item-card} Every output file
:link: reference/outputs
:link-type: doc
:class-card: mc-card-small

Folders at a glance and file descriptions.
:::

:::{grid-item-card} Something went wrong
:link: reference/troubleshooting
:link-type: doc
:class-card: mc-card-small

FAQ, common errors, upgrading from 1.x.
:::
::::

## Quick taste

```bash
nextflow run tnmquann/microcosm -r 2.0.0 -profile singularity \
    --input samples.csv --outdir project_A --enable_readbase \
    --sourmash_database /db/gtdb-rs226-reps.k31.rocksdb --sourmash_taxonomy_csv /db/gtdb-rs226.lineages.csv \
    --yacht_database /db/yacht_k31.json --enable_rgi_bwt --readbased_postprocess all
```

Then open `project_A/MultiQC/project_A_multiqc_report.html`.

```{toctree}
:caption: Getting started
:hidden:

getting-started/quickstart
getting-started/installation
getting-started/databases
```

```{toctree}
:caption: Using microcosm
:hidden:

usage/builder
usage/inputs
usage/running
usage/cookbook
```

```{toctree}
:caption: Pipeline
:hidden:

pipeline/overview
pipeline/preprocessing
pipeline/read-based
pipeline/assembly
pipeline/binning
pipeline/bacterial
pipeline/report
```

```{toctree}
:caption: Reference
:hidden:

reference/parameters
reference/outputs
reference/troubleshooting
```

```{toctree}
:caption: About
:hidden:

about/changelog
about/roadmap
about/citation
```
