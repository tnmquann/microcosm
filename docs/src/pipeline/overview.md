# Pipeline overview

```{image} ../../images/microcosm_currentbeta_cleaned_metromap_light.png
:alt: microcosm workflow
:class: mc-metromap
:align: center
```

microcosm has three routes that share the same read preprocessing and the same report. A run takes one route.

```{mermaid}
flowchart TD
  R[Paired FASTQ] --> P[Preprocessing<br/>FastQC · fastp · hostile]
  P -->|--enable_readbase| RB[Read-based profiling<br/>sourmash + YACHT · sylph · RGI bwt]
  RB --> PP[Post-processing<br/>standard · phyloseq · taxburst]
  P -->|default| A[Assembly<br/>MEGAHIT · metaSPAdes]
  AI[Pre-computed assemblies<br/>--assembly_input] --> M
  A --> M[Read mapping, once<br/>Bowtie2 / BWA-MEM2]
  M --> C[Contig analysis<br/>Pyrodigal/Prokka · geNomad · Tiara · skani]
  M --> B[Binning + refinement<br/>MetaBAT2 · MaxBin2 · COMEBin · SemiBin2 · Binette]
  B --> Q[Bin QC, annotation, taxonomy<br/>CheckM2 · BUSCO · GUNC · Prokka/Bakta · sourmash/GTDB-Tk]
  Q -->|--enable_bacterial_downstream| BA
  G[Genome FASTA<br/>--genome_input] --> BA[Bacterial genome analysis<br/>typing · AMR · Klebsiella · Salmonella]
  PP --> MQ[(MultiQC report)]
  C --> MQ
  Q --> MQ
  BA --> MQ
```

## Stages

| # | Stage | Tools | Route | Page |
| :-: | :--- | :--- | :--- | :--- |
| 1 | **Read QC and decontamination** | [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/), [fastp](https://github.com/OpenGene/fastp), [hostile](https://github.com/bede/hostile) | read, assembly | [Preprocessing](preprocessing.md) |
| 2 | **Read-based profiling** | [sourmash](https://github.com/sourmash-bio/sourmash) + [YACHT](https://github.com/KoslickiLab/YACHT), [sylph](https://github.com/bluenote-1577/sylph) + [sylph-tax](https://github.com/bluenote-1577/sylph-tax), [RGI bwt](https://github.com/arpcard/rgi/blob/master/docs/rgi_bwt.rst) | read | [Read-based](read-based.md) |
| 3 | **Assembly and contigs** | [MEGAHIT](https://github.com/voutcn/megahit), [metaSPAdes](https://github.com/ablab/spades), [QUAST](https://github.com/ablab/quast), [SeqKit](https://github.com/shenwei356/seqkit), [samtools](https://github.com/samtools/samtools), [Pyrodigal](https://github.com/althonos/pyrodigal) / [Prokka](https://github.com/tseemann/prokka), [geNomad](https://github.com/apcamargo/genomad), [Tiara](https://github.com/ibe-uw/tiara), [skani](https://github.com/bluenote-1577/skani) | assembly | [Assembly](assembly.md) |
| 4 | **Read mapping** (once per assembly) | [Bowtie2](https://github.com/BenLangmead/bowtie2) (default) or [BWA-MEM2](https://github.com/bwa-mem2/bwa-mem2) | assembly | [Assembly](assembly.md) |
| 5 | **Binning and refinement** | [MetaBAT2](https://bitbucket.org/berkeleylab/metabat), [MaxBin2](https://sourceforge.net/projects/maxbin2/), [COMEBin](https://github.com/ziyewang/COMEBin), [SemiBin2](https://github.com/BigDataBiology/SemiBin), [CONCOCT](https://github.com/BinPro/CONCOCT), [VAMB](https://github.com/RasmussenLab/vamb); [Binette](https://github.com/genotoul-bioinfo/Binette) / [DAS Tool](https://github.com/cmks/DAS_Tool) | assembly | [Binning](binning.md) |
| 6 | **Bin quality** | [CheckM2](https://github.com/chklovski/CheckM2), [CheckM](https://github.com/Ecogenomics/CheckM), [BUSCO](https://busco.ezlab.org/), [GUNC](https://github.com/grp-bork/gunc), QUAST, `assembly_stats.py`, Tiara | assembly | [Binning](binning.md) |
| 7 | **Bin annotation and taxonomy** | [Prokka](https://github.com/tseemann/prokka) / [Bakta](https://github.com/oschwengers/bakta); sourmash / [GTDB-Tk](https://github.com/Ecogenomics/GTDBTk) | assembly | [Binning](binning.md) |
| 8 | **Bacterial genome analysis** | MLST, VirulenceFinder, PlasmidFinder, AMRFinderPlus, ABRicate, RGI, ResFinder, DeepARG, hAMRonization, ARGnorm, Kleborate, SISTR, SeqSero2, Mykrobe, GenoTyphi, Paratype | assembly, genome input | [Bacterial genomes](bacterial.md) |
| 9 | **Report** | [MultiQC](https://multiqc.info/), software versions, run summary | all | [MultiQC report](report.md) |

## Design choices

- **Containers for every step.** Each process has a pinned Singularity/Apptainer/Docker image and a Conda definition,
  so results are reproducible across machines.
- **Databases are optional where possible.** A missing database is downloaded at a fixed version and saved in
  `<outdir>/Databases/` for the next run; large references (sourmash, YACHT, sylph, GTDB-Tk) stay under your control.
- **Map once.** Reads are mapped to each assembly a single time; contig coverage and binning share the BAM.
- **Every step can be switched off** (`--skip_<step>`) and every tool accepts extra arguments (`--<tool>_args`).
- **Resumable.** `-resume` reuses finished steps after a failure or a parameter change.
