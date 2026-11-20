# Quick start

From nothing to your first report in six steps. Each step links to the page with the details.

## 1. Install Nextflow and a software profile

- [Nextflow](https://www.nextflow.io/docs/latest/install.html) **≥ 26.04** (needs Java 17 or later): `curl -s https://get.nextflow.io | bash`
- One way to get the tools: **Singularity/Apptainer** or **Docker** (recommended, `-profile singularity`, `-profile apptainer`, `-profile docker`), or **Conda/Mamba** (`-profile conda`).

Check: `nextflow -version` and `nextflow run tnmquann/microcosm -r 2.0.0 --help`.

Details: [Installation](installation.md).

## 2. Choose what you want to run

| I want to... | Route | Main switch |
| :--- | :--- | :--- |
| know **who is in my samples** (fast, from reads) | read-based | `--enable_readbase` |
| find **AMR genes in my reads** | read-based | `--enable_readbase --enable_rgi_bwt` |
| **recover genomes (MAGs)** and check their quality | assembly-based (default) | no switch needed |
| **type bacterial genomes and find AMR genes** in those MAGs | assembly-based + bacterial analysis | `--enable_bacterial_downstream` |
| **type my own genome FASTA files** (isolates, MAGs made elsewhere) | bacterial analysis only | `--genome_input genomes.csv` |

Overview of every step: [Overview](../pipeline/overview.md).

## 3. Get the databases you need

| Route | You must download | The pipeline downloads (when no path is given) |
| :--- | :--- | :--- |
| read-based | sourmash database + lineages, YACHT database (or a sylph database) | hostile host index, CARD (RGI) |
| assembly-based | sourmash database + lineages (bin classification) | hostile, geNomad, CheckM2, Bakta, GUNC, ... |
| bacterial analysis | sourmash (or GTDB-Tk) only for genomes without a lineage | AMRFinderPlus, CARD, ResFinder, DeepARG, VirulenceFinder, Paratype |

Downloaded databases are saved in `<outdir>/Databases/` so the next run can reuse them. Details: [Installation → Database Preparation](databases.md).

## 4. Describe your samples

Reads, one row per sample (`--input samples.csv`):

```text
sample_id,run_id,group,short_reads_1,short_reads_2,long_reads
sample1,1,0,data/sample1_1.fq.gz,data/sample1_2.fq.gz,
sample2,1,0,data/sample2_1.fq.gz,data/sample2_2.fq.gz,
```

or a folder of paired FASTQ files (`--input reads/ --input_format directory`), or genome FASTA files (`--genome_input genomes.csv` with `sample_id,fasta[,lineage]`). Details: [Usage → Input Specifications](../usage/inputs.md).

## 5. Run

```bash
nextflow run tnmquann/microcosm -r 2.0.0 -profile singularity \
    --input samples.csv --outdir results \
    --enable_readbase \
    --sourmash_database /db/gtdb-rs226-reps.k31.rocksdb \
    --sourmash_taxonomy_csv /db/gtdb-rs226.lineages.csv \
    --yacht_database /db/yacht_k31.json
```

- Add `-resume` to continue a stopped run without redoing finished steps.
- Limit the machine use with `--max_cpus 16 --max_memory 64.GB`.
- More commands, one per goal: [Cookbook](../usage/cookbook.md).

## 6. Look at the results

1. Open `results/MultiQC/results_multiqc_report.html` in a browser: one report with every step that ran.
2. Go to the folder of the step you care about (`Readbased_Analysis/`, `Assembly/`, `Binning/`, `Bacterial_Analysis/`, ...). What each file contains: [Outputs](../reference/outputs.md).
3. Something odd? [Troubleshooting & FAQ](../reference/troubleshooting.md).
