# Read preprocessing

Every route that starts from reads first checks, trims and decontaminates them.

These steps trim away the adapter sequences present in input reads, remove low-quality bases, and discard reads that are too short. It also removes host contaminants. FastQC is run for visualizing the general quality metrics of the sequencing runs before and after trimming.

## FastQC
The pipeline runs FastQC on raw reads (`Raw_reads`) and processed files at each step (`Remove_HostGenome`, `Trimming`) to assess input and output data quality during preprocessing. The pipeline performs trimming and adapter removal first, followed by host genome removal.

:::{dropdown} Output files
- `Preprocess/QC/`
  - `/Raw_reads/[sample_id]/[sample_id]_[1/2]_fastqc.html`: FastQC report, containing quality metrics for your untrimmed raw fastq files
  - `/Remove_HostGenome/[sample_id]/[sample_id]_[1/2]_fastqc.html`: FastQC report, containing quality metrics for your host-removed fastq files.
  - `/Trimming/[sample_id]/[sample_id]_[1/2]_fastqc.html`: FastQC report, containing quality metrics for your trimmed fastq files.
:::

[FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) gives general quality metrics about your sequenced reads. It provides information about the quality score distribution across your reads, per base sequence content (%A/T/G/C), adapter contamination and overrepresented sequences. For further reading and documentation see the [FastQC help pages](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/).

## fastp

[fastp](https://github.com/OpenGene/fastp) is an all-in-one fastq preprocessor for read/adapter trimming and quality control. It is used in this pipeline for trimming adapter sequences and discarding low-quality reads. Its output is in the results folder and part of the MultiQC report.

:::{dropdown} Output files
- `Preprocess/Trimming/[sample_id]/`
  - `report/[sample_id].fastp.html`: Interactive report
  - `report/[sample_id].fastp.json`: Report in json format
  - `[sample_id]*.fastq.gz`: Trimmed reads in FASTQ format (automatically saved by default).
:::

## Host read removal (hostile)

[hostile](https://github.com/bede/hostile) is a tool for removing host sequences from genomic data, improving precision by removing fewer microbial reads than current methods. It supports single or paired FASTQ input, downloads a human reference genome upon first use, and offers privacy through integer read headers. It operates efficiently, requiring 4GB of RAM for short reads and 13GB for long reads. More details and benchmarks are available in the associated paper. Users must provide the `--hostile_ref_name` parameter matching the reference genome name in `--hostile_reference` (if using a custom database) or `--hostile_index` (if allowing the pipeline to automatically select the database).

:::{dropdown} Output files
- `Preprocess/Remove_HostGenome/[sample_id]/`
  - `cleaned_reads/[sample_id]_[1/2].fastp.clean_[1/2].fastq.gz`: Trimmed and host removed reads in FASTQ format.
  - `report/[sample_id].json`: Contains information about the executed processes and the number of removed reads.
:::
