# Assembly and contig analysis

Without `--enable_readbase`, microcosm assembles each sample, maps its reads back **once** and characterises the contigs.

| Step | Default | Change with |
| :--- | :--- | :--- |
| Assembly | MEGAHIT and metaSPAdes | `--skip_megahit`, `--skip_spades`, `--assembly_input` (pre-computed assemblies) |
| Assembly QC | QUAST | `--skip_quast` |
| Read mapping (shared by contig coverage and binning) | Bowtie2 | `--binprepare_mappertool bwamem2` |
| Contig statistics and coverage | SeqKit, samtools coverage | `--skip_contig_analysis` |
| Gene prediction | Pyrodigal | `--annotation_tool prokka`, `--skip_contig_annotation` |
| Plasmids and viruses | geNomad | `--skip_genomad_e2e`, `--genomad_db` |
| Domain (bacteria / archaea / eukarya) | Tiara | `--skip_tiara_contigqc` |
| Closest reference per contig | skani | `--skip_skani_search false --skani_db <db>` |

## Assembly

Trimmed (short) reads are assembled with both megahit and SPAdes.

### MEGAHIT

[MEGAHIT](https://github.com/voutcn/megahit) is a single node assembler for large and complex metagenomics short reads.

:::{dropdown} Output files
- `Assembly/MEGAHIT/`
  - `[sample_id].contigs.fa.gz`: Compressed metagenome assembly in fasta format
  - `[sample_id].log`: Log file
  - `QC/[sample_id]/QUAST`: Directory containing QUAST files.
:::

### SPAdes

[SPAdes](http://cab.spbu.ru/software/spades/) was originally a single genome assembler that later added support for assembling metagenomes.

:::{dropdown} Output files
- `Assembly/SPAdes/`
  - `[sample_id].scaffolds.fa.gz`: Compressed assembled scaffolds in fasta format
  - `[sample_id].assembly.gfa.gz`: Compressed assembly graph in gfa format
  - `[sample_id].contigs.fa.gz`: Compressed assembled contigs in fasta format
  - `[sample_id].spades.log`: Log file
  - `QC/[sample_id]/QUAST`: Directory containing QUAST files.
:::

### Metagenome QC with QUAST

[QUAST](http://cab.spbu.ru/software/quast/) is a tool that evaluates metagenome assemblies by computing various metrics. The QUAST output is also included in the MultiQC report, as well as in the assembly directories themselves.

:::{dropdown} Output files
- `Assembly/[assembler]/QC/[sample_id]/QUAST/`
  - `report.*`: QUAST report in various formats, such as html, pdf, tex, tsv, or txt
  - `transposed_report.*`: QUAST report that has been transposed into wide format (tex, tsv, or txt)
  - `quast.log`: QUAST log file
  - `metaquast.log`: MetaQUAST log file
  - `icarus.html`: Icarus main menu with links to interactive viewers
  - `icarus_viewers/contig_size_viewer.html`: Diagram of contigs that are ordered from longest to shortest
  - `basic_stats/cumulative_plot.pdf`: Shows the growth of contig lengths (contigs are ordered from largest to shortest)
  - `basic_stats/GC_content_plot.pdf`: Shows the distribution of GC content in the contigs
  - `basic_stats/[assembler]-[sample_id]_GC_content_plot.pdf`: Histogram of the GC percentage for the contigs
  - `basic_stats/Nx_plot.pdf`: Plot of Nx values as x varies from 0 to 100%.
  - `predicted_genes/[assembler]-[sample_id].rna.gff`: Contig positions for rRNA genes in gff version 3 format
  - `predicted_genes/barrnap.log`: Barrnap log file (ribosomal RNA predictor)
:::

## Contigs analysis

### Protein-coding gene prediction
Protein-coding genes are predicted for each assembly.

#### Pyrodigal

[Pyrodigal](https://github.com/althonos/pyrodigal) annotates whole (meta-)genomes by identifying ORFs in a set of genomic DNA sequences. It produces the same results as [Prodigal](https://github.com/hyattpd/Prodigal) while being more resource-optimized, thus faster. Unlike Prodigal, Pyrodigal cannot produce output in GenBank format.

:::{dropdown} Output files
- `ContigsAnalysis/annotation/pyrodigal/`
  - `[assembler]-[sample_id].pyrodigal.gbk.gz`: Annotation in GBK format, containing both sequences and annotations.
  - `[assembler]-[sample_id].pyrodigal.faa.gz`: The protein translation file consists of all the proteins from all the sequences in multiple FASTA format.
  - `[assembler]-[sample_id].pyrodigal.fna.gz`: Nucleotide sequences of the predicted proteins using the DNA alphabet, not mRNA (so you will see 'T' in the output and not 'U').
  - `[assembler]-[sample_id].pyrodigal.score.gz`: Information about all potential genes (with scores).
:::

#### Prokka

[Prokka](https://github.com/tseemann/prokka) performs whole genome annotation to identify features of interest in a set of (meta-)genomic DNA sequences.

:::{dropdown} Output files
- `ContigsAnalysis/annotation/prokka/`
  - `[assembler]-[sample_id].contigs.prokka.gff`: annotation in GFF3 format, containing both sequences and annotations.
  - `[assembler]-[sample_id].contigs.prokka.gbk`: standard Genbank file derived from the master `.gff`.
  - `[assembler]-[sample_id].contigs.prokka.faa`: protein FASTA file of the translated CDS sequences.
  - `[assembler]-[sample_id].contigs.prokka.fna`: nucleotide FASTA file of the input contig sequences.
  - `[assembler]-[sample_id].contigs.prokka.ffn`: nucleotide FASTA file of all the prediction transcripts (CDS, rRNA, tRNA, tmRNA, misc_RNA).
  - `[assembler]-[sample_id].contigs.prokka.sqn`: Sequin file for GenBank submission.
  - `[assembler]-[sample_id].contigs.prokka.fsa`: FASTA file of the input contig sequences, used by “tbl2asn” to create the .sqn file.
  - `[assembler]-[sample_id].contigs.prokka.tbl`: Table file used by “tbl2asn” to create the .sqn file.
  - `[assembler]-[sample_id].contigs.prokka.log`: log file.
  - `[assembler]-[sample_id].contigs.prokka.txt`: statistics relating to the annotated features found.
  - `[assembler]-[sample_id].contigs.prokka.err`: unacceptable annotations - the NCBI discrepancy report.
  - `[assembler]-[sample_id].contigs.prokka.tsv`: tab-delimited feature table.
:::

### Contigs classification 

#### geNomad - Mobile genetic elements (MGEs) identification

[geNomad](https://github.com/apcamargo/genomad) identifies viruses and plasmids in sequencing data (isolates, metagenomes, and metatranscriptomes)

:::{dropdown} Output files
- `ContigsAnalysis/geNomad/[assembler]-[sample_id]*/`
  - `[assembler]-[sample_id]*_annotate`
    - `[assembler]-[sample_id]*_taxonomy.tsv`: Taxonomic assignment data
  - `[assembler]-[sample_id]*_aggregated_classification`
    - `[assembler]-[sample_id]*_aggregated_classification.tsv`: Sequence classification in tabular format
  - `[assembler]-[sample_id]*_find_proviruses`
    - `[assembler]-[sample_id]*_provirus.tsv`: Characteristics of proviruses identified by geNomad
  - `[assembler]-[sample_id]*_summary`
    - `[assembler]-[sample_id]*_virus_summary.tsv`: Virus classification summary file in tabular format
    - `[assembler]-[sample_id]*_plasmid_summary.tsv`: Plasmid classification summary file in tabular format
    - `[assembler]-[sample_id]*_viruses_genes.tsv`: Virus gene annotation data in tabular format
    - `[assembler]-[sample_id]*_plasmids_genes.tsv`: Plasmid gene annotation data in tabular format
    - `[assembler]-[sample_id]*_viruses.fna`: Virus nucleotide sequences in FASTA format
    - `[assembler]-[sample_id]*_plasmids.fna`: Plasmid nucleotide sequences in FASTA format
    - `[assembler]-[sample_id]*_viruses_proteins.faa`: Virus protein sequences in FASTA format
    - `[assembler]-[sample_id]*_plasmids_proteins.faa`: Plasmid protein sequences in FASTA format
  - `[assembler]-[sample_id]*.log`: Plain text log file detailing the steps executed by geNomad (annotate, find-proviruses, marker-classification, nn-classification, aggregated-classification and summary)
:::

#### skani - ANI and aligned fraction
[skani](https://github.com/bluenote-1577/skani) is a high-performance tool for calculating average nucleotide identity (ANI) and aligned fraction (AF) among DNA sequences such as contigs, MAGs, and complete genomes, optimized for comparisons above ~80% ANI. It employs an approximate mapping strategy without base-level alignment, achieving accuracy comparable to BLAST-based methods while being orders of magnitude faster. skani is particularly reliable for incomplete and medium-quality MAGs, reports genome aligned fractions, and supports rapid, memory-efficient searches against large genome databases.

:::{dropdown} Output files
- `ContigsAnalysis/skani/`
  - `[assembler]-[sample_id].skani.tsv`: TSV file containing ANI and AF as predicted by skani.
:::

#### Tiara - Domain classification

Tiara is a contig classifier that identifies the domain (prokarya, eukarya) of contigs within an assembly. This is used in this pipeline to rapidly and with few resources identify the most likely domain classification of each bin or unbin based on its contig identities.

:::{dropdown} Output files
- `Taxonomy/Tiara/`
  - `[assembler]-[sample_id].tiara.txt` - Tiara output classifications (with probabilities) for all contigs within the specified sample/group assembly
  - `log_[assembler]-[sample_id].txt` - log file detailing the parameters used by the Tiara model for contig classification.
- `Binning/Tiara/tiara_summary.tsv` - Summary of Tiara domain classification for all bins.
:::

Typically, you would use `tiara_summary.tsv` as the primary file to see which bins or unbins have been classified to which domains at a glance, whereas the files in `Taxonomy/Tiara` provide classifications for each contig.

### Contigs evaluations

#### seqkit_fx2tab
[seqkit_fx2tab](https://github.com/shenwei356/seqkit) is a cross-platform and ultrafast toolkit for FASTA/Q file manipulation. This option is convert FASTA/Q to tabular format, and provide various information, like sequence length, GC content/GC skew.

:::{dropdown} Output files
- `ContigsAnalysis/contig_stats/`
  - `[assembler]-[sample_id].seqkit.txt`: Text file in tabular format.
:::

#### samtools_coverage
[samtools_coverage](http://www.htslib.org/doc/samtools-coverage.html) computes coverage statistics for contigs based on read alignments in BAM format.

:::{dropdown} Output files
- `ContigsAnalysis/contig_stats/`
  - `[assembler]-[sample_id].coverage.txt`: Text file in tabular format.
:::

#### Mapping reads against the assemblies
The trimmed, host-removed reads of each sample are mapped **once** to its assembly with [Bowtie2](https://github.com/BenLangmead/bowtie2) (default) or [BWA-MEM2](https://github.com/bwa-mem2/bwa-mem2) (`--binprepare_mappertool bwamem2`). The same sorted BAM is used for the contig coverage above and for the binning depths below, so the BAM files are written only once.

:::{dropdown} Output files
- `Binning/mapping/sorted_bam/[assembler]/[sample_id]/` (`--save_assembly_mapped_reads`, default `true`)
  - `[assembler]-[sample_id].sorted.bam`: reads mapped to the assembly, sorted by position.
  - `[assembler]-[sample_id].sorted.bam.bai`: BAM index.
- `Binning/mapping/bowtie2/log/[assembler]/[sample_id]/*.log`: Bowtie2 mapping logs (also in the MultiQC report).
:::
