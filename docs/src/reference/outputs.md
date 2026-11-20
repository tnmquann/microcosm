# Outputs

All paths are relative to `--outdir`. Each pipeline page describes the files of its steps: [preprocessing](../pipeline/preprocessing.md), [read-based](../pipeline/read-based.md), [assembly](../pipeline/assembly.md), [binning](../pipeline/binning.md), [bacterial genomes](../pipeline/bacterial.md), [MultiQC](../pipeline/report.md).

## Output folders at a glance

Click a route to see its folders.

:::{dropdown} Read-based subworkflow
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
:::

:::{dropdown} Assembly-based subworkflow
```text
outdir/
├── Databases/
│   ├── genomad/
│   ├── Bakta/
│   ├── CheckM/
│   ├── CheckM2/
│   ├── GUNC/
│   ├── BUSCO/
│   ├── AMRFinderPlus/   # bacterial downstream downloads (--save_<tool>_db)
│   ├── RGI/
│   ├── ResFinder/
│   ├── DeepARG/
│   ├── VirulenceFinder/
│   └── ParaType/
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
│   └── skani/
├── Taxonomy/
│   ├── sourmash/
│   ├── GTDB-Tk/         (--bin_classifier gtdbtk)
│   └── Tiara/
├── Bacterial_Analysis/  (--genome_input or --enable_bacterial_downstream; see above)
├── MultiQC/
└── pipeline_info/
```
:::

:::{dropdown} Bacterial genome analysis
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
:::

### Where are the main results?

| Result | Folder |
| :--- | :--- |
| **Report of the whole run** | `MultiQC/<prefix>_multiqc_report.html` |
| Taxonomic profiles (sourmash + YACHT) | `Readbased_Analysis/Sourmash-YACHT/final_results/` |
| Taxonomic profiles (sylph) | `Readbased_Analysis/sylph/` |
| phyloseq / taxburst / standard tables | `Readbased_Analysis/Post_processing/` |
| AMR genes in reads | `Readbased_Analysis/rgi_bwt/` |
| Assemblies and their QC | `Assembly/<assembler>/` |
| Contig annotation, plasmids/viruses, domains | `ContigsAnalysis/` |
| Bins (MAGs) and their QC | `Binning/<binner or refiner>/`, `Binning/QC/` |
| Bin taxonomy | `Taxonomy/` |
| Typing and AMR of bacterial genomes | `Bacterial_Analysis/` (one `summary/<prefix>_<tool>.tsv` per tool) |
| Tool versions, execution report | `pipeline_info/` |
| Downloaded databases (for reuse) | `Databases/` |

`<prefix>` is the last folder name of `--outdir`.
