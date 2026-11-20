# Bacterial genome analysis

Bacterial genomes are typed and screened for antimicrobial resistance (AMR). The genomes are the bacterial bins of an
assembly run (`--enable_bacterial_downstream`) or your own FASTA files (`--genome_input`, see [Inputs](../usage/inputs.md)).

```{mermaid}
flowchart TD
  G[Genome / MAG] --> L{GTDB lineage}
  L -->|given in CSV| R
  L -->|blank| C[sourmash or GTDB-Tk] --> R{d__Bacteria?}
  R -->|no| X[excluded / unresolved]
  R -->|yes| T[General typing + AMR]
  T --> K[Klebsiella / Escherichia: Kleborate]
  T --> S[Salmonella: SISTR, SeqSero2]
  S -->|Typhi| GT[Mykrobe + GenoTyphi]
  S -->|Paratyphi A| PT[Paratype]
```

| Genome | Tools |
| :--- | :--- |
| **Every bacterium** | MLST, VirulenceFinder, PlasmidFinder · AMRFinderPlus (organism from the lineage), ABRicate, RGI, ResFinder + PointFinder + DisinFinder, DeepARG → hAMRonization → ARGnorm |
| ***Klebsiella*, *Raoultella*, *Escherichia*** | + Kleborate v3 (`kpsc`, `kosc` or `escherichia` preset from the lineage) |
| ***Salmonella*** | + SISTR, SeqSero2; serovar Typhi → Mykrobe + GenoTyphi; Paratyphi A → Paratype |

:::{dropdown} Customising the steps
| What | How |
| :--- | :--- |
| Turn a tool off | `--skip_mlst`, `--skip_virulencefinder`, `--skip_plasmidfinder`, `--skip_kleborate`, `--skip_sistr`, `--skip_seqsero2`, `--skip_genotyphi`, `--skip_paratype`, `--skip_amrfinderplus`, `--skip_abricate`, `--skip_rgi`, `--skip_resfinder`, `--skip_deeparg`, `--skip_argnorm` |
| Tool options | `--<tool>_args`; defaults are part of the value: `amrfinderplus_args='--plus'`, `resfinder_args='--acquired --disinfectant'`, `rgi_main_args='--clean -g PYRODIGAL'`, `seqsero2_args='-m k -t 4'` |
| Lineage rules | `--kleborate_preset_map`, `--amrfinderplus_organism_map` (TSV; first matching `lineage_pattern` wins; GTDB suffixes such as `_A` are ignored) |
| Per-lineage steps | `--bacterial_decision_table` (see the [Cookbook](../usage/cookbook.md#per-lineage-decision-table)) |
| Databases | all optional; downloaded and saved to `Databases/<Tool>/` when not given ([Databases](../getting-started/databases.md#bacterial-genome-analysis-databases)) |

SISTR alone routes *Salmonella* follow-up, so `--skip_sistr` also stops GenoTyphi and Paratype.
:::

Bacterial bins (`--enable_bacterial_downstream`) or genome FASTAs (`--genome_input`) resolved as `d__Bacteria` are typed and screened for AMR. Each tool writes one file per genome (`[genome_id]` = bin name or `sample_id` of `--genome_input`) and one stacked summary per tool (`summary/[prefix]_[tool].tsv`, `genome_id` first). Folders of tools that did not run are absent.

## Routing

:::{dropdown} Output files
- `Bacterial_Analysis/routing/`
  - `bacterial_routing_manifest.tsv`: genomes sent to the downstream tools, with their lineage, the matched `--bacterial_decision_table` row and the steps/parameters chosen.
  - `excluded_genomes.tsv`: genomes not analysed (non-bacterial lineage, e.g. archaea).
  - `unresolved_genomes.tsv`: genomes without a usable lineage, with the reason given by the classifier (e.g. `no_resolved_lineage_from_BIN_CLASSIFICATION (gtdbtk: No bacterial or archaeal marker)`).
:::

## General typing

[MLST](https://github.com/tseemann/mlst), [VirulenceFinder](https://bitbucket.org/genomicepidemiology/virulencefinder) and [PlasmidFinder](https://bitbucket.org/genomicepidemiology/plasmidfinder).

:::{dropdown} Output files
- `Bacterial_Analysis/general_typing/`
  - `mlst/[genome_id].mlst.tsv`: scheme, sequence type and alleles.
  - `virulencefinder/[genome_id].virulencefinder.{tsv,txt,json,standard.json}`, `*-virulence_genes.fsa`, `*-hit_in_genome_seq.fsa`: virulence genes found.
  - `plasmidfinder/[genome_id].plasmidfinder.{tsv,txt,json}`, `*-plasmid_seqs.fsa`, `*-hit_in_genome_seq.fsa`: plasmid replicons found.
  - `summary/[prefix]_{mlst,virulencefinder,plasmidfinder}.tsv`.
:::

## Antimicrobial resistance

[AMRFinderPlus](https://github.com/ncbi/amr), [ABRicate](https://github.com/tseemann/abricate), [RGI](https://github.com/arpcard/rgi) (contig mode), [ResFinder](http://genepi.food.dtu.dk/resfinder) with PointFinder and DisinFinder, and [DeepARG](https://github.com/gaarangoa/deeparg). Each report is converted to one format by [hAMRonization](https://github.com/pha4ge/hAMRonization) and normalised to ARO terms by [ARGnorm](https://github.com/BigDataBiology/argNorm) (RGI is already in ARO terms).

:::{dropdown} Output files
- `Bacterial_Analysis/amr/`
  - `amrfinderplus/[genome_id].amrfinderplus.tsv`
  - `abricate/[genome_id].abricate.tsv`
  - `rgi/[genome_id].rgi.{txt,json}`
  - `resfinder/[genome_id].resfinder.json`, `resfinder/[genome_id].resfinder/` (ResFinder, PointFinder and DisinFinder tables, `pheno_table.txt`, hit sequences)
  - `deeparg/[genome_id].deeparg.mapping.{ARG,potential.ARG}`, `*.align.daa*`
  - `hamronization/[genome_id].[tool]_hamronized.tsv`, `argnorm/[genome_id].[tool]_hamronized_argnorm.tsv`
  - `hamronization_summary/amr_summary.tsv`, `amr_summary_argnorm.tsv`: all genomes and tools in one table.
  - `summary/[prefix]_{amrfinderplus,abricate,rgi,resfinder,pointfinder,disinfinder,deeparg}.tsv`
:::

## *Klebsiella* / *Escherichia*

[Kleborate](https://github.com/klebgenomics/Kleborate) v3 with the preset chosen from the lineage (`--kleborate_preset_map`).

:::{dropdown} Output files
- `Bacterial_Analysis/klebsiella/`
  - `kleborate/[genome_id].kleborate.txt`
  - `summary/[prefix]_kleborate.tsv`
:::

## *Salmonella*

[SISTR](https://github.com/phac-nml/sistr_cmd) and [SeqSero2](https://github.com/denglab/SeqSero2) for every *Salmonella*; the SISTR serovar sends Typhi to [Mykrobe](https://github.com/Mykrobe-tools/mykrobe) + [GenoTyphi](https://github.com/typhoidgenomics/genotyphi) and Paratyphi A to [Paratype](https://github.com/CHRF-Genomics/Paratype).

:::{dropdown} Output files
- `Bacterial_Analysis/salmonella/`
  - `sistr/[genome_id].tab`, `[genome_id]-allele.{fasta,json}`, `[genome_id]-cgmlst.csv`
  - `seqsero2/[genome_id].seqsero2_result.tsv`
  - `genotyphi/[genome_id].mykrobe.json`, `[genome_id].genotype_predictResults.tsv`
  - `paratype/[genome_id].tsv`
  - `routing/[genome_id].salmonella_{routing,evidence}.tsv`: serovar route chosen and the evidence used.
  - `summary/[prefix]_{sistr,seqsero2,genotyphi,paratype,salmonella_typing}.tsv`
:::

## Downloaded databases

Databases downloaded because no path was given are saved for later runs (turn off with `--save_<tool>_db false`):

:::{dropdown} Output files
- `Databases/{AMRFinderPlus,RGI,ResFinder,DeepARG,VirulenceFinder,ParaType}/`
:::
