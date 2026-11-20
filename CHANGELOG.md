# Changelog

## 2.0.0 - 2026-11-20

First release under the name **microcosm** (formerly **metaflow**; the GitHub repository moved from
`tnmquann/metaflow` to `tnmquann/microcosm`, old links redirect). Requires **Nextflow ≥ 26.04** (strict syntax).

### Breaking changes

- **Rename metaflow → microcosm.** The main workflow is `workflow/microcosm.nf` and every process is now named
  `MICROCOSM:...` instead of `METAFLOW:...`: update custom configs that select processes by their full name, and expect
  `-resume` not to reuse tasks of runs made with 1.x.
- **Read-based output file renamed:** `<prefix>.processed_metaflow.csv` → `<prefix>.processed_microcosm.csv`;
  `post_processing.py --input_format raw_metaflow` → `raw_microcosm`.
- **Nextflow 26.04 strict syntax** ([#83](https://github.com/tnmquann/microcosm/issues/83)): the pipeline needs
  Nextflow ≥ 26.04; the custom `check_max` helper is replaced by native `resourceLimits` (set limits with
  `process.resourceLimits` in a custom config).
- `--postprocess_options` is renamed `--readbased_postprocess_options` and now defaults to `--min_coverage 0.05`.

### Added

#### Read-based analysis ([#84](https://github.com/tnmquann/microcosm/issues/84))
- **Post-processing modes** `--readbased_postprocess all|phyloseq|taxburst|rgi_bwt` for batch and single-sketch runs,
  published below `Readbased_Analysis/Post_processing/` (`<batch-prefix>_<mode>/` in batch mode, `<sample_id>/<mode>/`
  plus all-sample aggregates in single-sketch mode). `all` also writes the standard processed CSV; RGI-BWT results are
  included when `--enable_rgi_bwt` is set (otherwise a warning, and the other modes continue).
- `--readbased_postprocess_options`: extra `post_processing.py` options, tokenized safely (module-managed arguments
  cannot be overridden); directory input accepts both `_1/_2` and `_R1/_R2` read pairs.
- **sylph profiler:** `--readbased_taxotools microcosm|sylph|microcosm,sylph` chooses sourmash + YACHT (default), sylph,
  or both. New parameters `--sylph_database` (one or several `*.syldb`/`*.syl2db`), `--sylph_args`, `--enable_sylphtax`
  with `--sylph_taxonomy` (sylph-tax `taxprof`, one metadata file per database), `--sylph_datatype`
  (`taxonomic_abundance` or `sequence_abundance`). Outputs in `Readbased_Analysis/sylph/` (profiles, `.sylphmpa`,
  clade × sample table, phyloseq tables with `--readbased_postprocess phyloseq|all`). RGI-BWT also runs with sylph alone.
- Several sourmash databases (`--sourmash_database db1.zip,db2.zip` with one taxonomy each in
  `--sourmash_taxonomy_csv`); YACHT, which takes one database, is skipped in that case.
- Per-sample YACHT in single-sketch mode; `--enable_copymergedseqs`, `--enable_copysketch`,
  `--enable_copysketch_singles` and `--enable_copyintermediate` control publishing of intermediate files.

#### Bacterial genome downstream analysis ([#86](https://github.com/tnmquann/microcosm/issues/86))
- **Standalone genomes:** `--genome_input` CSV (`sample_id,fasta`, optional `run_id,group,assembler,lineage`). A valid
  GTDB lineage is used as given; a blank lineage is classified by the pipeline. Relative FASTA paths are resolved
  against the CSV location.
- **Bins from the assembly route** go through the same layer with `--enable_bacterial_downstream`.
- Only genomes resolved as `d__Bacteria` are analysed; the others are listed in
  `Bacterial_Analysis/routing/{excluded,unresolved}_genomes.tsv`, with the classifier's reason for unresolved genomes
  (e.g. `no_resolved_lineage_from_BIN_CLASSIFICATION (gtdbtk: No bacterial or archaeal marker)`).
- **Optional GTDB-Tk classifier:** `--bin_classifier gtdbtk --gtdbtk_db <R226>` (`classify_wf`, results in
  `Taxonomy/GTDB-Tk/`; `--gtdbtk_args`, `--gtdbtk_pplacer_cpus`, `--gtdbtk_pplacer_scratch`). sourmash stays the default.
- **General typing:** MLST, VirulenceFinder, PlasmidFinder.
- **AMR:** AMRFinderPlus (`--plus`, organism chosen from the lineage), ABRicate, RGI main (contig mode), ResFinder with
  PointFinder and DisinFinder, DeepARG; every report converted with hAMRonization, summarised, and normalised with
  ARGnorm.
- **Klebsiella / Escherichia:** Kleborate v3 with the preset chosen from the lineage (`kpsc`, `kosc`, `escherichia`).
- **Salmonella:** SISTR and SeqSero2; the SISTR serovar routes `Typhi` to Mykrobe + GenoTyphi and `Paratyphi A` to
  Paratype.
- **Configurable steps:** `--skip_<tool>` for every tool, `--<tool>_args` for every tool, lineage rule tables
  (`--kleborate_preset_map`, `--amrfinderplus_organism_map`; GTDB placeholder suffixes such as `_A` are ignored when
  matching), and an optional per-lineage decision table `--bacterial_decision_table` (y/n step columns plus per-lineage
  `<tool>_args`, `amrfinderplus_organism` and `resfinder_species` columns; example in
  `assets/bacterial_decision_table_example.tsv`).
- **Databases:** every database parameter is optional; a missing database is downloaded at a fixed version (database
  revisions selectable with `--virulencefinder_db_revision`, `--resfinder_db_revision`, `--paratype_revision`) and saved
  to `${outdir}/Databases/<Tool>/` (`--save_amrfinderplus_db`, `--save_rgi_db`, `--save_resfinder_db`,
  `--save_deeparg_db`, `--save_virulencefinder_db`, `--save_paratype_db`, default `true`).
- **Outputs** in `${outdir}/Bacterial_Analysis/` (`routing/`, `amr/`, `general_typing/`, `klebsiella/`, `salmonella/`),
  one file per genome and one `summary/<prefix>_<tool>.tsv` per tool.
- Console task tags name the genome / bin being analysed (e.g. `VIRULENCEFINDER_TYPING (MEGAHIT-COMEBin-sample.0)`).

#### MultiQC report
- New microcosm-styled **MultiQC report** `${outdir}/MultiQC/<prefix>_multiqc_report.html` (`--skip_multiqc`,
  `--multiqc_config`, `--multiqc_extra_config`, `--multiqc_logo`, `--multiqc_replace_names`, `--multiqc_sample_names`).
- Sections: read QC (FastQC raw / after host removal, fastp, hostile); microcosm top taxa (tab per rank, following
  `--readbased_postprocess_options`); sylph profile and sylph-tax; RGI bwt (table + AMR gene families); QUAST and Bowtie2
  for assemblies; contigs (SeqKit + samtools coverage, Pyrodigal or Prokka, geNomad, Tiara, skani); bins (Tiara domains,
  QUAST, assembly_stats, CheckM2, CheckM, BUSCO, GUNC, Prokka / Bakta, GTDB-Tk or sourmash gather, bin taxonomy and
  summary); bacterial genomes (overview: species, MLST, SISTR serovar, AMR genes per tool; Kleborate: ST, scores, K/O
  locus, virulence loci; Salmonella: SISTR/SeqSero2 serovars, cgMLST, follow-up route, GenoTyphi/Paratype genotype; RGI
  main table + gene families); run summary (parameters that differ from the defaults) and software versions.
- Tools without a MultiQC module are summarised by `bin/py_scripts/multiqc/tool_multiqc.py`, one task per tool, so a
  skipped or failing tool never affects the other sections.
- All software versions of a run are collected in `${outdir}/pipeline_info/microcosm_software_versions.yml`.

#### Other
- Bakta and Prokka annotation options exposed as parameters (`--annotation_bakta_*`, `--annotation_prokka_*`,
  `--annotation_pyrodigal_usespecialstopcharacter`, `--prokka_bins_with_compliance`).
- `--skani_search_options`, `--exclude_unbins_from_postbinning`, `--help`, `--help_full`, `--show_hidden`.
- `test_data/`: small reusable inputs (eight RefSeq bacterial genomes with `genome_input` CSVs and decision tables, a
  skani database, a GTDB-Tk failure case, two sourmash databases, a read simulation script), documented in
  `test_data/README.md`.

### Changed

- Strict-syntax refactor of all local modules and subworkflows (nf-core-style layout, `channel.*` factories, explicit
  closure parameters, helper scripts staged as inputs, topic-based version outputs).
- nf-schema plugin 2.7.2; schema updated with all new parameters.
- Every local module has a pinned container for Docker / Singularity / Apptainer and a Conda definition; short
  `biocontainers/` image names resolve on quay.io.
- Read-based post-processing applies `--min_coverage 0.05` by default.
- Copyright: 2026 TON NGOC MINH QUAN.

### Fixed

- Random Bakta crash `AttributeError: 'str' object has no attribute 'decode'` (upstream bakta#424).
- Reads were mapped twice per assembly (BWA-MEM2 for contig analysis, `--binprepare_mappertool` for binning) and two
  BAM sets were published. Reads are now mapped once with `--binprepare_mappertool` (`bowtie2` default, or `bwamem2`);
  contig coverage and binning share the sorted BAM, published once in `Binning/mapping/sorted_bam/`
  (`ContigsAnalysis/mapped_reads/` is gone).
- MultiQC put Prokka results of Binette bins (no `.fa` in the name) in "Prokka: contigs". Contig Prokka files are
  now named `<assembler>-<sample>.contigs.prokka.*`, and every other Prokka file is reported under "Prokka: bins".
- Clean tool versions in `microcosm_software_versions.yml`: Kleborate (no `tput` warning), `sed` in sourmash fastmultigather, `cat` in read merging.
- `ASSEMBLY_STATS` null handling; `BIN_CLASSIFICATION` version export; `DOMAIN_CLASSIFICATION` for MAGs given with
  `--postbinning_input`.
- `HOSTILE_CLEAN` outputs sometimes not recognised; `--spades_options` not applied; `SAMTOOLS_SORT` input errors.
- Duplicate input declarations in the read-based finalization module; `SEQKIT_FX2TAB` `ext.suffix` usage.

### Dependencies

| Tool | 1.0.2 | 2.0.0 |
|---|---|---|
| Nextflow | ≥ 25.04.6 | ≥ 26.04 |
| nf-schema | - | 2.7.2 |
| Bakta | 1.11.4 | 1.12.0 (db download 1.11.4) |
| CheckM | 1.2.3 | 1.2.5 |
| csvtk | 0.31.0 | 0.37.0 |
| fastp | 0.24.0 | 1.1.0 |
| geNomad | 1.11.0 | 1.12.0 |
| GTDB-Tk | 2.5.2 | 2.5.2 (pinned: 2.7.2 needs GTDB R232) |
| GUNC | 1.0.6 | 1.1.0 |
| hostile | 2.0.1 | 2.0.2 |
| MultiQC | 1.29 | 1.35 |
| Pyrodigal | 3.6.3 | 3.7.1 |
| RGI | 6.0.5 | 6.0.8 |
| skani | 0.3.0 | 0.2.2 (current nf-core `skani/search` module) |
| sylph / sylph-tax | - | 1.0.0 / 1.9.2 |
| AMRFinderPlus | - | 4.2.7 |
| ABRicate | - | 1.0.1 |
| ResFinder | - | 4.7.2 |
| DeepARG | - | 1.0.4 |
| hAMRonization | - | 1.1.9 |
| ARGnorm | - | 0.8.0 |
| MLST | - | 2.25.0 |
| VirulenceFinder | - | 3.2.1 |
| PlasmidFinder | - | 2.1.6 |
| Kleborate | - | 3.2.4 |
| SISTR | - | 1.1.1 |
| SeqSero2 | - | 1.2.1 |
| Mykrobe / GenoTyphi | - | 0.11.0 / 1.9.1 |

**Full changelog**: [1.0.2...2.0.0](https://github.com/tnmquann/microcosm/compare/1.0.2...2.0.0)
