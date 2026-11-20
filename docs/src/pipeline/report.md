# MultiQC report

Every run ends with one [MultiQC](https://multiqc.info/) report, `MultiQC/<prefix>_multiqc_report.html`
(`<prefix>` is the last folder name of `--outdir`). A section appears only when its tool ran.

| Part | Sections |
| :--- | :--- |
| **Reads** | FastQC (raw and after host removal), fastp, hostile |
| **Read-based profiles** | microcosm top taxa (tab per rank; needs `--readbased_postprocess all`), sylph profile, sylph-tax, RGI bwt (table + AMR gene families) |
| **Assemblies and contigs** | QUAST, Bowtie2, contig summary (SeqKit + coverage), Pyrodigal or Prokka, geNomad, Tiara, skani |
| **Bins** | Tiara domains, QUAST, assembly_stats, CheckM2, CheckM, BUSCO, GUNC, Prokka / Bakta, GTDB-Tk or sourmash, bin taxonomy and summary |
| **Bacterial genomes** | overview (species, MLST, serovar, AMR genes per tool), Kleborate, *Salmonella* typing, RGI main |
| **Run** | parameters that differ from the defaults, software versions |

| Option | Effect |
| :--- | :--- |
| `--skip_multiqc` | no report |
| `--multiqc_config`, `--multiqc_extra_config` | applied after the default `assets/multiqc_config.yml` |
| `--multiqc_logo` | your logo (light and dark mode) |
| `--multiqc_replace_names`, `--multiqc_sample_names` | rename samples (TSV) |

Every tool version used in the run is listed in `pipeline_info/microcosm_software_versions.yml`.

[MultiQC](https://multiqc.info/) collects the results of the run into one HTML report (`--skip_multiqc` to turn off). Sections appear only for the tools that ran: read QC (FastQC, fastp, hostile), microcosm top taxa, sylph and sylph-tax, RGI bwt, QUAST and Bowtie2 for assemblies, contigs (SeqKit + coverage, Pyrodigal or Prokka, geNomad, Tiara, skani), bins (Tiara, QUAST, assembly_stats, CheckM2, CheckM, BUSCO, GUNC, Prokka / Bakta, GTDB-Tk or sourmash, bin taxonomy and summary), bacterial genomes (overview, Kleborate, Salmonella typing, RGI main), the parameters that differ from the defaults and the software versions.

:::{dropdown} Output files
- `MultiQC/`
  - `[prefix]_multiqc_report.html`: the report.
  - `[prefix]_multiqc_report_data/`: the parsed data behind every section.
- `pipeline_info/`
  - `microcosm_software_versions.yml`: versions of every tool used in the run.
  - Nextflow execution report, timeline and trace.
:::
