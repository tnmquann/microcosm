# Troubleshooting & FAQ

## Upgrading from metaflow 1.x

- Use Nextflow ≥ 26.04.
- Custom configs selecting `METAFLOW:...` processes must use `MICROCOSM:...`.
- `<prefix>.processed_metaflow.csv` is now `<prefix>.processed_microcosm.csv`.
- `--postprocess_options` is now `--readbased_postprocess_options`.
- Resource limits: the old `*_max_cpu`, `*_max_mem`, `*_max_time` parameters are gone; use `--max_cpus`, `--max_memory`, `--max_time` or `process.resourceLimits`.
- Bacterial genome results are in `${outdir}/Bacterial_Analysis/`.
- Reads are mapped to the assemblies once (`--binprepare_mappertool`, Bowtie2 by default): `ContigsAnalysis/mapped_reads/` is gone; the BAM files are in `Binning/mapping/sorted_bam/`.


## Why is the proportion of unclassified reads so high?

:::{note}
**Understanding Unclassified Reads**
A high unclassified proportion (>50%) represents the estimated DNA mass in your sample that `sourmash gather` could not match to any reference genome.

Accurate relative abundance estimation relies on the reference database having near-perfect coverage of the strains in your sample. In metagenomics, this is often theoretically ideal but practically difficult to achieve.
:::

**Common causes:**

1.  **Unweighted results:**
    *   *Issue:* Calculating k-mer proportions based only on presence/absence.
    *   *Impact:* Sequencing errors generate unique k-mers that don't match references, inflating the unclassified count.

2.  **Weighted results issues:**
    *   **Strain-level bias:** If your sample contains a specific strain absent from the database (or only a distant relative is present), its DNA mass cannot be accurately assigned.
    *   **Low ANI (Average Nucleotide Identity):** Low genetic similarity between sample and reference reduces k-mer matches, underestimating the identified DNA mass (`bp_match`).
    *   **High scaled parameter ($\sigma$):** A large scaled value (e.g., 10,000) undersamples the data. Rare organisms may be missed entirely, or their abundance estimated unreliably.
    *   **Inappropriate k-mer size:** A very large $k$ (e.g., 51) makes the search highly specific but "brittle." Minor evolutionary divergences prevent k-mers from matching, reducing sensitivity.

3.  **Assembler limitations:**
    *   Low-abundance organisms may not be fully recovered during assembly.
    *   Assemblers often collapse strain variants into a single consensus, losing specific variation data.

## Error: Incompatible sourmash/YACHT database versions

If the pipeline fails due to database incompatibility, verify that your `sourmash` and `YACHT` files share the exact same attributes:
*   **Database Type:** (e.g., GTDB, RefSeq)
*   **k-mer Size:** (e.g., 31, 51)
*   **Scaled Value:** (e.g., 1000, 2000)

*Solution:* You must regenerate the databases with matching parameters or download compatible pre-computed versions.

## How is Relative Abundance calculated?

In the read-based sub-workflow, relative abundance is computed in `${projectDir}/bin/py_scripts/read_based/process_sourmash_results.py` using the formula:

$$ \boxed{\text{relative\\_abundance} = \left(\frac{\text{unique\\_intersect\\_bp}}{\text{scaled}} \times \text{median\\_abund}\right) \times \frac{1}{\text{total\\_weighted\\_hashes}}} $$

:::{dropdown} <strong>Why use <code>median_abund</code> instead of <code>mean_abund</code>?</strong>
We use `median_abund` because it is robust against outliers. In complex metagenomes, coverage depth can vary wildly across a genome due to sequencing bias or conserved regions. The median provides a more stable representation of the "typical" abundance of a species than the mean.
:::

## `Max. contig ID length` error when annotating MAGs with Prokka
[This issue](https://github.com/tseemann/prokka/issues/76) typically shows up when Prokka is given assemblies with very long contig identifiers, which is common for MAGs produced by metaSPAdes.
* Option 1 (recommended): Add `--annotation_prokka_retaincontigheaders` so Prokka keeps the original contig names instead of trying to rewrite them. This is the safest fix and works well in most pipelines.
* Option 2 (generally not preferred): Switch the assembler from metaSPAdes to MEGAHIT. metaSPAdes tends to generate extremely long contig IDs, which is what usually triggers this error.
* Option 3: Use Bakta for annotation (`--bin_annotation_tool bakta`), or skip the bin annotation step entirely with `--skip_bin_annotation` and run annotation manually later if you need more con

## Reporting other issues

If you encounter an error not covered here, please reach out via the GitHub repository issues page or contact us at `quantnm@oucru.org` or `minhquan.tdn.ct1619@gmail.com`. To help us debug effectively, please include:

*   **Versions:** microcosm and Nextflow versions (`pipeline_info/microcosm_software_versions.yml` lists all tool versions).
*   **Environment:** Operating System (OS).
*   **Error Log:** The complete error message from the terminal and any relevant `.log` files.
*   **Context:** The exact command run (including parameters) and steps taken prior to failure.
*   **Data:** A description of the input data (or a sample file, if shareable).

We are committed to resolving your issues promptly!
