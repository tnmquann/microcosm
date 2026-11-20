# Running the pipeline

## Running the pipeline

For read-based sub-workflow, the typical command for running the pipeline is as follows:

```bash
nextflow run tnmquann/microcosm -r 2.0.0 \
    -profile conda \
    --input /path/to/your/samples.csv \
    --input_format csv \
    --outdir /path/to/output/directory \
    --enable_readbase \
    --sourmash_database /path/to/your/sourmash_database \
    --sourmash_taxonomy_csv /path/to/your/lineages.csv \
    --yacht_database /path/to/your/yacht_database.json \
    --sourmash_ksize 31
```

For assembly-based sub-workflow (the default when `--enable_readbase` is not given), the typical command is:

```bash
nextflow run tnmquann/microcosm -r 2.0.0 \
    -profile conda \
    --input /path/to/your/samples.csv \
    --input_format csv \
    --outdir /path/to/output/directory \
    --sourmash_database /path/to/your/sourmash_database \
    --sourmash_taxonomy_csv /path/to/your/lineages.csv
```

For your own genome FASTA files (bacterial genome analysis only), use `--genome_input` instead of `--input` (see [Genome Input](inputs.md#genome-input)).

This will launch the pipeline with the `conda` configuration profile. See below for more information about profiles.

Note that the pipeline will create the following files in your working directory:

```
work # Directory containing the nextflow working files
<OUTDIR> # Finished results in specified location (defined with --outdir)
.nextflow_log # Log files from Nextflow
# Other nextflow hidden files, eg. history of pipeline runs and old logs.
```

If you wish to repeatedly use the same parameters for multiple runs, rather than specifying each flag in the command, you can specify these in a params file.

Pipeline settings can be provided in a `yaml` or `json` file via `-params-file <file>`.

Do not use `-c <file>` to specify parameters as this will result in errors. Custom config files specified with `-c` must only be used for [tuning process resource specifications](https://nf-co.re/docs/usage/configuration#tuning-workflow-resources), other infrastructural tweaks (such as output directories), or module arguments (args).

The above pipeline run specified with a params file in yaml format:

```bash
nextflow run main.nf -profile conda -params-file params.yaml
```

with:

```yaml
input: './samplesheet.csv'
input_format: 'csv'
outdir: './results/'
<...>
```

You can also generate such `YAML`/`JSON` files via [nf-core/launch](https://nf-co.re/launch).


## Updating the pipeline

When you run the above command, Nextflow automatically pulls the pipeline code from GitHub and stores it as a cached version. When running the pipeline after this, it will always use the cached version if available - even if the pipeline has been updated since. To make sure that you’re running the latest version of the pipeline, make sure that you regularly update the cached version of the pipeline:
```bash
nextflow pull tnmquann/microcosm
```


## Version Control

Specify the pipeline version with `-r` to ensure reproducibility:
```bash
nextflow run tnmquann/microcosm -r 2.0.0
```
Find available versions on the [releases page](https://github.com/tnmquann/microcosm/releases).


## Reproducibility

- Use parameter files for consistent runs.
- Do not include cluster-specific paths or institutional profiles if sharing.
- For reproducible assembly results, use `--spades_fix_cpus` or `--spadeshybrid_fix_cpus` for SPAdes, and `--megahit_fix_cpu_1` for MEGAHIT.
- MetaBAT2 uses a fixed seed by default.
- For BUSCO, download required lineage datasets in advance and specify with `--busco_db` for reproducible results.


## Nextflow arguments

### `-profile`
Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments.

Several generic profiles are bundled with the pipeline which instruct the pipeline to use software packaged using different methods (Docker, Singularity, Podman, Shifter, Charliecloud, Apptainer, Conda) - see below.

The pipeline also dynamically loads configurations from [https://github.com/nf-core/configs](https://github.com/nf-core/configs) when it runs, making multiple config profiles for various institutional clusters available at run time. For more information and to check if your system is supported, please see the [nf-core/configs documentation](https://github.com/nf-core/configs#documentation).

Note that multiple profiles can be loaded, for example: `-profile test,docker` - the order of arguments is important! They are loaded in sequence, so later profiles can overwrite earlier profiles.

If `-profile` is not specified, the pipeline will run locally and expect all software to be installed and available on the `PATH`. This is *not* recommended, since it can lead to different results on different machines dependent on the computer environment.

- `test`
  - A profile with a complete configuration for automated testing
  - Includes links to test data so needs no other parameters
- `docker`
  - A generic configuration profile to be used with [Docker](https://docker.com/)
- `singularity`
  - A generic configuration profile to be used with [Singularity](https://sylabs.io/docs/)
- `apptainer`
  - A generic configuration profile to be used with [Apptainer](https://apptainer.org/)
- `wave`
  - A generic configuration profile to enable [Wave](https://seqera.io/wave/) containers. Use together with one of the above (requires Nextflow `24.03.0-edge` or later).
- `conda`
  - A generic configuration profile to be used with [Conda](https://conda.io/docs/). Please only use Conda as a last resort i.e. when it’s not possible to run the pipeline with Docker, Singularity, Podman, Shifter, Charliecloud, or Apptainer.

### `-resume`
Specify this when restarting a pipeline. Nextflow will use cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously. For input to be considered the same, not only the names must be identical but the files’ contents as well. For more info about this parameter, see [this blog post](https://www.nextflow.io/blog/2019/demystifying-nextflow-resume.html).

You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.

### `-c`
Specify the path to a specific config file (this is a core Nextflow command). See the [nf-core website documentation](https://nf-co.re/usage/configuration) for more information.


## Resource Requests

Default resource settings per process label:

| Label               | CPU | Memory | Time     |
|---------------------|-----|--------|----------|
| Default             | 1   | 16 GB  | 8 hours  |
| process_single      | 1   | 6 GB   | 4 hours  |
| process_low         | 4   | 12 GB  | 4 hours  |
| process_medium      | 12  | 60 GB  | 16 hours |
| process_high        | 25  | 120 GB | 24 hours |
| process_long        | —   | —      | 72 hours |
| process_high_memory | —   | 256 GB | —        |

For most of the pipeline steps, if the job exits with any of the error codes specified [here](https://github.com/nf-core/rnaseq/blob/4c27ef5610c87db00c3c5a3eed10b1d161abf575/conf/base.config#L18) it will automatically be resubmitted with higher resources request (2 x original, then 3 x original). If it still fails after the third attempt then the pipeline execution is stopped.

All requests are capped by `--max_cpus`, `--max_memory` and `--max_time` (Nextflow `process.resourceLimits`). To change the resource requests, please see the [max resources](https://nf-co.re/docs/usage/configuration#max-resources) and [tuning workflow resources](https://nf-co.re/docs/usage/configuration#tuning-workflow-resources) section of the nf-core website.


## Custom containers & tool arguments

Override containers or conda environments for specific tools using the `process` scope and `withName` selector in a custom config file.

**Example:**
```groovy
process {
    withName: TOOL_NAME {
        container = 'quay.io/biocontainers/tool:version'
    }
}
```
Pass the config file with `-c /path/to/custom.config`.
Please see the [updating tool versions](https://nf-co.re/docs/usage/configuration#updating-tool-versions) section of the nf-core website for more information.


## nf-core/configs

In most cases, you will only need to create a custom config as a one-off but if you and others within your organisation are likely to be running nf-core pipelines regularly and need to use the same settings regularly it may be a good idea to request that your custom config file is uploaded to the `nf-core/configs` git repository. Before you do this please can you test that the config file works with your pipeline of choice using the `-c` parameter. For institutional or shared configs, request addition to the [nf-core/configs](https://github.com/nf-core/configs) repository.

See the main [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html) for more information about creating your own configuration files.


## Background execution

Nextflow handles job submissions and supervises the running jobs. The Nextflow process must run until the pipeline is finished. You can run Nextflow in the background with `-bg`, or use `screen`/`tmux` for detached sessions.


## Memory requirements

In some cases, the Nextflow Java virtual machines can start to request a large amount of memory. We recommend adding the following line to your environment to limit this (typically in `~/.bashrc` or `~/.bash_profile`):
```bash
NXF_OPTS='-Xms1g -Xmx4g'
```


## Bin refinement notes

### DASTool

DAS Tool may not always be able to refine bins due to insufficient recovery of enough single-copy genes. In these cases you will get a NOTE such as

```bash
[16/d330a6] NOTE: Process MICROCOSM:BINNING_REFINEMENT:DASTOOL_DASTOOL (test_minigut_sample2)
```

In this case, bin refinement has not necessarily failed but was unable to complete the refinement. You will therefore not expect to find any output files in the `Binning/DASTool/` results directory for that particular sample.

If you are regularly getting such errors, you can try reducing the `--refine_bins_dastool_threshold` value, which will modify the scoring threshold defined in the [DAS Tool publication](https://www.nature.com/articles/s41564-018-0171-1).

### Binette

Binette requires CheckM2 database to proceed. When using this tool, you need to provide the additional parameter `--checkm2_db`.

Binette may not always be able to refine bins due to insufficient recovery of enough single-copy genes. In these cases you will get a NOTE such as

```bash
[16/d330a6] NOTE: Process MICROCOSM:BINNING_REFINEMENT:BINETTE_BINETTE (test_minigut_sample2)
```

In this case, bin refinement has not necessarily failed but was unable to complete the refinement. You will therefore not expect to find any output files in the `Binning/Binette/` results directory for that particular sample.

If you are regularly getting such errors, you can try reducing the `--refine_bins_binette_min_completeness` or`--refine_bins_binette_contamination_weight`value.


The output folders are described in [Outputs](../reference/outputs.md).
