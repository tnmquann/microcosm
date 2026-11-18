process MICROCOSM_MULTIQC_TABLES {
    tag "${processed_csv.name}"
    label 'process_single'

    conda "conda-forge::python=3.11"
    // Standard-library Python helper.
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.11' :
        'quay.io/biocontainers/python:3.11' }"

    input:
    path processed_csv
    path multiqc_script

    output:
    path "microcosm_*_mqc.tsv", emit: tables
    tuple val("${task.process}"), val('python'), eval("python --version 2>&1 | sed 's/Python //'"), emit: versions_python, topic: versions

    script:
    def args = task.ext.args ?: ''
    """
    python ${multiqc_script} --input ${processed_csv} --outdir . ${args}
    """

    stub:
    """
    touch microcosm_species_mqc.tsv
    """
}
