process BACTERIAL_MULTIQC_TABLE {
    tag "bacterial_overview"
    label 'process_single'

    conda "conda-forge::python=3.11"
    // Standard-library Python helper.
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.11' :
        'quay.io/biocontainers/python:3.11' }"

    input:
    path routing_manifest
    path mlst_summary, stageAs: 'mlst/*'
    path salmonella_summary, stageAs: 'salmonella/*'
    path amr_summary, stageAs: 'amr/*'
    path table_script

    output:
    path "bacterial_overview_mqc.tsv", emit: table
    tuple val("${task.process}"), val('python'), eval("python --version 2>&1 | sed 's/Python //'"), emit: versions_python, topic: versions

    script:
    """
    python ${table_script} \\
        --routing ${routing_manifest} \\
        ${mlst_summary ? "--mlst ${mlst_summary}" : ''} \\
        ${salmonella_summary ? "--salmonella ${salmonella_summary}" : ''} \\
        ${amr_summary ? "--amr ${amr_summary}" : ''} \\
        --output bacterial_overview_mqc.tsv
    """

    stub:
    """
    touch bacterial_overview_mqc.tsv
    """
}
