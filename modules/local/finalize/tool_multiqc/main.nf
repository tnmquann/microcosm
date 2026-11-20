process TOOL_MULTIQC {
    tag "${tool}"
    label 'process_single'

    conda "conda-forge::python=3.11"
    // Standard-library Python helper.
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.11' :
        'quay.io/biocontainers/python:3.11' }"

    input:
    // One task per tool: a skipped or failing tool never touches the other MultiQC sections.
    // in*/: tools such as GUNC write the same file name for every bin set.
    tuple val(tool), path(inputs, stageAs: 'in*/*')
    path multiqc_script

    output:
    path "*_mqc.json", emit: tables, optional: true
    tuple val("${task.process}"), val('python'), eval("python --version 2>&1 | sed 's/Python //'"), emit: versions_python, topic: versions

    script:
    """
    python ${multiqc_script} --tool ${tool} ${task.ext.args ?: ''} --outdir . in*
    """

    stub:
    """
    touch ${tool}_mqc.json
    """
}
