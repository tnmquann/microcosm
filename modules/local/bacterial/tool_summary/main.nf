process BACTERIAL_TOOL_SUMMARY {
    tag "${module}:${tool}"
    label 'process_single'

    conda "conda-forge::python=3.11"
    // Standard-library Python helper.
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.11' :
        'quay.io/biocontainers/python:3.11' }"

    input:
    tuple val(module), val(tool), val(genome_ids), path(tables, stageAs: 'tables/??/*')
    val prefix
    path summary_script

    output:
    tuple val(module), val(tool), path("${prefix}_${tool}.tsv"), emit: summary
    tuple val("${task.process}"), val('python'), eval("python --version 2>&1 | sed 's/Python //'"), emit: versions_python, topic: versions

    script:
    // mlst writes no header: FILE SCHEME ST allele...
    def header_args = tool == 'mlst' ? '--no-header --columns file,scheme,st,alleles' : ''
    def manifest_args = [genome_ids, [tables].flatten()].transpose()
        .collect { genome_id, table -> "'${genome_id}' '${table}'" }
        .join(' ')
    """
    printf '%s\\t%s\\n' ${manifest_args} > manifest.tsv

    python ${summary_script} \\
        --manifest manifest.tsv \\
        --output ${prefix}_${tool}.tsv \\
        ${header_args}
    """

    stub:
    """
    printf 'genome_id\\n' > ${prefix}_${tool}.tsv
    """
}
