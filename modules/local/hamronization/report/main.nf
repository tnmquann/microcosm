/*
 * Local consolidation of the nf-core/modules hAMRonization report wrappers
 * for AMRFinderPlus, ABRicate, RGI, DeepARG and ResFinder.
 *
 * The installed upstream modules remain unchanged. This adapter keeps their
 * command semantics while using a Python version compatible with
 * hAMRonization 1.1.9 and a deterministic versions.yml/stub contract.
 */
process HAMRONIZATION_REPORT_LOCAL {
    tag "${meta.id}:${tool}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/hamronization:1.1.9--pyhdfd78af_0'
        : 'quay.io/biocontainers/hamronization:1.1.9--pyhdfd78af_0'}"

    input:
    tuple val(meta), path(report)
    val tool
    val format
    val software_version
    val reference_db_version
    val include_input_file_name

    output:
    tuple val(meta), path("*.json"), optional: true, emit: json
    tuple val(meta), path("*.tsv"), optional: true, emit: tsv
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def input_file_arg = include_input_file_name
        ? "--input_file_name ${meta.genome_id ?: prefix}"
        : ''
    // ResFinder JSON carries its own software/database versions; hamronize
    // resfinder rejects the version options (playground/modular/6_hamr.sh).
    def version_args = software_version
        ? "--analysis_software_version ${software_version} --reference_database_version ${reference_db_version}"
        : ''
    if (!(tool in ['amrfinderplus', 'abricate', 'rgi', 'deeparg', 'resfinder'])) {
        error("Unsupported hAMRonization report tool: ${tool}")
    }
    if (!(format in ['tsv', 'json'])) {
        error("Unsupported hAMRonization report format: ${format}")
    }

    """
    hamronize \\
        ${tool} \\
        ${report} \\
        ${args} \\
        --format ${format} \\
        ${version_args} \\
        ${input_file_arg} \\
        --output ${prefix}.${format}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hamronization: \$(hamronize --version 2>&1 | sed 's/hamronize //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    if (!(tool in ['amrfinderplus', 'abricate', 'rgi', 'deeparg', 'resfinder'])) {
        error("Unsupported hAMRonization report tool: ${tool}")
    }
    if (!(format in ['tsv', 'json'])) {
        error("Unsupported hAMRonization report format: ${format}")
    }

    """
    touch ${prefix}.${format}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hamronization: "1.1.9"
    END_VERSIONS
    """
}
