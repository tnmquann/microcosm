/*
 * Local adaptation of nf-core/modules argnorm.
 *
 * The installed upstream process is retained unchanged under modules/nf-core.
 * This process uses the same ARGnorm 0.8.0 command contract but emits a regular
 * versions.yml and a suffix-correct stub output, allowing deterministic local
 * stub validation without executing `argnorm --version` outside its container.
 */
process ARGNORM_LOCAL {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/argnorm:0.8.0--pyhdfd78af_0'
        : 'quay.io/biocontainers/argnorm:0.8.0--pyhdfd78af_0'}"

    input:
    tuple val(meta), path(input_tsv)
    val tool
    val db

    output:
    tuple val(meta), path("*.tsv"), emit: tsv
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}.argnorm.tsv"
    def db_args = db ? "--db ${db}" : ''
    if (!prefix.endsWith('.tsv')) {
        error("ARGNORM_LOCAL output prefix must end with .tsv: ${prefix}")
    }
    if (!tool) {
        error('ARGNORM_LOCAL tool was not provided.')
    }
    if (tool == 'abricate' && !db) {
        error('ARGNORM_LOCAL requires an ABRicate database identifier.')
    }

    """
    argnorm \
        ${tool} \
        -i ${input_tsv} \
        -o ${prefix} \
        ${db_args} \
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        argnorm: \$(argnorm --version)
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}.argnorm.tsv"
    if (!prefix.endsWith('.tsv')) {
        error("ARGNORM_LOCAL output prefix must end with .tsv: ${prefix}")
    }
    if (!tool) {
        error('ARGNORM_LOCAL tool was not provided.')
    }
    if (tool == 'abricate' && !db) {
        error('ARGNORM_LOCAL requires an ABRicate database identifier.')
    }

    """
    touch ${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        argnorm: "stub_version"
    END_VERSIONS
    """
}
