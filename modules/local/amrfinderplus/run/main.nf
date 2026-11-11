/*
 * Local adaptation of nf-core/modules amrfinderplus/run.
 *
 * The installed upstream AMRFinderPlus 4.2.7 process remains unchanged. This
 * process pins AMRFinderPlus 3.12.8 to support the explicitly validated 3.12.x
 * database format and uses deterministic local version/stub outputs.
 */
process AMRFINDERPLUS_RUN_LOCAL {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/ncbi-amrfinderplus:3.12.8--h283d18e_0'
        : 'quay.io/biocontainers/ncbi-amrfinderplus:3.12.8--h283d18e_0'}"

    input:
    tuple val(meta), path(fasta)
    path db

    output:
    tuple val(meta), path("${prefix}.tsv"), emit: report
    tuple val(meta), path("${prefix}-mutations.tsv"), emit: mutation_report, optional: true
    env 'VER', emit: tool_version
    env 'DBVER', emit: db_version
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def is_compressed_fasta = fasta.getName().endsWith(".gz")
    def is_compressed_db = db.getName().endsWith(".gz")
    def organism_param = meta.containsKey("organism")
        ? "--organism ${meta.organism} --mutation_all ${prefix}-mutations.tsv"
        : ""
    def fasta_name = fasta.getName().replace(".gz", "")
    def fasta_param = meta.get("is_proteins", false) ? "-p" : "-n"

    """
    if [ "${is_compressed_fasta}" == "true" ]; then
        gzip -c -d ${fasta} > ${fasta_name}
    fi

    if [ "${is_compressed_db}" == "true" ]; then
        mkdir amrfinderdb
        tar xzvf ${db} -C amrfinderdb
    else
        mv ${db} amrfinderdb
    fi

    amrfinder \\
        ${fasta_param} ${fasta_name} \\
        ${organism_param} \\
        ${args} \\
        --database amrfinderdb \\
        --threads ${task.cpus} \\
        > ${prefix}.tsv

    VER=\$(amrfinder --version)
    DBVER=\$(amrfinder --database amrfinderdb --database_version 2>&1 | grep -oE '[0-9]{4}-[0-9]{2}-[0-9]{2}\\.[0-9]+' | tail -1)

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        amrfinderplus: \$(echo \$VER)
        amrfinderplus-database: \$(echo \$DBVER)
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"

    """
    echo ${args}

    touch ${prefix}.tsv
    ${meta.containsKey("organism") ? "touch ${prefix}-mutations.tsv" : ""}

    VER=3.12.8
    DBVER=2024-07-22.1

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        amrfinderplus: \$(echo \$VER)
        amrfinderplus-database: \$(echo \$DBVER)
    END_VERSIONS
    """
}
