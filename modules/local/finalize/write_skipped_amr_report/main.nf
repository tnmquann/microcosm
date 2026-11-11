process WRITE_SKIPPED_AMR_REPORT {
    tag "${meta.sample_id}:${meta.genome_id}:${tool}"
    label 'process_single'

    input:
    tuple val(meta), path(fasta)
    val tool
    val suffix
    val reason

    output:
    tuple val(meta), path("${meta.genome_id}.${tool}.skipped.${suffix}"), emit: report
    path "versions.yml", emit: versions

    script:
    """
    cat <<-END_SKIP > "${meta.genome_id}.${tool}.skipped.${suffix}"
    ${tool} skipped: ${reason}
    END_SKIP

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        write_skipped_amr_report: 1.0.0
    END_VERSIONS
    """
}
