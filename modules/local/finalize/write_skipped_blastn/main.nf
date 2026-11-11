process WRITE_SKIPPED_BLASTN {
    tag "${meta.sample_id}:${meta.genome_id}"
    label 'process_single'

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("${meta.genome_id}.blastn.skipped.txt"), emit: txt
    path "versions.yml", emit: versions

    script:
    """
    cat <<-END_SKIP > "${meta.genome_id}.blastn.skipped.txt"
    BLASTN skipped: no --bacterial_blast_reference or --bacterial_blast_db was supplied.
    END_SKIP

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        write_skipped_blastn: 1.0.0
    END_VERSIONS
    """
}
