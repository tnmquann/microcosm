process FILTER_BACTERIAL_BLAST_HITS {
    tag "${meta.sample_id}:${meta.genome_id}"
    label 'process_low'

    input:
    tuple val(meta), path(blast_results)
    val minimum_identity
    val minimum_coverage
    val minimum_overlap_bp
    path filter_script

    output:
    tuple val(meta), path("${meta.genome_id}.blastn_filtered_hits.tsv"), emit: filtered_hits
    tuple val(meta), path("${meta.genome_id}.blastn_summary.tsv"), emit: summary
    path "versions.yml", emit: versions

    script:
    """
    python ${filter_script} \
        --blast-results ${blast_results} \
        --min-identity ${minimum_identity} \
        --min-coverage ${minimum_coverage} \
        --overlap-bp ${minimum_overlap_bp} \
        --output "${meta.genome_id}"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        filter_blasthits.py: \$(python ${filter_script} --version)
        python: \$(python --version 2>&1 | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    """
    python ${filter_script} \
        --blast-results ${blast_results} \
        --min-identity ${minimum_identity} \
        --min-coverage ${minimum_coverage} \
        --overlap-bp ${minimum_overlap_bp} \
        --output "${meta.genome_id}"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        filter_blasthits.py: \$(python ${filter_script} --version)
        python: \$(python --version 2>&1 | sed 's/Python //')
    END_VERSIONS
    """
}
