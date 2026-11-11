process SUMMARIZE_GENERAL_TYPING {
    tag "${meta.sample_id}:${meta.genome_id}"
    label 'process_low'

    input:
    tuple val(meta), path(virulencefinder), path(plasmidfinder), path(mlst), path(blastn_filtered)
    path summary_script

    output:
    tuple val(meta), path("${meta.genome_id}.general_typing_evidence.tsv"), emit: evidence
    tuple val(meta), path("${meta.genome_id}.mlst_summary.tsv"), emit: mlst_summary
    path "versions.yml", emit: versions

    script:
    """
    python ${summary_script} \
        --sample-id "${meta.sample_id}" \
        --genome-id "${meta.genome_id}" \
        --virulencefinder ${virulencefinder} \
        --plasmidfinder ${plasmidfinder} \
        --mlst ${mlst} \
        --blastn-filtered ${blastn_filtered} \
        --output "${meta.genome_id}"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        summarize_gene_hits.py: \$(python ${summary_script} --version)
        python: \$(python --version 2>&1 | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    """
    python ${summary_script} \
        --sample-id "${meta.sample_id}" \
        --genome-id "${meta.genome_id}" \
        --virulencefinder ${virulencefinder} \
        --plasmidfinder ${plasmidfinder} \
        --mlst ${mlst} \
        --blastn-filtered ${blastn_filtered} \
        --output "${meta.genome_id}"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        summarize_gene_hits.py: \$(python ${summary_script} --version)
        python: \$(python --version 2>&1 | sed 's/Python //')
    END_VERSIONS
    """
}
