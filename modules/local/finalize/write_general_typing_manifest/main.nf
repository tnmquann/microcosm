process WRITE_GENERAL_TYPING_MANIFEST {
    tag "${meta.sample_id}:${meta.genome_id}"
    label 'process_low'

    input:
    tuple val(meta), path(virulencefinder_tsv), path(plasmidfinder_tsv), path(mlst_tsv), path(blastn_output), path(blastn_filtered_output)
    val blast_reference_type
    val blast_reference_label
    path manifest_script

    output:
    tuple val(meta), path("${meta.genome_id}.general_typing_manifest.tsv"), emit: manifest
    path "versions.yml", emit: versions

    script:
    def blast_evidence = blast_reference_type == 'not_configured' ? blastn_output : blastn_filtered_output
    """
    python ${manifest_script} \
        --sample-id "${meta.sample_id}" \
        --genome-id "${meta.genome_id}" \
        --source-route "${meta.source_route}" \
        --virulencefinder-tsv ${virulencefinder_tsv} \
        --plasmidfinder-tsv ${plasmidfinder_tsv} \
        --mlst-tsv ${mlst_tsv} \
        --blastn-output ${blastn_output} \
        --blastn-filtered-output ${blast_evidence} \
        --blast-reference-type "${blast_reference_type}" \
        --blast-reference-label "${blast_reference_label}" \
        --output "${meta.genome_id}.general_typing_manifest.tsv"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        write_general_typing_manifest.py: \$(python ${manifest_script} --version)
        python: \$(python --version 2>&1 | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    def blast_evidence = blast_reference_type == 'not_configured' ? blastn_output : blastn_filtered_output
    """
    python ${manifest_script} \
        --sample-id "${meta.sample_id}" \
        --genome-id "${meta.genome_id}" \
        --source-route "${meta.source_route}" \
        --virulencefinder-tsv ${virulencefinder_tsv} \
        --plasmidfinder-tsv ${plasmidfinder_tsv} \
        --mlst-tsv ${mlst_tsv} \
        --blastn-output ${blastn_output} \
        --blastn-filtered-output ${blast_evidence} \
        --blast-reference-type "${blast_reference_type}" \
        --blast-reference-label "${blast_reference_label}" \
        --output "${meta.genome_id}.general_typing_manifest.tsv"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        write_general_typing_manifest.py: \$(python ${manifest_script} --version)
        python: \$(python --version 2>&1 | sed 's/Python //')
    END_VERSIONS
    """
}
