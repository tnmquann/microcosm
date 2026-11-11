process NORMALIZE_GENOME_CLASSIFICATIONS {
    tag "genome classifications"
    label 'process_low'

    input:
    path classification_csvs
    path normalizer_script

    output:
    path "resolved_genome_lineages.tsv", emit: resolved_tsv
    path "unresolved_classification_rows.tsv", emit: unresolved_tsv
    path "versions.yml", emit: versions

    script:
    """
    python ${normalizer_script} \
        --classification ${classification_csvs} \
        --resolved-output resolved_genome_lineages.tsv \
        --unresolved-output unresolved_classification_rows.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        normalize_genome_classifications.py: \$(python ${normalizer_script} --version)
        python: \$(python --version 2>&1 | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    """
    python ${normalizer_script} \
        --classification ${classification_csvs} \
        --resolved-output resolved_genome_lineages.tsv \
        --unresolved-output unresolved_classification_rows.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        normalize_genome_classifications.py: \$(python ${normalizer_script} --version)
        python: \$(python --version 2>&1 | sed 's/Python //')
    END_VERSIONS
    """
}
