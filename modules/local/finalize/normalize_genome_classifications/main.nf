process NORMALIZE_GENOME_CLASSIFICATIONS {
    tag "genome classifications"
    label 'process_low'

    conda "conda-forge::python=3.11"
    // Standard-library Python helper.
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.11' :
        'quay.io/biocontainers/python:3.11' }"

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
