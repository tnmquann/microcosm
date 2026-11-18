process VALIDATE_GENOME_INPUT {
    tag "${genome_csv}"
    label 'process_low'

    conda "conda-forge::python=3.11"
    // Standard-library Python helper.
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.11' :
        'quay.io/biocontainers/python:3.11' }"

    input:
    path genome_csv
    val classification_available
    path validator_script

    output:
    path "validated_genomes.csv", emit: validated_csv
    path "versions.yml", emit: versions

    script:
    def classification_flag = classification_available ? '--classification-available' : ''
    """
    python ${validator_script} \
        --input ${genome_csv} \
        --output validated_genomes.csv \
        ${classification_flag}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        validate_genome_input.py: \$(python ${validator_script} --version)
        python: \$(python --version 2>&1 | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    def classification_flag = classification_available ? '--classification-available' : ''
    """
    python ${validator_script} \
        --input ${genome_csv} \
        --output validated_genomes.csv \
        ${classification_flag}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        validate_genome_input.py: \$(python ${validator_script} --version)
        python: \$(python --version 2>&1 | sed 's/Python //')
    END_VERSIONS
    """
}
