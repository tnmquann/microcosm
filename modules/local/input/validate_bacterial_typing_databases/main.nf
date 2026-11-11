process VALIDATE_BACTERIAL_TYPING_DATABASES {
    tag "bacterial general-typing databases"
    label 'process_low'

    input:
    val validate_trigger
    val virulencefinder_db
    val blast_reference
    val blast_db
    tuple val(blast_db_label), val(blast_min_identity), val(blast_min_coverage), path(validator_script)

    output:
    path "general_typing_database_validation.tsv", emit: validation
    path "general_typing_databases.ready", emit: ready
    path "versions.yml", emit: versions

    script:
    def reference_arg = blast_reference ? "--blast-reference \"${blast_reference}\"" : ''
    def database_arg = blast_db ? "--blast-db \"${blast_db}\"" : ''
    def identity_arg = blast_min_identity.toString() ? "--blast-min-identity ${blast_min_identity}" : ''
    def coverage_arg = blast_min_coverage.toString() ? "--blast-min-coverage ${blast_min_coverage}" : ''
    """
    python ${validator_script} \
        --virulencefinder-db "${virulencefinder_db}" \
        ${reference_arg} \
        ${database_arg} \
        --blast-db-label "${blast_db_label}" \
        ${identity_arg} \
        ${coverage_arg} \
        --validation-output general_typing_database_validation.tsv

    touch general_typing_databases.ready

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        validate_general_typing_databases.py: \$(python ${validator_script} --version)
        python: \$(python --version 2>&1 | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    def reference_arg = blast_reference ? "--blast-reference \"${blast_reference}\"" : ''
    def database_arg = blast_db ? "--blast-db \"${blast_db}\"" : ''
    def identity_arg = blast_min_identity.toString() ? "--blast-min-identity ${blast_min_identity}" : ''
    def coverage_arg = blast_min_coverage.toString() ? "--blast-min-coverage ${blast_min_coverage}" : ''
    """
    python ${validator_script} \
        --virulencefinder-db "${virulencefinder_db}" \
        ${reference_arg} \
        ${database_arg} \
        --blast-db-label "${blast_db_label}" \
        ${identity_arg} \
        ${coverage_arg} \
        --validation-output general_typing_database_validation.tsv

    touch general_typing_databases.ready

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        validate_general_typing_databases.py: \$(python ${validator_script} --version)
        python: \$(python --version 2>&1 | sed 's/Python //')
    END_VERSIONS
    """
}
