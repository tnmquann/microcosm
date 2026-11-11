process VALIDATE_BACTERIAL_AMR_DATABASES {
    tag "bacterial AMR databases"
    label 'process_low'

    input:
    val validate_trigger
    val amrfinderplus_db
    val abricate_db
    val abricate_db_id
    val rgi_db
    val include_rgi_wildcard
    val resfinder_db
    val pointfinder_db
    val deeparg_db
    val deeparg_db_version
    path validator_script

    output:
    path "amr_database_validation.tsv", emit: validation
    path "amrfinder_supported_organisms.tsv", emit: amrfinder_organisms
    path "pointfinder_supported_species.tsv", emit: pointfinder_species
    path "amr_databases.ready", emit: ready
    path "versions.yml", emit: versions

    script:
    def wildcard_flag = include_rgi_wildcard ? '--include-rgi-wildcard' : ''
    """
    python ${validator_script} \
        --amrfinderplus-db "${amrfinderplus_db}" \
        --abricate-db "${abricate_db ?: ''}" \
        --abricate-db-id "${abricate_db_id}" \
        --rgi-db "${rgi_db}" \
        ${wildcard_flag} \
        --resfinder-db "${resfinder_db}" \
        --pointfinder-db "${pointfinder_db}" \
        --deeparg-db "${deeparg_db}" \
        --deeparg-db-version "${deeparg_db_version}" \
        --validation-output amr_database_validation.tsv \
        --amrfinder-organisms-output amrfinder_supported_organisms.tsv \
        --pointfinder-species-output pointfinder_supported_species.tsv

    touch amr_databases.ready

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        validate_amr_databases.py: \$(python ${validator_script} --version)
        python: \$(python --version 2>&1 | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    def wildcard_flag = include_rgi_wildcard ? '--include-rgi-wildcard' : ''
    """
    python ${validator_script} \
        --amrfinderplus-db "${amrfinderplus_db}" \
        --abricate-db "${abricate_db ?: ''}" \
        --abricate-db-id "${abricate_db_id}" \
        --rgi-db "${rgi_db}" \
        ${wildcard_flag} \
        --resfinder-db "${resfinder_db}" \
        --pointfinder-db "${pointfinder_db}" \
        --deeparg-db "${deeparg_db}" \
        --deeparg-db-version "${deeparg_db_version}" \
        --validation-output amr_database_validation.tsv \
        --amrfinder-organisms-output amrfinder_supported_organisms.tsv \
        --pointfinder-species-output pointfinder_supported_species.tsv

    touch amr_databases.ready

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        validate_amr_databases.py: \$(python ${validator_script} --version)
        python: \$(python --version 2>&1 | sed 's/Python //')
    END_VERSIONS
    """
}
