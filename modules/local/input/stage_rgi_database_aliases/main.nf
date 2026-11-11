process STAGE_RGI_DATABASE_ALIASES {
    tag "RGI CARD/WildCARD aliases"
    label 'process_low'

    input:
    val rgi_database
    path validation_ready

    output:
    path "card_db", emit: card
    path "wildcard_db", emit: wildcard
    path "versions.yml", emit: versions

    script:
    """
    ln -s "${rgi_database}" card_db
    ln -s "${rgi_database}" wildcard_db

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        coreutils: \$(ln --version | head -n 1 | awk '{print \$NF}')
    END_VERSIONS
    """

    stub:
    """
    ln -s "${rgi_database}" card_db
    ln -s "${rgi_database}" wildcard_db

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        coreutils: "stub_version"
    END_VERSIONS
    """
}
