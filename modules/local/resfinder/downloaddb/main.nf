process RESFINDER_DOWNLOADDB_LOCAL {
    tag "CGE databases ${revision}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/resfinder:4.7.2--pyhdfd78af_0' :
        'quay.io/biocontainers/resfinder:4.7.2--pyhdfd78af_0' }"

    input:
    val revision

    output:
    path 'resfinder_db', emit: resfinder_db
    path 'pointfinder_db', emit: pointfinder_db
    path 'disinfinder_db', emit: disinfinder_db
    path 'versions.yml', emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // The three repositories cloned by 5_amr.sh, pinned to one CGE release tag.
    // No KMA index: assemblies are searched with BLAST (-ifa).
    """
    for db in resfinder_db pointfinder_db disinfinder_db; do
        mkdir \$db
        wget -q -O \$db.tar.gz https://bitbucket.org/genomicepidemiology/\$db/get/${revision}.tar.gz
        tar -xzf \$db.tar.gz -C \$db --strip-components 1
        rm \$db.tar.gz
        if [ ! -s \$db/config ]; then
            echo "ERROR: \$db has no config file" >&2
            exit 1
        fi
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        resfinder_db: \$(cat resfinder_db/VERSION)
        pointfinder_db: \$(cat pointfinder_db/VERSION)
        disinfinder_db: \$(cat disinfinder_db/VERSION)
    END_VERSIONS
    """

    stub:
    """
    mkdir resfinder_db pointfinder_db disinfinder_db
    touch resfinder_db/config pointfinder_db/config disinfinder_db/config

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        resfinder_db: ${revision}
    END_VERSIONS
    """
}
