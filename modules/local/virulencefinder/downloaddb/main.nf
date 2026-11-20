process VIRULENCEFINDER_DOWNLOADDB_LOCAL {
    tag "virulencefinder_db ${revision}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/virulencefinder%3A3.2.1--pyhdfd78af_0' :
        'quay.io/biocontainers/virulencefinder:3.2.1--pyhdfd78af_0' }"

    input:
    val revision

    output:
    path 'virulencefinder_db', emit: db
    path 'versions.yml', emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Same source as the packaged download-virulence-db.sh, pinned to a revision.
    // No KMA index: assemblies are searched with BLAST (-ifa), as in 2_general_bac.sh.
    """
    mkdir virulencefinder_db
    wget -q -O virulencefinder_db.tar.gz \\
        https://bitbucket.org/genomicepidemiology/virulencefinder_db/get/${revision}.tar.gz
    tar -xzf virulencefinder_db.tar.gz -C virulencefinder_db --strip-components 1
    rm virulencefinder_db.tar.gz
    for required in config notes.txt; do
        if [ ! -s virulencefinder_db/\$required ]; then
            echo "ERROR: VirulenceFinder database lacks \$required" >&2
            exit 1
        fi
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        virulencefinder_db: ${revision}
    END_VERSIONS
    """

    stub:
    """
    mkdir virulencefinder_db
    printf 'virulence_ecoli\\tvirulence_ecoli\\tstub\\n' > virulencefinder_db/config
    touch virulencefinder_db/notes.txt virulencefinder_db/virulence_ecoli.fsa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        virulencefinder_db: ${revision}
    END_VERSIONS
    """
}
