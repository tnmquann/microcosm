process GTDBTK_TAXGENOME {
    tag "${meta.id}"
    label 'process_single'

    conda "conda-forge::python=3.11"
    // Standard-library Python helper.
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.11' :
        'quay.io/biocontainers/python:3.11' }"

    input:
    tuple val(meta), path(summaries)
    path converter_script

    output:
    tuple val(meta), path("*.gtdbtk_classification.csv"), emit: genome_classification
    path "versions.yml", emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    python ${converter_script} --summary ${summaries} --output ${prefix}.gtdbtk_classification.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gtdbtk_to_taxgenome.py: \$(python ${converter_script} --version)
        python: \$(python --version 2>&1 | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "query_name,status,rank,lineage" > ${prefix}.gtdbtk_classification.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gtdbtk_to_taxgenome.py: \$(python ${converter_script} --version)
        python: \$(python --version 2>&1 | sed 's/Python //')
    END_VERSIONS
    """
}
