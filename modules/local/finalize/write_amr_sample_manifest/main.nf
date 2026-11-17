process WRITE_AMR_SAMPLE_MANIFEST {
    tag "${meta.sample_id}:${meta.genome_id}"
    label 'process_low'

    // Standard-library Python helper.
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.11' :
        'quay.io/biocontainers/python:3.11' }"

    input:
    tuple val(meta), path(amrfinder_report), path(abricate_report), path(rgi_report), path(resfinder_report), path(deeparg_report)
    val abricate_normalizable
    path manifest_script
    path database_validation

    output:
    tuple val(meta), path("${meta.genome_id}.amr_manifest.tsv"), emit: manifest
    path "versions.yml", emit: versions

    script:
    def pointfinder_arg = meta.pointfinder_species ?: ''
    def normalizable_flag = abricate_normalizable ? '--abricate-normalizable' : ''
    """
    python ${manifest_script} \
        --sample-id "${meta.sample_id}" \
        --genome-id "${meta.genome_id}" \
        --source-route "${meta.source_route}" \
        --amrfinder-report ${amrfinder_report} \
        --abricate-report ${abricate_report} \
        --rgi-report ${rgi_report} \
        --resfinder-report ${resfinder_report} \
        --deeparg-report ${deeparg_report} \
        --pointfinder-species "${pointfinder_arg}" \
        ${normalizable_flag} \
        --database-validation ${database_validation} \
        --output "${meta.genome_id}.amr_manifest.tsv"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        write_amr_manifest.py: \$(python ${manifest_script} --version)
        python: \$(python --version 2>&1 | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    def pointfinder_arg = meta.pointfinder_species ?: ''
    def normalizable_flag = abricate_normalizable ? '--abricate-normalizable' : ''
    """
    python ${manifest_script} \
        --sample-id "${meta.sample_id}" \
        --genome-id "${meta.genome_id}" \
        --source-route "${meta.source_route}" \
        --amrfinder-report ${amrfinder_report} \
        --abricate-report ${abricate_report} \
        --rgi-report ${rgi_report} \
        --resfinder-report ${resfinder_report} \
        --deeparg-report ${deeparg_report} \
        --pointfinder-species "${pointfinder_arg}" \
        ${normalizable_flag} \
        --database-validation ${database_validation} \
        --output "${meta.genome_id}.amr_manifest.tsv"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        write_amr_manifest.py: \$(python ${manifest_script} --version)
        python: \$(python --version 2>&1 | sed 's/Python //')
    END_VERSIONS
    """
}
