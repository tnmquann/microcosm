process PARSE_SALMONELLA_TYPING {
    tag "${meta.sample_id}:${meta.genome_id}"
    label 'process_low'

    // Standard-library Python helper.
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.11' :
        'quay.io/biocontainers/python:3.11' }"

    input:
    tuple val(meta), path(seqsero2_tsv), path(sistr_tsv)
    path parser_script

    output:
    tuple val(meta), path("${meta.genome_id}.salmonella_evidence.tsv"), emit: evidence
    tuple val(meta), path("${meta.genome_id}.salmonella_routing.tsv"), emit: routing
    path "versions.yml", emit: versions

    script:
    """
    python ${parser_script} \\
        --sample-id "${meta.sample_id}" \\
        --genome-id "${meta.genome_id}" \\
        --seqsero2-tsv ${seqsero2_tsv} \\
        --sistr-tsv ${sistr_tsv} \\
        --evidence-output "${meta.genome_id}.salmonella_evidence.tsv" \\
        --routing-output "${meta.genome_id}.salmonella_routing.tsv"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        parse_salmonella_typing.py: \$(python ${parser_script} --version)
        python: \$(python --version 2>&1 | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    """
    cat <<-END_EVIDENCE > "${meta.genome_id}.salmonella_evidence.tsv"
    sample_id	genome_id	tool	predicted_serovar_or_lineage	normalized_prediction	confidence_or_status	raw_result_path
    ${meta.sample_id}	${meta.genome_id}	seqsero2			empty_stub	${seqsero2_tsv}
    ${meta.sample_id}	${meta.genome_id}	sistr			empty_stub	${sistr_tsv}
    END_EVIDENCE

    cat <<-END_ROUTING > "${meta.genome_id}.salmonella_routing.tsv"
    routing_key	sample_id	genome_id	sistr_serovar	seqsero2_prediction	route_genotyphi	route_paratype	tool_conflict	routing_reason
    ${meta.sample_id}::${meta.genome_id}	${meta.sample_id}	${meta.genome_id}			false	false	false	sistr_serovar_missing
    END_ROUTING

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        parse_salmonella_typing.py: \$(python ${parser_script} --version)
        python: \$(python --version 2>&1 | sed 's/Python //')
    END_VERSIONS
    """
}
