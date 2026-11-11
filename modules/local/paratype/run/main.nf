process PARATYPE_RUN {
    tag "${meta.sample_id}:${meta.genome_id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path(fasta)
    path paratype_bundle

    output:
    tuple val(meta), path("*.paratype.txt"), emit: report
    path 'versions.yml', emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.genome_id}"
    def args = task.ext.args ?: ''
    """
    python ${paratype_bundle}/paratype.py \\
        --id "${meta.genome_id}" \\
        --mode fasta \\
        --fasta ${fasta} \\
        --threads ${task.cpus} \\
        --output "${prefix}.paratype.txt" \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        paratype: user_supplied_bundle
        python: \$(python --version 2>&1 | sed 's/Python //')
        bwa: \$(bwa 2>&1 | awk '/Version/{print \$2; exit}')
        samtools: \$(samtools --version | awk 'NR==1{print \$2}')
        bcftools: \$(bcftools --version | awk 'NR==1{print \$2}')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.genome_id}"
    """
    cat <<-END_REPORT > ${prefix}.paratype.txt
    sample_id\tgenotype\tstatus
    ${meta.genome_id}\tstub\tnot_a_biological_result
    END_REPORT

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        paratype: user_supplied_bundle_stub
    END_VERSIONS
    """
}
