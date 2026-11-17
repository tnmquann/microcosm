process PARATYPE_RUN {
    tag "${meta.sample_id}:${meta.genome_id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/python_biopython_bwa_samtools_pruned:d893df7a7ad0a807' :
        'community.wave.seqera.io/library/python_biopython_bwa_samtools_pruned:07d8f98513940f99' }"

    input:
    tuple val(meta), path(fasta)
    path paratype_bundle

    output:
    tuple val(meta), path("*.tsv"), emit: report
    path 'versions.yml', emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.genome_id}"
    def args = task.ext.args ?: ''
    """
    # paratype.py finds its reference next to its own realpath and runs
    # 'bwa index' there: use a task-local copy so the bundle is never written.
    cp -rL ${paratype_bundle} paratype_local
    paratype_revision=\$(git -C paratype_local rev-parse HEAD 2>/dev/null || echo user_supplied_bundle)

    python paratype_local/paratype.py \\
        --id "${meta.genome_id}" \\
        --mode fasta \\
        --fasta ${fasta} \\
        --threads ${task.cpus} \\
        --output "${prefix}.tsv" \\
        ${args}

    # paratype.py runs bwa/samtools/bcftools through os.system and ignores failures.
    if [ ! -s "${prefix}.tsv" ]; then
        echo "ERROR: Paratype wrote no result for ${meta.genome_id}" >&2
        exit 1
    fi
    rm -rf paratype_local

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        paratype: \${paratype_revision}
        python: \$(python --version 2>&1 | sed 's/Python //')
        bwa: \$(bwa 2>&1 | awk '/Version/{print \$2; exit}')
        samtools: \$(samtools --version | awk 'NR==1{print \$2}')
        bcftools: \$(bcftools --version | awk 'NR==1{print \$2}')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.genome_id}"
    """
    printf 'sample_id\\tgenotype\\tstatus\\n%s\\tstub\\tnot_a_biological_result\\n' "${meta.genome_id}" > ${prefix}.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        paratype: stub
    END_VERSIONS
    """
}
