process PARATYPE_DOWNLOAD_LOCAL {
    tag "Paratype ${revision}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/python_biopython_bwa_samtools_pruned:d893df7a7ad0a807' :
        'community.wave.seqera.io/library/python_biopython_bwa_samtools_pruned:07d8f98513940f99' }"

    input:
    val revision

    output:
    path 'paratype_bundle', emit: bundle
    path 'versions.yml', emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    git clone --quiet https://github.com/CHRF-Genomics/Paratype.git paratype_bundle
    git -C paratype_bundle checkout --quiet ${revision}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        paratype: \$(git -C paratype_bundle rev-parse HEAD)
    END_VERSIONS
    """

    stub:
    """
    mkdir paratype_bundle
    touch paratype_bundle/paratype.py paratype_bundle/SParatyphiAKU12601.fasta \\
        paratype_bundle/SParatyphiA_genotype_specific_alleles.txt \\
        paratype_bundle/SParatyphiA_gene_mutation_codons.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        paratype: ${revision}
    END_VERSIONS
    """
}
