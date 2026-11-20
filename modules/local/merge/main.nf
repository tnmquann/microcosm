process MERGE_PAIREDENDSEQS {
    tag { meta.id }
    label 'process_low'
    
    conda "conda-forge::pigz=2.3.4"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/52/52ccce28d2ab928ab862e25aae26314d69c8e38bd41ca9431c67ef05221348aa/data'
        : 'community.wave.seqera.io/library/coreutils_grep_gzip_lbzip2_pruned:838ba80435a629f8'}"

    input:
    tuple val(meta), path(reads)  // reads is now a list of two files

    output:
    tuple val(meta), path("${meta.id}/${meta.id}.fastq.gz"), emit: merged_seqs
    path "versions.yml", emit: versions

    script:
    def prefix = meta.id
    """
    mkdir -p ${meta.id}
    cat ${reads.join(' ')} > ${meta.id}/${prefix}.fastq.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cat: \$(cat --version 2>&1 | head -n 1 | sed 's/^.* //')
    END_VERSIONS
    """

    stub:
    def prefix = meta.id
    """
    mkdir -p ${meta.id}
    touch ${meta.id}/${prefix}.fastq.gz
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cat: "stub_version"
    END_VERSIONS
    """
}