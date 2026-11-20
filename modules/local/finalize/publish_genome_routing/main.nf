process PUBLISH_GENOME_ROUTING {
    tag "genome routing reports"
    label 'process_low'
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/52/52ccce28d2ab928ab862e25aae26314d69c8e38bd41ca9431c67ef05221348aa/data'
        : 'community.wave.seqera.io/library/coreutils_grep_gzip_lbzip2_pruned:838ba80435a629f8'}"

    input:
    path bacterial_report
    path excluded_report
    path unresolved_report

    output:
    path "bacterial_routing_manifest.tsv", emit: bacterial
    path "excluded_genomes.tsv", emit: excluded
    path "unresolved_genomes.tsv", emit: unresolved
    path "versions.yml", emit: versions

    script:
    """
    cp ${bacterial_report} bacterial_routing_manifest.tsv
    cp ${excluded_report} excluded_genomes.tsv
    cp ${unresolved_report} unresolved_genomes.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        coreutils: \$(cp --version | head -n 1 | awk '{print \$NF}')
    END_VERSIONS
    """

    stub:
    """
    cp ${bacterial_report} bacterial_routing_manifest.tsv
    cp ${excluded_report} excluded_genomes.tsv
    cp ${unresolved_report} unresolved_genomes.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        coreutils: "stub_version"
    END_VERSIONS
    """
}
