process PUBLISH_GENOME_ROUTING {
    tag "genome routing reports"
    label 'process_low'

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
