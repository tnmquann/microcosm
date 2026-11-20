/*
 * Local adaptation of nf-core/modules hamronization/summarize.
 *
 * The command contract is unchanged; the local environment pins Python 3.11
 * for hAMRonization 1.1.9 and the process emits a regular versions.yml.
 */
process HAMRONIZATION_SUMMARIZE_LOCAL {
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/hamronization:1.1.9--pyhdfd78af_0'
        : 'quay.io/biocontainers/hamronization:1.1.9--pyhdfd78af_0'}"

    input:
    path reports
    val format

    output:
    path "hamronization_combined_report.json", optional: true, emit: json
    path "hamronization_combined_report.tsv", optional: true, emit: tsv
    path "hamronization_combined_report.html", optional: true, emit: html
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def outformat = format == 'interactive' ? 'html' : format
    if (!(format in ['tsv', 'json', 'interactive'])) {
        error("Unsupported hAMRonization summary format: ${format}")
    }

    """
    hamronize \\
        summarize \\
        ${reports.join(' ')} \\
        -t ${format} \\
        ${args} \\
        -o hamronization_combined_report.${outformat}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hamronization: \$(hamronize --version 2>&1 | sed 's/hamronize //')
    END_VERSIONS
    """

    stub:
    def outformat = format == 'interactive' ? 'html' : format
    if (!(format in ['tsv', 'json', 'interactive'])) {
        error("Unsupported hAMRonization summary format: ${format}")
    }

    """
    touch hamronization_combined_report.${outformat}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hamronization: "1.1.9"
    END_VERSIONS
    """
}
