process BIN_MULTIQC_TABLES {
    tag "bins"
    label 'process_single'

    conda "conda-forge::python=3.11"
    // Standard-library Python helper.
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.11' :
        'quay.io/biocontainers/python:3.11' }"

    input:
    path taxgenome, stageAs: 'taxgenome/??/*'
    path qc_summaries, stageAs: 'qc/*'
    path quast_summary, stageAs: 'quast/*'
    path multiqc_script

    output:
    path "bin_*_mqc.json", emit: tables, optional: true
    path "report_bins.tsv", emit: quast_report, optional: true
    tuple val("${task.process}"), val('python'), eval("python --version 2>&1 | sed 's/Python //'"), emit: versions_python, topic: versions

    script:
    def args = task.ext.args ?: ''
    // [x].flatten(): a single staged file is a Path, and Path.join() would join its name parts
    def tax_args = taxgenome ? "--taxgenome ${[taxgenome].flatten().join(' ')}" : ''
    def qc_args = qc_summaries ? "--qc ${[qc_summaries].flatten().join(' ')}" : ''
    def quast_args = quast_summary ? "--quast ${[quast_summary].flatten().join(' ')}" : ''
    """
    python ${multiqc_script} ${tax_args} ${qc_args} ${quast_args} --outdir . ${args}
    """

    stub:
    """
    touch bin_summary_mqc.json
    """
}
