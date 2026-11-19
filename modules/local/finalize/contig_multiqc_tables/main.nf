process CONTIG_MULTIQC_TABLES {
    tag "contigs"
    label 'process_single'

    conda "conda-forge::python=3.11"
    // Standard-library Python helper.
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.11' :
        'quay.io/biocontainers/python:3.11' }"

    input:
    path seqkit, stageAs: 'seqkit/*'
    path coverage, stageAs: 'coverage/*'
    path faa, stageAs: 'faa/*'
    path plasmid, stageAs: 'plasmid/*'
    path virus, stageAs: 'virus/*'
    path tiara, stageAs: 'tiara/*'
    path skani, stageAs: 'skani/*'
    path multiqc_script

    output:
    path "contig_*_mqc.json", emit: tables, optional: true
    tuple val("${task.process}"), val('python'), eval("python --version 2>&1 | sed 's/Python //'"), emit: versions_python, topic: versions

    script:
    // [x].flatten(): a single staged file is a Path, and Path.join() would join its name parts
    def options = [['--seqkit', seqkit], ['--coverage', coverage], ['--faa', faa], ['--plasmid', plasmid], ['--virus', virus], ['--tiara', tiara], ['--skani', skani]]
        .findAll { option -> option[1] }
        .collect { option -> "${option[0]} ${[option[1]].flatten().join(' ')}" }
        .join(' ')
    """
    python ${multiqc_script} ${options} --outdir .
    """

    stub:
    """
    touch contig_summary_mqc.json
    """
}
