process SYLPH_POSTPROCESS {
    tag "${prefix}"
    label 'process_single'

    conda "conda-forge::python=3.11"
    // Standard-library Python helper.
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.11' :
        'quay.io/biocontainers/python:3.11' }"

    input:
    tuple val(sample_ids), path(sylphmpa, stageAs: 'sylphtax/??/*')
    val prefix
    val datatype
    val phyloseq
    path postprocess_script

    output:
    path "results/${prefix}_sylphtax_${datatype}.tsv", emit: merged
    path "results/phyloseq", emit: phyloseq, optional: true
    tuple val("${task.process}"), val('python'), eval("python --version 2>&1 | sed 's/Python //'"), emit: versions_python, topic: versions

    script:
    def manifest_args = [sample_ids, [sylphmpa].flatten()].transpose()
        .collect { sample_id, report -> "'${sample_id}' '${report}'" }
        .join(' ')
    """
    printf '%s\\t%s\\n' ${manifest_args} > manifest.tsv

    python ${postprocess_script} \\
        --manifest manifest.tsv \\
        --datatype ${datatype} \\
        --prefix ${prefix} \\
        --outdir results \\
        ${phyloseq ? '--phyloseq' : ''}
    """

    stub:
    """
    mkdir -p results
    printf 'clade_name\\n' > results/${prefix}_sylphtax_${datatype}.tsv
    ${phyloseq ? 'mkdir -p results/phyloseq && touch results/phyloseq/tax_table.csv results/phyloseq/otu_table.csv results/phyloseq/sample_data.csv' : ''}
    """
}
