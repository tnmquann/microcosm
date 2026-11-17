process RESFINDER_RUN_LOCAL {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/resfinder:4.7.2--pyhdfd78af_0' :
        'quay.io/biocontainers/resfinder:4.7.2--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(fasta)
    path db_res
    path db_point
    path db_disinf

    output:
    tuple val(meta), path("${prefix}.json"), emit: json
    tuple val(meta), path("${prefix}"), emit: results
    tuple val("${task.process}"), val('resfinder'), eval("python -m resfinder --version 2>&1 | tail -1"), emit: versions_resfinder, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}.resfinder"
    // PointFinder species key as in 5_amr.sh: GTDB species, else genus, lower case with '_'.
    def species_key = (meta.species ?: '').toString().replaceFirst('^s__', '').trim().toLowerCase().replaceAll('\\s+', '_')
    def genus_key = (meta.genus ?: '').toString().replaceFirst('^g__', '').trim().toLowerCase().replaceAll('\\s+', '_')
    """
    point_args=""
    if [ -n "${species_key}" ] && [ -d "${db_point}/${species_key}" ]; then
        point_args="--point --species ${species_key}"
    elif [ -n "${genus_key}" ] && [ -d "${db_point}/${genus_key}" ]; then
        point_args="--point --species ${genus_key}"
    fi

    python -m resfinder \\
        -ifa ${fasta} \\
        -o ${prefix} \\
        -j ${prefix}.json \\
        --kma_threads ${task.cpus} \\
        --acquired \\
        --disinfectant \\
        \$point_args \\
        --db_path_point ${db_point} \\
        --db_path_disinf ${db_disinf} \\
        --db_path_res ${db_res} \\
        ${args}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}.resfinder"
    """
    mkdir ${prefix}
    touch ${prefix}/ResFinder_results_tab.txt
    echo '{"type": "software_result", "software_name": "ResFinder", "seq_regions": {}, "seq_variations": {}, "phenotypes": {}}' > ${prefix}.json
    """
}
