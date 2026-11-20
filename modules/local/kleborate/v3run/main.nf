process KLEBORATE_V3_RUN {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/kleborate_ezclermont_python_biopython:fe7e9ceb6272d003' :
        'community.wave.seqera.io/library/kleborate_ezclermont_python_biopython:facf2befd09bd8ec' }"

    input:
    tuple val(meta), val(preset), val(module_list), path(assembly)

    output:
    tuple val(meta), path("*.txt"), emit: txt, optional: true
    tuple val("${task.process}"), val('kleborate'), eval("kleborate --version 2>&1 | grep -oE '[0-9]+(\\.[0-9]+)+' | head -1"), emit: versions_kleborate, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def preset_arg = preset ? "--preset ${preset}" : ''
    """
    kleborate \\
        --outdir kleborate_results \\
        --assemblies ${assembly} \\
        ${preset_arg} \\
        --modules ${module_list} \\
        ${args}

    # Kleborate v3 names its table after the preset/modules; keep exactly one.
    shopt -s nullglob
    results=(kleborate_results/*_output.txt)
    if [ \${#results[@]} -gt 1 ]; then
        echo "ERROR: ambiguous Kleborate outputs: \${results[*]}" >&2
        exit 1
    elif [ \${#results[@]} -eq 1 ]; then
        mv "\${results[0]}" ${prefix}.txt
    else
        echo "Kleborate output was not found for ${prefix}" >&2
    fi
    rm -rf kleborate_results
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    printf 'strain\\n%s\\n' "${meta.id}" > ${prefix}.txt
    """
}
