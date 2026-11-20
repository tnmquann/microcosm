process VIRULENCEFINDER_RUN {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/virulencefinder%3A3.2.1--pyhdfd78af_0' :
        'quay.io/biocontainers/virulencefinder:3.2.1--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(assembly)
    path db

    output:
    tuple val(meta), path("*.standard.json"), emit: json
    tuple val(meta), path("*.legacy.json"), emit: legacy_json
    tuple val(meta), path("*.tsv"), emit: tsv
    tuple val(meta), path("*.txt"), emit: txt
    tuple val(meta), path("*-hit_in_genome_seq.fsa"), emit: genome_seq
    tuple val(meta), path("*-virulence_genes.fsa"), emit: virulence_seq
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def is_compressed = assembly.getName().endsWith('.gz')
    def fasta_name = is_compressed ? assembly.getName().replaceFirst(/\.gz$/, '') : assembly
    """
    if [ "${is_compressed}" = "true" ]; then
        gzip -c -d ${assembly} > ${fasta_name}
    fi

    mkdir virulencefinder_results
    python -m virulencefinder \
        -ifa ${fasta_name} \
        -o virulencefinder_results \
        -j ${prefix}.standard.json \
        -p ${db} \
        -d all \
        -x \
        -q \
        ${args}

    mv virulencefinder_results/data.json ${prefix}.legacy.json
    mv virulencefinder_results/results_tab.tsv ${prefix}.tsv
    mv virulencefinder_results/results.txt ${prefix}.txt
    mv virulencefinder_results/Hit_in_genome_seq.fsa ${prefix}-hit_in_genome_seq.fsa
    mv virulencefinder_results/Virulence_genes.fsa ${prefix}-virulence_genes.fsa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        virulencefinder: \$(python -m virulencefinder --version 2>&1)
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo '{}' > ${prefix}.standard.json
    echo '{}' > ${prefix}.legacy.json
    printf 'Database\\tVirulence factor\\n' > ${prefix}.tsv
    printf 'VirulenceFinder Results\\n\\nNo hit found\\n' > ${prefix}.txt
    touch ${prefix}-hit_in_genome_seq.fsa
    touch ${prefix}-virulence_genes.fsa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        virulencefinder: 3.2.1
    END_VERSIONS
    """
}
