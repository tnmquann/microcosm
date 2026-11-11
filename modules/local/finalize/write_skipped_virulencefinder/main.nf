process WRITE_SKIPPED_VIRULENCEFINDER {
    tag "${meta.sample_id}:${meta.genome_id}"
    label 'process_single'

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("${meta.genome_id}.virulencefinder.skipped.standard.json"), emit: json
    tuple val(meta), path("${meta.genome_id}.virulencefinder.skipped.legacy.json"), emit: legacy_json
    tuple val(meta), path("${meta.genome_id}.virulencefinder.skipped.tsv"), emit: tsv
    tuple val(meta), path("${meta.genome_id}.virulencefinder.skipped.txt"), emit: txt
    tuple val(meta), path("${meta.genome_id}.virulencefinder.skipped-hit_in_genome_seq.fsa"), emit: genome_seq
    tuple val(meta), path("${meta.genome_id}.virulencefinder.skipped-virulence_genes.fsa"), emit: virulence_seq
    path "versions.yml", emit: versions

    script:
    """
    cat > "${meta.genome_id}.virulencefinder.skipped.standard.json" <<-END_JSON
    {"status":"skipped","reason":"no --virulencefinder_db supplied and no automatic VirulenceFinder database setup is configured"}
    END_JSON
    cp "${meta.genome_id}.virulencefinder.skipped.standard.json" "${meta.genome_id}.virulencefinder.skipped.legacy.json"
    printf 'status\\treason\\n' > "${meta.genome_id}.virulencefinder.skipped.tsv"
    printf 'skipped\\tno --virulencefinder_db supplied and no automatic VirulenceFinder database setup is configured\\n' >> "${meta.genome_id}.virulencefinder.skipped.tsv"
    printf 'VirulenceFinder skipped: no --virulencefinder_db supplied and no automatic VirulenceFinder database setup is configured.\\n' > "${meta.genome_id}.virulencefinder.skipped.txt"
    touch "${meta.genome_id}.virulencefinder.skipped-hit_in_genome_seq.fsa"
    touch "${meta.genome_id}.virulencefinder.skipped-virulence_genes.fsa"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        write_skipped_virulencefinder: 1.0.0
    END_VERSIONS
    """
}
