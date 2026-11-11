#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { FILTER_BACTERIAL_BLAST_HITS } from '../main'

workflow {
    meta = [
        id: 'BLAST_FILTER_001',
        sample_id: 'BLAST_FILTER_001',
        genome_id: 'BLAST_FILTER_001',
        source_route: 'standalone'
    ]
    FILTER_BACTERIAL_BLAST_HITS(
        channel.of(tuple(meta, file("${projectDir}/data/hits.tsv", checkIfExists: true))),
        90,
        80,
        1,
        file(
            "${params.pipeline_root}/bin/py_scripts/assembly_based/filter_blasthits.py",
            checkIfExists: true
        )
    )
    FILTER_BACTERIAL_BLAST_HITS.out.filtered_hits.view { output_meta, output ->
        "FILTERED_BLAST ${output_meta.genome_id} ${output.name}"
    }
}
