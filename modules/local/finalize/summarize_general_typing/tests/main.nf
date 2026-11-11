#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { SUMMARIZE_GENERAL_TYPING } from '../main'

workflow {
    meta = [
        id: 'GENERAL_SUMMARY_001',
        sample_id: 'GENERAL_SUMMARY_001',
        genome_id: 'GENERAL_SUMMARY_001',
        source_route: 'standalone'
    ]
    input_ch = channel.of(
        tuple(
            meta,
            file("${projectDir}/data/virulence.tsv", checkIfExists: true),
            file("${projectDir}/data/plasmid.tsv", checkIfExists: true),
            file("${projectDir}/data/mlst.tsv", checkIfExists: true),
            file("${projectDir}/data/blast.tsv", checkIfExists: true)
        )
    )
    SUMMARIZE_GENERAL_TYPING(
        input_ch,
        file(
            "${params.pipeline_root}/bin/py_scripts/assembly_based/summarize_gene_hits.py",
            checkIfExists: true
        )
    )
    SUMMARIZE_GENERAL_TYPING.out.evidence.view { output_meta, output ->
        "GENERAL_TYPING_EVIDENCE ${output_meta.genome_id} ${output.name}"
    }
}
