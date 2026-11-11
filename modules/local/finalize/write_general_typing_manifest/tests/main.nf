#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { WRITE_GENERAL_TYPING_MANIFEST } from '../main'

workflow {
    meta = [
        id: 'sample',
        sample_id: 'sample',
        genome_id: 'genome',
        source_route: 'standalone'
    ]
    reports = channel.of(tuple(
        meta,
        file("${projectDir}/data/virulencefinder.tsv", checkIfExists: true),
        file("${projectDir}/data/plasmidfinder.tsv", checkIfExists: true),
        file("${projectDir}/data/mlst.tsv", checkIfExists: true),
        file("${projectDir}/data/blastn.txt", checkIfExists: true)
    ))
    WRITE_GENERAL_TYPING_MANIFEST(
        reports,
        'user_supplied_reference_fasta',
        'unit_test_reference',
        file(
            "${params.pipeline_root}/bin/py_scripts/assembly_based/write_general_typing_manifest.py",
            checkIfExists: true
        )
    )
    WRITE_GENERAL_TYPING_MANIFEST.out.manifest.view { outputMeta, output ->
        "WROTE_GENERAL_TYPING_MANIFEST ${outputMeta.genome_id} ${output.name}"
    }
}
