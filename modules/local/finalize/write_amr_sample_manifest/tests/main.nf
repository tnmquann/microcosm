#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { WRITE_AMR_SAMPLE_MANIFEST } from '../main'

workflow {
    meta = [
        id: 'sample',
        sample_id: 'sample',
        genome_id: 'genome',
        source_route: 'standalone',
        pointfinder_species: 'Escherichia coli'
    ]
    reports = channel.of(tuple(
        meta,
        file("${projectDir}/data/amrfinder.tsv", checkIfExists: true),
        file("${projectDir}/data/abricate.txt", checkIfExists: true),
        file("${projectDir}/data/rgi.txt", checkIfExists: true),
        file("${projectDir}/data/resfinder.json", checkIfExists: true),
        file("${projectDir}/data/deeparg.mapping.ARG", checkIfExists: true)
    ))
    WRITE_AMR_SAMPLE_MANIFEST(
        reports,
        true,
        file("${params.pipeline_root}/bin/py_scripts/assembly_based/write_amr_manifest.py", checkIfExists: true),
        file("${projectDir}/data/amr_database_validation.tsv", checkIfExists: true)
    )
    WRITE_AMR_SAMPLE_MANIFEST.out.manifest.view { outputMeta, output ->
        "WROTE_AMR_MANIFEST ${outputMeta.genome_id} ${output.name}"
    }
}
