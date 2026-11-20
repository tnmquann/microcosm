#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { ARGNORM_LOCAL } from '../main'

workflow {
    ARGNORM_LOCAL(
        channel.of(tuple([id: 'sample'], file("${projectDir}/data/input.tsv", checkIfExists: true))),
        'amrfinderplus',
        'ncbi'
    )
    ARGNORM_LOCAL.out.tsv.view { meta, output ->
        "ARGNORM_LOCAL_TEST ${meta.id} ${output.name}"
    }
}
