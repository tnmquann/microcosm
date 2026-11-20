#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { PUBLISH_GENOME_ROUTING } from '../main'

workflow {
    PUBLISH_GENOME_ROUTING(
        file("${projectDir}/data/bacterial.tsv", checkIfExists: true),
        file("${projectDir}/data/excluded.tsv", checkIfExists: true),
        file("${projectDir}/data/unresolved.tsv", checkIfExists: true)
    )
    PUBLISH_GENOME_ROUTING.out.bacterial.view { output ->
        "PUBLISHED_GENOME_ROUTING ${output.name}"
    }
}
