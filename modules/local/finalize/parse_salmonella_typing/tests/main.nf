#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { PARSE_SALMONELLA_TYPING } from '../main'

workflow {
    PARSE_SALMONELLA_TYPING(
        channel.of(tuple(
            [id: 'SALMONELLA_TEST', sample_id: 'SALMONELLA_TEST', genome_id: 'SALMONELLA_TEST'],
            file("${projectDir}/data/seqsero2_typhi.tsv", checkIfExists: true),
            file("${projectDir}/data/sistr_typhi.tsv", checkIfExists: true)
        )),
        file(
            "${projectDir}/../../../../../bin/py_scripts/bacterial_analysis/parse_salmonella_typing.py",
            checkIfExists: true
        )
    )
    PARSE_SALMONELLA_TYPING.out.routing.view { meta, routing ->
        "SALMONELLA_ROUTING_TEST ${meta.sample_id} ${routing.name}"
    }
}
