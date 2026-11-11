#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { PARATYPE_RUN } from '../main'

workflow {
    PARATYPE_RUN(
        channel.of(tuple(
            [
                id: 'PARATYPE_TEST',
                sample_id: 'PARATYPE_TEST',
                genome_id: 'PARATYPE_TEST'
            ],
            file(
                "${params.pipeline_root}/test_data/genome_input/bacterial.fa",
                checkIfExists: true
            )
        )),
        file(
            "${params.pipeline_root}/test_data/general_typing/paratype_bundle",
            checkIfExists: true
        )
    )

    PARATYPE_RUN.out.report.view { meta, report ->
        "PARATYPE_RUN_STUB ${meta.sample_id} ${report.name}"
    }
    PARATYPE_RUN.out.versions.view { versions ->
        "PARATYPE_RUN_STUB_VERSIONS ${versions.name}"
    }
}
