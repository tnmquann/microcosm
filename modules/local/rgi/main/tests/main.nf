#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { RGI_MAIN_LOCAL } from '../main'

workflow {
    RGI_MAIN_LOCAL(
        channel.of(
            tuple(
                [id: 'RGI_TEST'],
                file(
                    "${projectDir}/../../../../../test_data/genome_input/bacterial.fa",
                    checkIfExists: true
                )
            )
        ),
        file(
            "${projectDir}/../../../../../test_data/amr_databases/rgi",
            checkIfExists: true
        ),
        channel.value([])
    )

    RGI_MAIN_LOCAL.out.tsv.view { meta, report ->
        "RGI_LOCAL_TEST ${meta.id} ${report.name}"
    }
    RGI_MAIN_LOCAL.out.versions.view { versions ->
        "RGI_LOCAL_VERSIONS ${versions.name}"
    }
}
