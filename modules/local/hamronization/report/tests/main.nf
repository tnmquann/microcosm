#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { HAMRONIZATION_REPORT_LOCAL } from '../main'

workflow {
    HAMRONIZATION_REPORT_LOCAL(
        channel.of(
            tuple(
                [id: 'HAM_REPORT_TEST'],
                file("${projectDir}/data/input.tsv", checkIfExists: true)
            )
        ),
        'abricate',
        'tsv',
        '1.0.1',
        '2021-Mar-27',
        false
    )

    HAMRONIZATION_REPORT_LOCAL.out.tsv.view { meta, report ->
        "HAM_REPORT_LOCAL_TEST ${meta.id} ${report.name}"
    }
    HAMRONIZATION_REPORT_LOCAL.out.versions.view { versions ->
        "HAM_REPORT_LOCAL_VERSIONS ${versions.name}"
    }
}
