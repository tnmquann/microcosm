#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { HAMRONIZATION_SUMMARIZE_LOCAL } from '../main'

workflow {
    HAMRONIZATION_SUMMARIZE_LOCAL(
        channel.of(
            [
                file("${projectDir}/data/report.tsv", checkIfExists: true)
            ]
        ),
        'tsv'
    )

    HAMRONIZATION_SUMMARIZE_LOCAL.out.tsv.view { report ->
        "HAM_SUMMARIZE_LOCAL_TEST ${report.name}"
    }
    HAMRONIZATION_SUMMARIZE_LOCAL.out.versions.view { versions ->
        "HAM_SUMMARIZE_LOCAL_VERSIONS ${versions.name}"
    }
}
