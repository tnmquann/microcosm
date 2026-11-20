#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { VIRULENCEFINDER_RUN } from '../main'

workflow {
    VIRULENCEFINDER_RUN(
        channel.of(tuple(
            [id: 'VIRULENCE_TEST', sample_id: 'VIRULENCE_TEST', genome_id: 'VIRULENCE_TEST'],
            file(
                "${projectDir}/../../../../test_data/genome_input/bacterial.fa",
                checkIfExists: true
            )
        )),
        file(
            "${projectDir}/../../../../test_data/general_typing/virulencefinder",
            checkIfExists: true
        )
    )

    VIRULENCEFINDER_RUN.out.tsv.view { meta, report ->
        "VIRULENCEFINDER_LOCAL_TEST ${meta.id} ${report.name}"
    }
    VIRULENCEFINDER_RUN.out.versions.view { versions ->
        "VIRULENCEFINDER_LOCAL_VERSIONS ${versions.name}"
    }
}
