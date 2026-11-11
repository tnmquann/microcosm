#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { AMRFINDERPLUS_RUN_LOCAL } from '../main'

workflow {
    AMRFINDERPLUS_RUN_LOCAL(
        channel.of(
            tuple(
                [id: 'AMRFINDER_TEST', organism: 'Escherichia'],
                file(
                    "${projectDir}/../../../../../test_data/genome_input/bacterial.fa",
                    checkIfExists: true
                )
            )
        ),
        file(
            "${projectDir}/../../../../../test_data/amr_databases/amrfinder",
            checkIfExists: true
        )
    )

    AMRFINDERPLUS_RUN_LOCAL.out.report.view { meta, report ->
        "AMRFINDER_LOCAL_TEST ${meta.id} ${report.name}"
    }
    AMRFINDERPLUS_RUN_LOCAL.out.mutation_report.view { meta, report ->
        "AMRFINDER_LOCAL_MUTATION_TEST ${meta.id} ${report.name}"
    }
    AMRFINDERPLUS_RUN_LOCAL.out.versions.view { versions ->
        "AMRFINDER_LOCAL_VERSIONS ${versions.name}"
    }
}
