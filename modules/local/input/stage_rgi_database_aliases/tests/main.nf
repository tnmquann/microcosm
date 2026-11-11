#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { STAGE_RGI_DATABASE_ALIASES } from '../main'

workflow {
    STAGE_RGI_DATABASE_ALIASES(
        "${params.pipeline_root}/test_data/amr_databases/rgi",
        file("${projectDir}/data/amr_databases.ready", checkIfExists: true)
    )
    STAGE_RGI_DATABASE_ALIASES.out.card.view { output ->
        "STAGED_RGI_DATABASE ${output.name}"
    }
}
