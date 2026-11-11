#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { VALIDATE_BACTERIAL_AMR_DATABASES } from '../main'

workflow {
    def db_root = "${params.pipeline_root}/test_data/amr_databases"
    VALIDATE_BACTERIAL_AMR_DATABASES(
        true,
        "${db_root}/amrfinder",
        "${db_root}/abricate",
        'ncbi',
        "${db_root}/rgi",
        true,
        "${db_root}/resfinder",
        "${db_root}/pointfinder",
        "${db_root}/deeparg",
        '1',
        file("${params.pipeline_root}/bin/py_scripts/assembly_based/validate_amr_databases.py", checkIfExists: true)
    )
    VALIDATE_BACTERIAL_AMR_DATABASES.out.validation.view { output ->
        "VALIDATED_AMR_DATABASES ${output.name}"
    }
}
