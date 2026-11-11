#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { VALIDATE_BACTERIAL_TYPING_DATABASES } from '../main'

workflow {
    def data_root = "${params.pipeline_root}/test_data/general_typing"
    VALIDATE_BACTERIAL_TYPING_DATABASES(
        true,
        file("${data_root}/virulencefinder", checkIfExists: true),
        file("${data_root}/blast_reference.fa", checkIfExists: true),
        [],
        tuple(
            'unit_test_reference',
            90,
            80,
            file(
                "${params.pipeline_root}/bin/py_scripts/assembly_based/validate_general_typing_databases.py",
                checkIfExists: true
            )
        )
    )
    VALIDATE_BACTERIAL_TYPING_DATABASES.out.validation.view { output ->
        "VALIDATED_GENERAL_TYPING_DATABASES ${output.name}"
    }
}
