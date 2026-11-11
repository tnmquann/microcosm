#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { VALIDATE_GENOME_INPUT } from '../main'

workflow {
    VALIDATE_GENOME_INPUT(
        file("${params.pipeline_root}/test_data/genome_input/valid_provided.csv", checkIfExists: true),
        false,
        file("${params.pipeline_root}/bin/py_scripts/assembly_based/validate_genome_input.py", checkIfExists: true)
    )
    VALIDATE_GENOME_INPUT.out.validated_csv.view { output ->
        "VALIDATED_GENOME_INPUT ${output.name}"
    }
}
