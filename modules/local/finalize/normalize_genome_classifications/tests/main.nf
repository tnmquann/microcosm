#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { NORMALIZE_GENOME_CLASSIFICATIONS } from '../main'

workflow {
    classifications = channel.of([
        file("${params.pipeline_root}/subworkflows/local/genome_input_bacterial_routing/tests/data/assembly_classification.csv", checkIfExists: true)
    ])
    NORMALIZE_GENOME_CLASSIFICATIONS(
        classifications,
        file("${params.pipeline_root}/bin/py_scripts/bacterial_analysis/normalize_genome_classifications.py", checkIfExists: true)
    )
    NORMALIZE_GENOME_CLASSIFICATIONS.out.resolved_tsv.view { output ->
        "NORMALIZED_GENOME_CLASSIFICATIONS ${output.name}"
    }
}
