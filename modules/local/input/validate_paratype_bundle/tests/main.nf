#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { VALIDATE_PARATYPE_BUNDLE } from '../main'

workflow {
    VALIDATE_PARATYPE_BUNDLE(
        file(
            "${params.pipeline_root}/test_data/general_typing/paratype_bundle",
            checkIfExists: true
        ),
        file(
            "${params.pipeline_root}/bin/py_scripts/assembly_based/validate_paratype_bundle.py",
            checkIfExists: true
        )
    )

    VALIDATE_PARATYPE_BUNDLE.out.ready.view { ready ->
        "VALIDATED_PARATYPE_BUNDLE ${ready.name}"
    }
    VALIDATE_PARATYPE_BUNDLE.out.versions.view { versions ->
        "VALIDATED_PARATYPE_BUNDLE_VERSIONS ${versions.name}"
    }
}
