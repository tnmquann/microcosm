#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { GENOME_INPUT_AND_BACTERIAL_ROUTING } from '../main'

workflow {
    test_data_dir = file(params.test_data_dir, checkIfExists: true)
    assembly_genomes_ch = channel.of(
        tuple(
            [
                id: 'ASSEMBLY_SAMPLE',
                run_id: 'RUN_ASSEMBLY',
                group: 'GROUP_ASSEMBLY',
                assembler: 'SPAdes',
                binner: 'MetaBAT2',
                bin_id: 'assembly_genome'
            ],
            file("${test_data_dir}/assembly_genome.fa", checkIfExists: true)
        )
    )
    assembly_classifications_ch = channel.of(
        tuple(
            [id: 'ASSEMBLY_SAMPLE'],
            file("${test_data_dir}/assembly_classification.csv", checkIfExists: true)
        )
    )

    GENOME_INPUT_AND_BACTERIAL_ROUTING(
        assembly_genomes_ch,
        assembly_classifications_ch
    )

    GENOME_INPUT_AND_BACTERIAL_ROUTING.out.bacterial_genomes.view { meta, genome ->
        "BACTERIAL ${meta.sample_id} ${meta.genome_id} ${meta.domain} ${meta.genus} ${genome.name}"
    }
}
