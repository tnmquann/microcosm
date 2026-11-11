#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { SALMONELLA_TYPING } from '../main'

workflow {
    def common = [
        run_id: 'WP4_TEST', group: 'WP4_TEST', assembler: 'external',
        lineage_source: 'provided', source_route: 'standalone',
        classification_status: 'provided', domain: 'd__Bacteria'
    ]
    ch_genomes = channel.of(
        tuple(common + [
            id: 'SALMONELLA_001', sample_id: 'SALMONELLA_001', genome_id: 'SALMONELLA_001',
            genus: 'g__Salmonella', species: 's__Salmonella enterica'
        ], file("${projectDir}/../../../../test_data/genome_input/bacterial.fa", checkIfExists: true)),
        tuple(common + [
            id: 'KLEB_001', sample_id: 'KLEB_001', genome_id: 'KLEB_001',
            genus: 'g__Klebsiella', species: 's__Klebsiella pneumoniae'
        ], file("${projectDir}/../../../../test_data/genome_input/missing_lineage.fa", checkIfExists: true)),
        tuple(common + [
            id: 'ECOLI_001', sample_id: 'ECOLI_001', genome_id: 'ECOLI_001',
            genus: 'g__Escherichia', species: 's__Escherichia coli'
        ], file("${projectDir}/../../../../test_data/genome_input/bacterial.fa", checkIfExists: true))
    )

    SALMONELLA_TYPING(ch_genomes)
    SALMONELLA_TYPING.out.evidence.view { meta, report ->
        "SALMONELLA_EVIDENCE ${meta.sample_id} ${report.name}"
    }
    SALMONELLA_TYPING.out.routing.view { meta, report ->
        "SALMONELLA_ROUTING ${meta.sample_id} ${report.name}"
    }
    SALMONELLA_TYPING.out.versions.collect().view { records ->
        "SALMONELLA_VERSION_RECORDS ${records.size()}"
    }
}
