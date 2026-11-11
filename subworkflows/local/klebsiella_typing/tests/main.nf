#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { KLEBSIELLA_TYPING } from '../main'

workflow {
    def common = [
        run_id: 'WP4_TEST', group: 'WP4_TEST', assembler: 'external',
        lineage_source: 'provided', source_route: 'standalone',
        classification_status: 'provided', domain: 'd__Bacteria'
    ]
    ch_genomes = channel.of(
        tuple(common + [
            id: 'KLEB_001', sample_id: 'KLEB_001', genome_id: 'KLEB_001',
            genus: 'g__Klebsiella', species: 's__Klebsiella pneumoniae'
        ], file("${projectDir}/../../../../test_data/genome_input/bacterial.fa", checkIfExists: true)),
        tuple(common + [
            id: 'ECOLI_001', sample_id: 'ECOLI_001', genome_id: 'ECOLI_001',
            genus: 'g__Escherichia', species: 's__Escherichia coli'
        ], file("${projectDir}/../../../../test_data/genome_input/missing_lineage.fa", checkIfExists: true))
    )

    KLEBSIELLA_TYPING(ch_genomes)
    KLEBSIELLA_TYPING.out.report.view { meta, report ->
        "KLEBORATE_REPORT ${meta.sample_id} ${meta.genus} ${report.name}"
    }
    KLEBSIELLA_TYPING.out.versions.collect().view { records ->
        "KLEBORATE_VERSION_RECORDS ${records.size()}"
    }
}
