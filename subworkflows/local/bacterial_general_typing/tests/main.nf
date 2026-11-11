#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { BACTERIAL_GENERAL_TYPING } from '../main'

workflow {
    def common = [
        run_id: 'RUN_TYPING',
        group: 'GROUP_TYPING',
        assembler: 'external',
        lineage_source: 'provided',
        source_route: 'standalone',
        classification_status: 'provided',
        domain: 'd__Bacteria'
    ]
    ch_genomes = channel.of(
        tuple(
            common + [
                id: 'ECOLI_TYPING_001',
                sample_id: 'ECOLI_TYPING_001',
                genome_id: 'ECOLI_TYPING_001',
                lineage: 'd__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae;g__Escherichia;s__Escherichia coli',
                genus: 'g__Escherichia',
                species: 's__Escherichia coli'
            ],
            file("${projectDir}/../../../../test_data/genome_input/bacterial.fa", checkIfExists: true)
        ),
        tuple(
            common + [
                id: 'PSEUDO_TYPING_001',
                sample_id: 'PSEUDO_TYPING_001',
                genome_id: 'PSEUDO_TYPING_001',
                lineage: 'd__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria;o__Pseudomonadales;f__Pseudomonadaceae;g__Pseudomonas;s__Pseudomonas aeruginosa',
                genus: 'g__Pseudomonas',
                species: 's__Pseudomonas aeruginosa'
            ],
            file("${projectDir}/../../../../test_data/genome_input/missing_lineage.fa", checkIfExists: true)
        )
    )

    BACTERIAL_GENERAL_TYPING(ch_genomes)

    BACTERIAL_GENERAL_TYPING.out.manifests.view { meta, manifest ->
        "GENERAL_TYPING_MANIFEST ${meta.sample_id} ${meta.genome_id} ${manifest.name}"
    }
    BACTERIAL_GENERAL_TYPING.out.database_validation.view { validation ->
        "GENERAL_TYPING_DATABASES ${validation.name}"
    }
    BACTERIAL_GENERAL_TYPING.out.versions.collect().view { records ->
        "GENERAL_TYPING_VERSION_RECORDS ${records.size()}"
    }
}
