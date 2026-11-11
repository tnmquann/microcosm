#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { BACTERIAL_AMR_ANALYSIS } from '../main'

workflow {
    def common = [
        run_id: 'RUN_AMR',
        group: 'GROUP_AMR',
        assembler: 'external',
        lineage_source: 'provided',
        source_route: 'standalone',
        classification_status: 'provided',
        domain: 'd__Bacteria'
    ]
    ch_genomes = channel.of(
        tuple(
            common + [
                id: 'ECOLI_001',
                sample_id: 'ECOLI_001',
                genome_id: 'ECOLI_001',
                lineage: 'd__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae;g__Escherichia;s__Escherichia coli',
                genus: 'g__Escherichia',
                species: 's__Escherichia coli'
            ],
            file("${projectDir}/../../../../test_data/genome_input/bacterial.fa", checkIfExists: true)
        ),
        tuple(
            common + [
                id: 'PSEUDO_001',
                sample_id: 'PSEUDO_001',
                genome_id: 'PSEUDO_001',
                lineage: 'd__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria;o__Pseudomonadales;f__Pseudomonadaceae;g__Pseudomonas;s__Pseudomonas aeruginosa',
                genus: 'g__Pseudomonas',
                species: 's__Pseudomonas aeruginosa'
            ],
            file("${projectDir}/../../../../test_data/genome_input/missing_lineage.fa", checkIfExists: true)
        )
    )

    BACTERIAL_AMR_ANALYSIS(ch_genomes)

    BACTERIAL_AMR_ANALYSIS.out.manifests.view { meta, manifest ->
        "AMR_MANIFEST ${meta.sample_id} ${meta.genome_id} ${meta.pointfinder_species ?: 'SKIPPED'} ${manifest.name}"
    }
    BACTERIAL_AMR_ANALYSIS.out.database_validation.view { validation ->
        "AMR_DATABASES ${validation.name}"
    }
    BACTERIAL_AMR_ANALYSIS.out.versions.collect().view { records ->
        "AMR_VERSION_RECORDS ${records.size()}"
    }
}
