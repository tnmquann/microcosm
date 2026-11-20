#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { KLEBSIELLA_TYPING } from '../main'

// Lineage routing check. With -stub-run the genomes are fixtures; for a real
// run pass --genomes_dir with <name>.fasta files (see tests/nextflow.config).
workflow {
    def common = [
        run_id: 'KLEB_TEST', group: 'KLEB_TEST', assembler: 'external',
        lineage_source: 'provided', source_route: 'standalone',
        classification_status: 'provided', domain: 'd__Bacteria'
    ]
    def ent = 'd__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae'
    def cases = [
        ['kpneumoniae_HS11286', "${ent};g__Klebsiella;s__Klebsiella pneumoniae"],
        ['koxytoca_KCTC1686', "${ent};g__Klebsiella;s__Klebsiella michiganensis"],
        ['ecoli_K12', "${ent};g__Escherichia;s__Escherichia coli"],
        ['kaerogenes_KCTC2190', "${ent};g__Klebsiella;s__Klebsiella aerogenes"],
        ['styphi_CT18', "${ent};g__Salmonella;s__Salmonella enterica"]
    ]
    ch_genomes = channel.fromList(cases).map { name, lineage ->
        def genome = params.genomes_dir
            ? file("${params.genomes_dir}/${name}.fasta", checkIfExists: true)
            : file("${projectDir}/../../../../test_data/genome_input/bacterial.fa", checkIfExists: true)
        tuple(common + [id: name, sample_id: name, genome_id: name, lineage: lineage], genome)
    }

    KLEBSIELLA_TYPING(ch_genomes)
    KLEBSIELLA_TYPING.out.report.view { meta, report ->
        "KLEBORATE_REPORT ${meta.genome_id} preset=${meta.kleborate_preset ?: 'none'} ${report.name}"
    }
}
