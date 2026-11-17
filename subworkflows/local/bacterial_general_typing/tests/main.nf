#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { BACTERIAL_GENERAL_TYPING } from '../main'

// -stub-run checks wiring; --genomes_dir runs the real tools on <name>.fasta.
workflow {
    def common = [
        run_id: 'GEN_TEST', group: 'GEN_TEST', assembler: 'external',
        lineage_source: 'provided', source_route: 'standalone',
        classification_status: 'provided', domain: 'd__Bacteria'
    ]
    def ent = 'd__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae'
    def cases = [
        ['ecoli_K12', "${ent};g__Escherichia;s__Escherichia coli"],
        ['kpneumoniae_HS11286', "${ent};g__Klebsiella;s__Klebsiella pneumoniae"]
    ]
    ch_genomes = channel.fromList(cases).map { name, lineage ->
        def genome = params.genomes_dir
            ? file("${params.genomes_dir}/${name}.fasta", checkIfExists: true)
            : file("${projectDir}/../../../../test_data/genome_input/bacterial.fa", checkIfExists: true)
        tuple(common + [id: name, sample_id: name, genome_id: name, lineage: lineage], genome)
    }

    BACTERIAL_GENERAL_TYPING(ch_genomes)
    BACTERIAL_GENERAL_TYPING.out.mlst_tsv.view { meta, tsv -> "GENERAL_MLST ${meta.genome_id} ${tsv.name} ${tsv.readLines() ? tsv.readLines()[-1] : 'empty'}" }
    BACTERIAL_GENERAL_TYPING.out.virulencefinder_tsv.view { meta, tsv -> "GENERAL_VIRULENCEFINDER ${meta.genome_id} ${tsv.name} rows=${Math.max(tsv.readLines().size() - 1, 0)}" }
    BACTERIAL_GENERAL_TYPING.out.plasmidfinder_tsv.view { meta, tsv -> "GENERAL_PLASMIDFINDER ${meta.genome_id} ${tsv.name} rows=${Math.max(tsv.readLines().size() - 1, 0)}" }
}
