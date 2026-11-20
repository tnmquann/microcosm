#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { SALMONELLA_TYPING } from '../main'

// With --genomes_dir (real run) the SISTR serovar routes Typhi CT18 to
// Mykrobe/GenoTyphi, Paratyphi A ATCC 9150 to Paratype and Typhimurium LT2 to
// neither; K. pneumoniae never enters. -stub-run checks wiring only.
workflow {
    def common = [
        run_id: 'SAL_TEST', group: 'SAL_TEST', assembler: 'external',
        lineage_source: 'provided', source_route: 'standalone',
        classification_status: 'provided', domain: 'd__Bacteria'
    ]
    def ent = 'd__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae'
    def cases = [
        ['styphi_CT18', "${ent};g__Salmonella;s__Salmonella enterica"],
        ['sparatyphiA_ATCC9150', "${ent};g__Salmonella;s__Salmonella enterica"],
        ['styphimurium_LT2', "${ent};g__Salmonella;s__Salmonella enterica"],
        ['kpneumoniae_HS11286', "${ent};g__Klebsiella;s__Klebsiella pneumoniae"]
    ]
    ch_genomes = channel.fromList(cases).map { name, lineage ->
        def genome = params.genomes_dir
            ? file("${params.genomes_dir}/${name}.fasta", checkIfExists: true)
            : file("${projectDir}/../../../../test_data/genome_input/bacterial.fa", checkIfExists: true)
        tuple(common + [id: name, sample_id: name, genome_id: name, lineage: lineage], genome)
    }

    SALMONELLA_TYPING(ch_genomes)
    SALMONELLA_TYPING.out.routing.view { meta, routing ->
        def row = routing.splitCsv(header: true, sep: '\t')[0]
        "SALMONELLA_ROUTING ${meta.genome_id} sistr=${row.sistr_serovar} seqsero2=${row.seqsero2_prediction} genotyphi=${row.route_genotyphi} paratype=${row.route_paratype}"
    }
    SALMONELLA_TYPING.out.genotyphi_tsv.view { meta, tsv -> "SALMONELLA_GENOTYPHI ${meta.genome_id} ${tsv.name}" }
    SALMONELLA_TYPING.out.paratype_report.view { meta, tsv -> "SALMONELLA_PARATYPE ${meta.genome_id} ${tsv.name}" }
}
