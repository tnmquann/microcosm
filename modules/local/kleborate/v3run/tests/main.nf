#!/usr/bin/env nextflow
// Real-data check: params.genomes_dir holds <name>.fasta files.
include { KLEBORATE_V3_RUN } from '../main'

workflow {
    def cases = [
        ['kpneumoniae_HS11286', 'kpsc', 'general__contig_stats'],
        ['koxytoca_KCTC1686', 'kosc', 'general__contig_stats'],
        ['ecoli_K12', 'escherichia', 'general__contig_stats'],
        ['kaerogenes_KCTC2190', '', 'enterobacterales__species,general__contig_stats']
    ]
    KLEBORATE_V3_RUN(channel.fromList(cases).map { name, preset, mods ->
        tuple([id: name, genome_id: name], preset, mods, file("${params.genomes_dir}/${name}.fasta", checkIfExists: true))
    })
    KLEBORATE_V3_RUN.out.txt.view { meta, txt -> "KLEBORATE_V3 ${meta.id} ${txt.name}" }
    KLEBORATE_V3_RUN.out.versions_kleborate.first().view { v -> "KLEBORATE_V3_VERSION ${v[2]}" }
}
