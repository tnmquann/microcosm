#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// Authors: @tnmquann (co-author)
// Maintainers: @tnmquann (co-maintainer)

include { KLEBORATE_V3_RUN as KLEBORATE_TYPING } from '../../../modules/local/kleborate/v3run/main'

// Kleborate v3 preset and module list per GTDB lineage (playground/modular/3_kleb.sh).
// Order matters: species complexes first, then Escherichia, then any other
// Klebsiella/Raoultella genome with the generic Enterobacterales modules.
def kleborateMode(lineage) {
    def text = lineage ?: ''
    def kpsc = ['pneumoniae', 'quasipneumoniae', 'variicola', 'quasivariicola', 'africana']
    def kosc = ['oxytoca', 'grimontii', 'michiganensis', 'pasteurii', 'huaxiensis', 'spallanzanii']
    if (kpsc.any { name -> text.contains("s__Klebsiella ${name}") }) {
        return ['kpsc', 'general__contig_stats']
    }
    if (kosc.any { name -> text.contains("s__Klebsiella ${name}") }) {
        return ['kosc', 'general__contig_stats']
    }
    if (text.contains('g__Escherichia')) {
        return ['escherichia', 'general__contig_stats']
    }
    if (text.contains('g__Klebsiella') || text.contains('g__Raoultella')) {
        return ['', 'enterobacterales__species,general__contig_stats']
    }
    return null
}

workflow KLEBSIELLA_TYPING {
    take:
    bacterial_genomes_ch

    main:
    ch_kleborate_input = bacterial_genomes_ch
        .map { meta, genome -> tuple(meta, kleborateMode(meta.lineage), genome) }
        .filter { _meta, mode, _genome -> mode != null }
        .map { meta, mode, genome ->
            def meta_new = meta.clone()
            meta_new.kleborate_preset = mode[0]
            tuple(meta_new, mode[0], mode[1], genome)
        }

    KLEBORATE_TYPING(ch_kleborate_input)

    emit:
    report = KLEBORATE_TYPING.out.txt
    summary = KLEBORATE_TYPING.out.txt
    versions = KLEBORATE_TYPING.out.versions_kleborate
}
