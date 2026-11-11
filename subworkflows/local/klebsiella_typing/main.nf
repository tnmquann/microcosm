#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// Authors: @tnmquann (co-author)
// Maintainers: @tnmquann (co-maintainer)

include { KLEBORATE as KLEBORATE_TYPING } from '../../../modules/nf-core/kleborate/main'

workflow KLEBSIELLA_TYPING {
    take:
    bacterial_genomes_ch

    main:
    ch_klebsiella_genomes = bacterial_genomes_ch
        .filter { meta, _genome ->
            meta.domain == 'd__Bacteria' && meta.genus == 'g__Klebsiella'
        }
        .map { meta, genome ->
            def meta_new = meta.clone()
            meta_new.klebsiella_typing_key = "${meta.sample_id}::${meta.genome_id}"
            tuple(meta_new, [genome])
        }

    KLEBORATE_TYPING(ch_klebsiella_genomes)

    emit:
    report = KLEBORATE_TYPING.out.txt
    summary = KLEBORATE_TYPING.out.txt
    versions = KLEBORATE_TYPING.out.versions_kleborate
}
