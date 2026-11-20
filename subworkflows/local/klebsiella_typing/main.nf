#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// Authors: @tnmquann (co-author)
// Maintainers: @tnmquann (co-maintainer)

include { KLEBORATE_V3_RUN as KLEBORATE_TYPING } from '../../../modules/local/kleborate/v3run/main'
include { readLineageMap; lookupLineage; stepWanted } from '../utils_bacterial_lineage/main'

// Kleborate v3 preset and module list per GTDB lineage (playground/modular/3_kleb.sh);
// the table is --kleborate_preset_map (default assets/kleborate_presets.tsv).
workflow KLEBSIELLA_TYPING {
    take:
    bacterial_genomes_ch

    main:
    pipeline_root = params.pipeline_root ?: projectDir
    preset_rows = readLineageMap(
        params.kleborate_preset_map ?: "${pipeline_root}/assets/kleborate_presets.tsv",
        ['preset', 'modules'],
        ['modules']
    )

    ch_kleborate_input = bacterial_genomes_ch
        .map { meta, genome ->
            tuple(meta, lookupLineage(preset_rows, meta.lineage), genome)
        }
        // kleborate=y in --bacterial_decision_table still needs a preset row (Kleborate needs its module list)
        .filter { meta, row, _genome -> row != null && stepWanted(meta, 'kleborate') }
        .map { meta, row, genome ->
            def meta_new = meta.clone()
            meta_new.kleborate_preset = row.preset?.trim() ?: ''
            tuple(meta_new, meta_new.kleborate_preset, row.modules.trim(), genome)
        }

    KLEBORATE_TYPING(ch_kleborate_input)

    emit:
    report = KLEBORATE_TYPING.out.txt
    summary = KLEBORATE_TYPING.out.txt
    versions = KLEBORATE_TYPING.out.versions_kleborate
}
