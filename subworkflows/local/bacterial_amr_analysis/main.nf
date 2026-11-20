#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// Authors: @tnmquann (co-author)
// Maintainers: @tnmquann (co-maintainer)

include { GUNZIP as GUNZIP_AMR } from '../../../modules/nf-core/gunzip/main'
include { AMRFINDERPLUS_UPDATE_LOCAL as AMRFINDERPLUS_UPDATE } from '../../../modules/local/amrfinderplus/update/main'
include { AMRFINDERPLUS_RUN as AMRFINDERPLUS_AMR } from '../../../modules/nf-core/amrfinderplus/run/main'
include { ABRICATE_RUN as ABRICATE_AMR } from '../../../modules/nf-core/abricate/run/main'
include { RGI_PREPARECARDDB as RGI_PREPARECARDDB_AMR } from '../../../modules/local/rgi/preparecarddb/main'
include { RGI_MAIN_LOCAL as RGI_AMR } from '../../../modules/local/rgi/main/main'
include { RESFINDER_DOWNLOADDB_LOCAL } from '../../../modules/local/resfinder/downloaddb/main'
include { RESFINDER_RUN_LOCAL as RESFINDER_AMR } from '../../../modules/local/resfinder/run/main'
include { DEEPARG_DOWNLOADDATA_LOCAL as DEEPARG_DOWNLOADDATA_AMR } from '../../../modules/local/deeparg/downloaddata/main'
include { DEEPARG_PREDICT as DEEPARG_AMR } from '../../../modules/nf-core/deeparg/predict/main'
include { HAMRONIZATION_REPORT_LOCAL as HAMRONIZATION_ABRICATE_AMR } from '../../../modules/local/hamronization/report/main'
include { HAMRONIZATION_REPORT_LOCAL as HAMRONIZATION_AMRFINDERPLUS_AMR } from '../../../modules/local/hamronization/report/main'
include { HAMRONIZATION_REPORT_LOCAL as HAMRONIZATION_DEEPARG_AMR } from '../../../modules/local/hamronization/report/main'
include { HAMRONIZATION_REPORT_LOCAL as HAMRONIZATION_RESFINDER_AMR } from '../../../modules/local/hamronization/report/main'
include { HAMRONIZATION_REPORT_LOCAL as HAMRONIZATION_RGI_AMR } from '../../../modules/local/hamronization/report/main'
include { HAMRONIZATION_SUMMARIZE_LOCAL as HAMRONIZATION_SUMMARIZE_AMR } from '../../../modules/local/hamronization/summarize/main'
include { ARGNORM_LOCAL as ARGNORM_SUMMARY_AMR } from '../../../modules/local/argnorm/main'
include { ARGNORM_LOCAL as ARGNORM_TOOL_AMR } from '../../../modules/local/argnorm/main'
include { readLineageMap; lookupLineage; stepWanted } from '../utils_bacterial_lineage/main'

// Steps of playground/modular/5_amr.sh (AMRFinderPlus, ABRicate, RGI, ResFinder,
// DeepARG) and 6_hamr.sh (hAMRonization per tool + summary, ARGnorm) for every
// d__Bacteria genome. Each tool has a --skip_<tool> switch; a database parameter
// that is not set is downloaded (only for tools that run).
workflow BACTERIAL_AMR_ANALYSIS {
    take:
    bacterial_genomes_ch

    main:
    ch_versions = channel.empty()
    ch_hamronized = channel.empty()
    ch_amrfinder_report = channel.empty()
    ch_abricate_report = channel.empty()
    ch_rgi_tsv = channel.empty()
    ch_rgi_json = channel.empty()
    ch_resfinder_json = channel.empty()
    ch_resfinder_results = channel.empty()
    ch_deeparg_arg = channel.empty()
    ch_deeparg_potential = channel.empty()
    ch_hamronization_summary = channel.empty()
    ch_normalized = channel.empty()
    ch_normalized_summary = channel.empty()
    pipeline_root = params.pipeline_root ?: projectDir

    // AMRFinderPlus --organism from --amrfinderplus_organism_map (default
    // assets/amrfinderplus_organisms.tsv, the 5_amr.sh table).
    organism_rows = params.skip_amrfinderplus ? [] : readLineageMap(
        params.amrfinderplus_organism_map ?: "${pipeline_root}/assets/amrfinderplus_organisms.tsv",
        ['organism'],
        ['organism']
    )

    // --bacterial_decision_table: only genomes with at least one AMR step
    ch_genome_splits = bacterial_genomes_ch
        .filter { meta, _genome -> ['amrfinderplus', 'abricate', 'rgi', 'resfinder', 'deeparg'].any { step -> stepWanted(meta, step) } }
        .map { meta, genome ->
            def meta_new = meta.clone()
            // --bacterial_decision_table amrfinderplus_organism wins over --amrfinderplus_organism_map
            meta_new.amrfinder_organism = meta.decision_params?.amrfinderplus_organism
                ?: lookupLineage(organism_rows, meta.lineage)?.organism?.trim() ?: ''
            tuple(meta_new, genome)
        }
        .branch { _meta, genome ->
            compressed: genome.getName().endsWith('.gz')
            uncompressed: true
        }
    GUNZIP_AMR(ch_genome_splits.compressed)
    ch_amr_genomes = GUNZIP_AMR.out.gunzip.mix(ch_genome_splits.uncompressed)
    ch_versions = ch_versions.mix(GUNZIP_AMR.out.versions_gunzip)

    // ---- 5_amr.sh: one task per genome and tool; 6_hamr.sh: hAMRonization per tool ----
    // Database downloads take the first genome as trigger: nothing is downloaded without one.
    if (!params.skip_amrfinderplus) {
        ch_amrfinderplus_genomes = ch_amr_genomes.filter { meta, _genome -> stepWanted(meta, 'amrfinderplus') }
        if (params.amrfinderplus_db) {
            ch_amrfinder_db = channel.value(file(params.amrfinderplus_db, checkIfExists: true))
        } else {
            AMRFINDERPLUS_UPDATE(ch_amrfinderplus_genomes.first().map { _genome -> true })
            ch_amrfinder_db = AMRFINDERPLUS_UPDATE.out.db
            ch_versions = ch_versions.mix(AMRFINDERPLUS_UPDATE.out.versions_amrfinder)
        }
        AMRFINDERPLUS_AMR(ch_amrfinderplus_genomes, ch_amrfinder_db)
        ch_amrfinder_report = AMRFINDERPLUS_AMR.out.report
        HAMRONIZATION_AMRFINDERPLUS_AMR(
            AMRFINDERPLUS_AMR.out.report, 'amrfinderplus', 'tsv',
            AMRFINDERPLUS_AMR.out.tool_version.first(), AMRFINDERPLUS_AMR.out.db_version.first(), true
        )
        ch_hamronized = ch_hamronized.mix(
            HAMRONIZATION_AMRFINDERPLUS_AMR.out.tsv.map { meta, tsv -> tuple(meta, 'amrfinderplus', tsv) }
        )
        ch_versions = ch_versions
            .mix(AMRFINDERPLUS_AMR.out.versions_amrfinderplus)
            .mix(AMRFINDERPLUS_AMR.out.versions_amrfinderplus_database)
            .mix(HAMRONIZATION_AMRFINDERPLUS_AMR.out.versions)
    }

    if (!params.skip_abricate) {
        ch_abricate_genomes = ch_amr_genomes.filter { meta, _genome -> stepWanted(meta, 'abricate') }
        ABRICATE_AMR(ch_abricate_genomes, [])
        ch_abricate_report = ABRICATE_AMR.out.report
        HAMRONIZATION_ABRICATE_AMR(
            ABRICATE_AMR.out.report, 'abricate', 'tsv',
            ABRICATE_AMR.out.versions_abricate.map { versions -> versions[2] }.first(),
            params.abricate_db_version, false
        )
        ch_hamronized = ch_hamronized.mix(
            HAMRONIZATION_ABRICATE_AMR.out.tsv.map { meta, tsv -> tuple(meta, 'abricate', tsv) }
        )
        ch_versions = ch_versions
            .mix(ABRICATE_AMR.out.versions_abricate)
            .mix(HAMRONIZATION_ABRICATE_AMR.out.versions)
    }

    if (!params.skip_rgi) {
        ch_rgi_genomes = ch_amr_genomes.filter { meta, _genome -> stepWanted(meta, 'rgi') }
        if (params.rgi_preparecarddb_dir) {
            ch_card_db = channel.value(file(params.rgi_preparecarddb_dir, checkIfExists: true))
        } else {
            RGI_PREPARECARDDB_AMR(ch_rgi_genomes.first().map { _genome -> [] })
            ch_card_db = RGI_PREPARECARDDB_AMR.out.db
            ch_versions = ch_versions.mix(RGI_PREPARECARDDB_AMR.out.versions)
        }
        RGI_AMR(ch_rgi_genomes, ch_card_db, [])
        ch_rgi_tsv = RGI_AMR.out.tsv
        ch_rgi_json = RGI_AMR.out.json
        HAMRONIZATION_RGI_AMR(
            RGI_AMR.out.tsv, 'rgi', 'tsv',
            RGI_AMR.out.tool_version.first(), RGI_AMR.out.db_version.first(), true
        )
        ch_hamronized = ch_hamronized.mix(
            HAMRONIZATION_RGI_AMR.out.tsv.map { meta, tsv -> tuple(meta, 'rgi', tsv) }
        )
        ch_versions = ch_versions
            .mix(RGI_AMR.out.versions)
            .mix(HAMRONIZATION_RGI_AMR.out.versions)
    }

    if (!params.skip_resfinder) {
        ch_resfinder_genomes = ch_amr_genomes.filter { meta, _genome -> stepWanted(meta, 'resfinder') }
        if (params.resfinder_db && params.pointfinder_db && params.disinfinder_db) {
            ch_resfinder_db = channel.value(file(params.resfinder_db, checkIfExists: true))
            ch_pointfinder_db = channel.value(file(params.pointfinder_db, checkIfExists: true))
            ch_disinfinder_db = channel.value(file(params.disinfinder_db, checkIfExists: true))
        } else {
            RESFINDER_DOWNLOADDB_LOCAL(
                ch_resfinder_genomes.first().map { _genome -> params.resfinder_db_revision }
            )
            ch_versions = ch_versions.mix(RESFINDER_DOWNLOADDB_LOCAL.out.versions)
            ch_resfinder_db = params.resfinder_db
                ? channel.value(file(params.resfinder_db, checkIfExists: true))
                : RESFINDER_DOWNLOADDB_LOCAL.out.resfinder_db
            ch_pointfinder_db = params.pointfinder_db
                ? channel.value(file(params.pointfinder_db, checkIfExists: true))
                : RESFINDER_DOWNLOADDB_LOCAL.out.pointfinder_db
            ch_disinfinder_db = params.disinfinder_db
                ? channel.value(file(params.disinfinder_db, checkIfExists: true))
                : RESFINDER_DOWNLOADDB_LOCAL.out.disinfinder_db
        }
        RESFINDER_AMR(ch_resfinder_genomes, ch_resfinder_db, ch_pointfinder_db, ch_disinfinder_db)
        ch_resfinder_json = RESFINDER_AMR.out.json
        ch_resfinder_results = RESFINDER_AMR.out.results
        // ResFinder JSON carries its own software/database versions.
        HAMRONIZATION_RESFINDER_AMR(RESFINDER_AMR.out.json, 'resfinder', 'tsv', '', '', false)
        ch_hamronized = ch_hamronized.mix(
            HAMRONIZATION_RESFINDER_AMR.out.tsv.map { meta, tsv -> tuple(meta, 'resfinder', tsv) }
        )
        ch_versions = ch_versions
            .mix(RESFINDER_AMR.out.versions_resfinder)
            .mix(HAMRONIZATION_RESFINDER_AMR.out.versions)
    }

    if (!params.skip_deeparg) {
        ch_deeparg_genomes = ch_amr_genomes.filter { meta, _genome -> stepWanted(meta, 'deeparg') }
        if (params.deeparg_db) {
            ch_deeparg_db = channel.value(file(params.deeparg_db, checkIfExists: true))
        } else {
            DEEPARG_DOWNLOADDATA_AMR(ch_deeparg_genomes.first().map { _genome -> true })
            ch_deeparg_db = DEEPARG_DOWNLOADDATA_AMR.out.db
            ch_versions = ch_versions.mix(DEEPARG_DOWNLOADDATA_AMR.out.versions_deeparg)
        }
        DEEPARG_AMR(
            ch_deeparg_genomes.map { meta, genome -> tuple(meta, genome, params.deeparg_model) },
            ch_deeparg_db
        )
        ch_deeparg_arg = DEEPARG_AMR.out.arg
        ch_deeparg_potential = DEEPARG_AMR.out.potential_arg
        HAMRONIZATION_DEEPARG_AMR(
            DEEPARG_AMR.out.potential_arg, 'deeparg', 'tsv',
            DEEPARG_AMR.out.versions_deeparg.map { versions -> versions[2] }.first(),
            params.deeparg_db_version, true
        )
        ch_hamronized = ch_hamronized.mix(
            HAMRONIZATION_DEEPARG_AMR.out.tsv.map { meta, tsv -> tuple(meta, 'deeparg', tsv) }
        )
        ch_versions = ch_versions
            .mix(DEEPARG_AMR.out.versions_deeparg)
            .mix(HAMRONIZATION_DEEPARG_AMR.out.versions)
    }

    HAMRONIZATION_SUMMARIZE_AMR(
        ch_hamronized.map { _meta, _tool, tsv -> tsv }.collect(),
        'tsv'
    )
    ch_hamronization_summary = HAMRONIZATION_SUMMARIZE_AMR.out.tsv
    ch_versions = ch_versions.mix(HAMRONIZATION_SUMMARIZE_AMR.out.versions)

    // ---- 6_hamr.sh: ARGnorm on hAMRonized output; RGI is not supported ----
    // Header-only tables are skipped: ARGnorm fails on them (pandas EmptyDataError).
    if (!params.skip_argnorm) {
        ARGNORM_SUMMARY_AMR(
            HAMRONIZATION_SUMMARIZE_AMR.out.tsv
                .filter { tsv -> tsv.readLines().size() > 1 }
                .map { tsv -> tuple([id: 'amr_summary'], tsv) },
            'hamronization', ''
        )
        ARGNORM_TOOL_AMR(
            ch_hamronized
                .filter { _meta, tool, tsv -> tool != 'rgi' && tsv.readLines().size() > 1 }
                .map { meta, tool, tsv ->
                    def meta_new = meta.clone()
                    meta_new.amr_tool = tool
                    tuple(meta_new, tsv)
                },
            'hamronization', ''
        )
        ch_normalized = ARGNORM_TOOL_AMR.out.tsv
        ch_normalized_summary = ARGNORM_SUMMARY_AMR.out.tsv
        ch_versions = ch_versions
            .mix(ARGNORM_SUMMARY_AMR.out.versions)
            .mix(ARGNORM_TOOL_AMR.out.versions)
    }

    emit:
    amrfinder_report = ch_amrfinder_report
    abricate_report = ch_abricate_report
    rgi_tsv = ch_rgi_tsv
    rgi_json = ch_rgi_json
    resfinder_json = ch_resfinder_json
    resfinder_results = ch_resfinder_results
    deeparg_arg = ch_deeparg_arg
    deeparg_potential = ch_deeparg_potential
    hamronized = ch_hamronized.map { meta, _tool, tsv -> tuple(meta, tsv) }
    hamronization_summary = ch_hamronization_summary
    normalized = ch_normalized
    normalized_summary = ch_normalized_summary
    versions = ch_versions
}
