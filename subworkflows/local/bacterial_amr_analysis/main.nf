#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// Authors: @tnmquann (co-author)
// Maintainers: @tnmquann (co-maintainer)

include { GUNZIP as GUNZIP_AMR } from '../../../modules/nf-core/gunzip/main'
include { AMRFINDERPLUS_UPDATE } from '../../../modules/nf-core/amrfinderplus/update/main'
include { AMRFINDERPLUS_RUN as AMRFINDERPLUS_AMR } from '../../../modules/nf-core/amrfinderplus/run/main'
include { ABRICATE_RUN as ABRICATE_AMR } from '../../../modules/nf-core/abricate/run/main'
include { RGI_PREPARECARDDB as RGI_PREPARECARDDB_AMR } from '../../../modules/local/rgi/preparecarddb/main'
include { RGI_MAIN_LOCAL as RGI_AMR } from '../../../modules/local/rgi/main/main'
include { RESFINDER_DOWNLOADDB_LOCAL } from '../../../modules/local/resfinder/downloaddb/main'
include { RESFINDER_RUN_LOCAL as RESFINDER_AMR } from '../../../modules/local/resfinder/run/main'
include { DEEPARG_DOWNLOADDATA as DEEPARG_DOWNLOADDATA_AMR } from '../../../modules/nf-core/deeparg/downloaddata/main'
include { DEEPARG_PREDICT as DEEPARG_AMR } from '../../../modules/nf-core/deeparg/predict/main'
include { HAMRONIZATION_REPORT_LOCAL as HAMRONIZATION_ABRICATE_AMR } from '../../../modules/local/hamronization/report/main'
include { HAMRONIZATION_REPORT_LOCAL as HAMRONIZATION_AMRFINDERPLUS_AMR } from '../../../modules/local/hamronization/report/main'
include { HAMRONIZATION_REPORT_LOCAL as HAMRONIZATION_DEEPARG_AMR } from '../../../modules/local/hamronization/report/main'
include { HAMRONIZATION_REPORT_LOCAL as HAMRONIZATION_RESFINDER_AMR } from '../../../modules/local/hamronization/report/main'
include { HAMRONIZATION_REPORT_LOCAL as HAMRONIZATION_RGI_AMR } from '../../../modules/local/hamronization/report/main'
include { HAMRONIZATION_SUMMARIZE_LOCAL as HAMRONIZATION_SUMMARIZE_AMR } from '../../../modules/local/hamronization/summarize/main'
include { ARGNORM_LOCAL as ARGNORM_SUMMARY_AMR } from '../../../modules/local/argnorm/main'
include { ARGNORM_LOCAL as ARGNORM_TOOL_AMR } from '../../../modules/local/argnorm/main'

// CGE tag shared by resfinder_db, pointfinder_db and disinfinder_db for ResFinder 4.7.2.
def resfinderDbRevision() {
    return 'resfinder-4.7.2'
}

// AMRFinderPlus --organism from the GTDB lineage (playground/modular/5_amr.sh).
def amrfinderOrganism(lineage) {
    def text = lineage ?: ''
    def species = [
        'Acinetobacter baumannii', 'Bordetella pertussis', 'Burkholderia cepacia',
        'Burkholderia mallei', 'Burkholderia pseudomallei', 'Citrobacter freundii',
        'Clostridioides difficile', 'Corynebacterium diphtheriae', 'Enterobacter asburiae',
        'Enterobacter cloacae', 'Enterococcus faecalis', 'Enterococcus faecium',
        'Haemophilus influenzae', 'Helicobacter pylori', 'Klebsiella oxytoca',
        'Klebsiella pneumoniae', 'Neisseria gonorrhoeae', 'Neisseria meningitidis',
        'Pseudomonas aeruginosa', 'Serratia marcescens', 'Staphylococcus aureus',
        'Staphylococcus epidermidis', 'Staphylococcus pseudintermedius',
        'Streptococcus agalactiae', 'Streptococcus pneumoniae', 'Streptococcus pyogenes',
        'Vibrio cholerae', 'Vibrio parahaemolyticus', 'Vibrio vulnificus'
    ]
    def genera = ['Campylobacter', 'Escherichia', 'Salmonella']
    def hit = species.find { name -> text.contains("s__${name}") }
    if (hit) {
        return hit.replace(' ', '_')
    }
    return genera.find { name -> text.contains("g__${name}") } ?: ''
}

// Steps of playground/modular/5_amr.sh (AMRFinderPlus, ABRicate, RGI, ResFinder,
// DeepARG) and 6_hamr.sh (hAMRonization per tool + summary, ARGnorm) for every
// d__Bacteria genome. A database parameter that is not set is downloaded.
workflow BACTERIAL_AMR_ANALYSIS {
    take:
    bacterial_genomes_ch

    main:
    ch_versions = channel.empty()

    ch_genome_splits = bacterial_genomes_ch
        .map { meta, genome ->
            def meta_new = meta.clone()
            meta_new.amrfinder_organism = amrfinderOrganism(meta.lineage)
            tuple(meta_new, genome)
        }
        .branch { _meta, genome ->
            compressed: genome.getName().endsWith('.gz')
            uncompressed: true
        }
    GUNZIP_AMR(ch_genome_splits.compressed)
    ch_amr_genomes = GUNZIP_AMR.out.gunzip.mix(ch_genome_splits.uncompressed)
    ch_versions = ch_versions.mix(GUNZIP_AMR.out.versions_gunzip)

    // ---- Databases: supplied, or downloaded into the run work directory ----
    if (params.amrfinderplus_db) {
        ch_amrfinder_db = channel.value(file(params.amrfinderplus_db, checkIfExists: true))
    } else {
        AMRFINDERPLUS_UPDATE()
        ch_amrfinder_db = AMRFINDERPLUS_UPDATE.out.db
        ch_versions = ch_versions.mix(AMRFINDERPLUS_UPDATE.out.versions_amrfinder)
    }
    if (params.rgi_preparecarddb_dir) {
        ch_card_db = channel.value(file(params.rgi_preparecarddb_dir, checkIfExists: true))
    } else {
        RGI_PREPARECARDDB_AMR([])
        ch_card_db = RGI_PREPARECARDDB_AMR.out.db
        ch_versions = ch_versions.mix(RGI_PREPARECARDDB_AMR.out.versions)
    }
    if (params.resfinder_db && params.pointfinder_db && params.disinfinder_db) {
        ch_resfinder_db = channel.value(file(params.resfinder_db, checkIfExists: true))
        ch_pointfinder_db = channel.value(file(params.pointfinder_db, checkIfExists: true))
        ch_disinfinder_db = channel.value(file(params.disinfinder_db, checkIfExists: true))
    } else {
        RESFINDER_DOWNLOADDB_LOCAL(
            ch_amr_genomes.first().map { _genome -> resfinderDbRevision() }
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
    if (params.deeparg_db) {
        ch_deeparg_db = channel.value(file(params.deeparg_db, checkIfExists: true))
    } else {
        DEEPARG_DOWNLOADDATA_AMR()
        ch_deeparg_db = DEEPARG_DOWNLOADDATA_AMR.out.db
        ch_versions = ch_versions.mix(DEEPARG_DOWNLOADDATA_AMR.out.versions_deeparg)
    }

    // ---- 5_amr.sh: one task per genome and tool ----
    AMRFINDERPLUS_AMR(ch_amr_genomes, ch_amrfinder_db)
    ABRICATE_AMR(ch_amr_genomes, [])
    RGI_AMR(ch_amr_genomes, ch_card_db, [])
    RESFINDER_AMR(ch_amr_genomes, ch_resfinder_db, ch_pointfinder_db, ch_disinfinder_db)
    DEEPARG_AMR(ch_amr_genomes.map { meta, genome -> tuple(meta, genome, 'SS') }, ch_deeparg_db)
    ch_versions = ch_versions
        .mix(AMRFINDERPLUS_AMR.out.versions_amrfinderplus)
        .mix(AMRFINDERPLUS_AMR.out.versions_amrfinderplus_database)
        .mix(ABRICATE_AMR.out.versions_abricate)
        .mix(RGI_AMR.out.versions)
        .mix(RESFINDER_AMR.out.versions_resfinder)
        .mix(DEEPARG_AMR.out.versions_deeparg)

    // ---- 6_hamr.sh: hAMRonization with the versions of the tools that ran ----
    HAMRONIZATION_ABRICATE_AMR(
        ABRICATE_AMR.out.report, 'abricate', 'tsv',
        ABRICATE_AMR.out.versions_abricate.map { versions -> versions[2] }.first(),
        params.abricate_db_version, false
    )
    HAMRONIZATION_AMRFINDERPLUS_AMR(
        AMRFINDERPLUS_AMR.out.report, 'amrfinderplus', 'tsv',
        AMRFINDERPLUS_AMR.out.tool_version.first(), AMRFINDERPLUS_AMR.out.db_version.first(), true
    )
    HAMRONIZATION_DEEPARG_AMR(
        DEEPARG_AMR.out.potential_arg, 'deeparg', 'tsv',
        DEEPARG_AMR.out.versions_deeparg.map { versions -> versions[2] }.first(),
        params.deeparg_db_version, true
    )
    HAMRONIZATION_RESFINDER_AMR(RESFINDER_AMR.out.json, 'resfinder', 'tsv', '', '', false)
    HAMRONIZATION_RGI_AMR(
        RGI_AMR.out.tsv, 'rgi', 'tsv',
        RGI_AMR.out.tool_version.first(), RGI_AMR.out.db_version.first(), true
    )
    ch_hamronized = HAMRONIZATION_ABRICATE_AMR.out.tsv.map { meta, tsv -> tuple(meta, 'abricate', tsv) }
        .mix(HAMRONIZATION_AMRFINDERPLUS_AMR.out.tsv.map { meta, tsv -> tuple(meta, 'amrfinderplus', tsv) })
        .mix(HAMRONIZATION_DEEPARG_AMR.out.tsv.map { meta, tsv -> tuple(meta, 'deeparg', tsv) })
        .mix(HAMRONIZATION_RESFINDER_AMR.out.tsv.map { meta, tsv -> tuple(meta, 'resfinder', tsv) })
        .mix(HAMRONIZATION_RGI_AMR.out.tsv.map { meta, tsv -> tuple(meta, 'rgi', tsv) })
    ch_versions = ch_versions
        .mix(HAMRONIZATION_ABRICATE_AMR.out.versions)
        .mix(HAMRONIZATION_AMRFINDERPLUS_AMR.out.versions)
        .mix(HAMRONIZATION_DEEPARG_AMR.out.versions)
        .mix(HAMRONIZATION_RESFINDER_AMR.out.versions)
        .mix(HAMRONIZATION_RGI_AMR.out.versions)

    HAMRONIZATION_SUMMARIZE_AMR(
        ch_hamronized.map { _meta, _tool, tsv -> tsv }.collect(),
        'tsv'
    )
    ch_versions = ch_versions.mix(HAMRONIZATION_SUMMARIZE_AMR.out.versions)

    // ---- 6_hamr.sh: ARGnorm on hAMRonized output; RGI is not supported ----
    // Header-only tables are skipped: ARGnorm fails on them (pandas EmptyDataError).
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
    ch_versions = ch_versions
        .mix(ARGNORM_SUMMARY_AMR.out.versions)
        .mix(ARGNORM_TOOL_AMR.out.versions)

    emit:
    amrfinder_report = AMRFINDERPLUS_AMR.out.report
    abricate_report = ABRICATE_AMR.out.report
    rgi_tsv = RGI_AMR.out.tsv
    rgi_json = RGI_AMR.out.json
    resfinder_json = RESFINDER_AMR.out.json
    resfinder_results = RESFINDER_AMR.out.results
    deeparg_arg = DEEPARG_AMR.out.arg
    deeparg_potential = DEEPARG_AMR.out.potential_arg
    hamronized = ch_hamronized.map { meta, _tool, tsv -> tuple(meta, tsv) }
    hamronization_summary = HAMRONIZATION_SUMMARIZE_AMR.out.tsv
    normalized = ARGNORM_TOOL_AMR.out.tsv
    normalized_summary = ARGNORM_SUMMARY_AMR.out.tsv
    versions = ch_versions
}
