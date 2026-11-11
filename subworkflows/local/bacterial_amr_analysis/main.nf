#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// Authors: @tnmquann (co-author)
// Maintainers: @tnmquann (co-maintainer)

include { VALIDATE_BACTERIAL_AMR_DATABASES } from '../../../modules/local/input/validate_bacterial_amr_databases/main'
include { STAGE_RGI_DATABASE_ALIASES } from '../../../modules/local/input/stage_rgi_database_aliases/main'
include { PYRODIGAL as PYRODIGAL_AMR } from '../../../modules/nf-core/pyrodigal/main'
include { GUNZIP as GUNZIP_DEEPARG_ORFS } from '../../../modules/nf-core/gunzip/main'
// Local wrappers are intentional here: their software/database/runtime contracts
// differ from the installed nf-core snapshots and are documented in DECISIONS.md.
include { AMRFINDERPLUS_RUN_LOCAL as AMRFINDERPLUS_AMR } from '../../../modules/local/amrfinderplus/run/main'
include { ABRICATE_RUN as ABRICATE_AMR } from '../../../modules/nf-core/abricate/run/main'
include { RGI_MAIN_LOCAL as RGI_AMR } from '../../../modules/local/rgi/main/main'
include { RESFINDER_RUN as RESFINDER_AMR } from '../../../modules/nf-core/resfinder/run/main'
include { DEEPARG_PREDICT as DEEPARG_AMR } from '../../../modules/nf-core/deeparg/predict/main'
// One generic local hAMRonization adapter is reused for supported tool formats
// instead of duplicating four near-identical process wrappers.
include { HAMRONIZATION_REPORT_LOCAL as HAMRONIZATION_AMRFINDERPLUS_AMR } from '../../../modules/local/hamronization/report/main'
include { HAMRONIZATION_REPORT_LOCAL as HAMRONIZATION_ABRICATE_AMR } from '../../../modules/local/hamronization/report/main'
include { HAMRONIZATION_REPORT_LOCAL as HAMRONIZATION_RGI_AMR } from '../../../modules/local/hamronization/report/main'
include { HAMRONIZATION_REPORT_LOCAL as HAMRONIZATION_DEEPARG_AMR } from '../../../modules/local/hamronization/report/main'
include { HAMRONIZATION_SUMMARIZE_LOCAL as HAMRONIZATION_SUMMARIZE_AMR } from '../../../modules/local/hamronization/summarize/main'
include { ARGNORM_LOCAL as ARGNORM_AMRFINDERPLUS_AMR } from '../../../modules/local/argnorm/main'
include { ARGNORM_LOCAL as ARGNORM_ABRICATE_AMR } from '../../../modules/local/argnorm/main'
include { ARGNORM_LOCAL as ARGNORM_DEEPARG_AMR } from '../../../modules/local/argnorm/main'
include { WRITE_AMR_SAMPLE_MANIFEST } from '../../../modules/local/finalize/write_amr_sample_manifest/main'
include {
    WRITE_SKIPPED_AMR_REPORT as WRITE_SKIPPED_AMRFINDERPLUS
    WRITE_SKIPPED_AMR_REPORT as WRITE_SKIPPED_RGI
    WRITE_SKIPPED_AMR_REPORT as WRITE_SKIPPED_RESFINDER
    WRITE_SKIPPED_AMR_REPORT as WRITE_SKIPPED_DEEPARG
} from '../../../modules/local/finalize/write_skipped_amr_report/main'

def unrankedTaxon(value, prefix) {
    if (!value) {
        return ''
    }
    def text = value.toString()
    return text.startsWith(prefix) ? text.substring(prefix.size()) : ''
}

def normalizedTaxon(value) {
    return value
        .toString()
        .toLowerCase()
        .replaceAll(/[^a-z0-9]+/, '_')
        .replaceAll(/^_+|_+$/, '')
}

def supportedTaxon(meta, supportedValues) {
    def supportedByKey = supportedValues.collectEntries { value ->
        [(normalizedTaxon(value)): value]
    }
    def species = unrankedTaxon(meta.species, 's__')
    def genus = unrankedTaxon(meta.genus, 'g__')
    def matched = [species, genus].find { candidate ->
        def key = normalizedTaxon(candidate)
        key && supportedByKey.containsKey(key)
    }
    return matched ? supportedByKey[normalizedTaxon(matched)] : ''
}

def keyedAmrOutput(outputChannel) {
    return outputChannel.map { meta, output ->
        tuple(meta.amr_key, meta, output)
    }
}

workflow BACTERIAL_AMR_ANALYSIS {
    take:
    bacterial_genomes_ch

    main:
    ch_versions = channel.empty()
    ch_amrfinder_report = channel.empty()
    ch_amrfinder_mutation_report = channel.empty()
    ch_rgi_tsv = channel.empty()
    ch_rgi_json = channel.empty()
    ch_resfinder_json = channel.empty()
    ch_resfinder_results = channel.empty()
    ch_pointfinder_results = channel.empty()
    ch_pointfinder_prediction = channel.empty()
    ch_deeparg_arg = channel.empty()
    ch_deeparg_potential = channel.empty()
    ch_hamronization_reports = channel.empty()
    ch_argnorm_outputs = channel.empty()
    pipeline_root = params.pipeline_root ?: projectDir
    amrfinderplus_db_path = params.amrfinderplus_db
        ? file(params.amrfinderplus_db, checkIfExists: true).toString()
        : ''
    abricate_db_path = params.abricate_db
        ? file(params.abricate_db, checkIfExists: true).toString()
        : ''
    rgi_db_path = params.rgi_preparecarddb_dir
        ? file(params.rgi_preparecarddb_dir, checkIfExists: true).toString()
        : ''
    resfinder_db_path = params.resfinder_db
        ? file(params.resfinder_db, checkIfExists: true).toString()
        : ''
    pointfinder_db_path = params.pointfinder_db
        ? file(params.pointfinder_db, checkIfExists: true).toString()
        : ''
    deeparg_db_path = params.deeparg_db
        ? file(params.deeparg_db, checkIfExists: true).toString()
        : ''
    validator_script = file(
        "${pipeline_root}/bin/py_scripts/assembly_based/validate_amr_databases.py",
        checkIfExists: true
    )
    ch_validate_trigger = bacterial_genomes_ch
        .map { _meta, _genome -> true }
        .ifEmpty(false)
        .first()
        .filter { value -> value }

    VALIDATE_BACTERIAL_AMR_DATABASES(
        ch_validate_trigger,
        amrfinderplus_db_path,
        abricate_db_path,
        params.abricate_db_id,
        rgi_db_path,
        params.rgi_include_wildcard,
        resfinder_db_path,
        pointfinder_db_path,
        deeparg_db_path,
        params.deeparg_db_version,
        validator_script
    )
    ch_versions = ch_versions.mix(VALIDATE_BACTERIAL_AMR_DATABASES.out.versions)

    ch_amrfinder_supported = VALIDATE_BACTERIAL_AMR_DATABASES.out.amrfinder_organisms
        .splitCsv(header: true, sep: '\t', strip: false)
        .map { row -> row.organism }
        .collect()
        .ifEmpty([])
        .map { values -> [values: values] }

    ch_pointfinder_supported = VALIDATE_BACTERIAL_AMR_DATABASES.out.pointfinder_species
        .splitCsv(header: true, sep: '\t', strip: false)
        .map { row -> row.database_name }
        .collect()
        .ifEmpty([])
        .map { values -> [values: values] }

    ch_amr_genomes = bacterial_genomes_ch
        .combine(ch_amrfinder_supported)
        .combine(ch_pointfinder_supported)
        .map { meta, genome, amrfinderSupport, pointfinderSupport ->
            def meta_new = meta.clone()
            meta_new.amr_key = "${meta.sample_id}::${meta.genome_id}"
            def amrfinder_organism = supportedTaxon(meta_new, amrfinderSupport.values)
            if (amrfinder_organism) {
                meta_new.organism = amrfinder_organism
            }
            meta_new.pointfinder_species = supportedTaxon(
                meta_new,
                pointfinderSupport.values
            )
            tuple(meta_new, genome)
        }

    has_amrfinder_db = params.amrfinderplus_db?.toString()?.trim()
    has_rgi_db = params.rgi_preparecarddb_dir?.toString()?.trim()
    has_resfinder_db = params.resfinder_db?.toString()?.trim()
    has_deeparg_db = params.deeparg_db?.toString()?.trim()

    ch_abricate_db = VALIDATE_BACTERIAL_AMR_DATABASES.out.ready
        .map { _ready -> params.abricate_db ? file(params.abricate_db, checkIfExists: true) : [] }

    if (has_amrfinder_db) {
        ch_amrfinder_db = VALIDATE_BACTERIAL_AMR_DATABASES.out.ready
            .map { _ready -> file(params.amrfinderplus_db, checkIfExists: true) }
        AMRFINDERPLUS_AMR(ch_amr_genomes, ch_amrfinder_db)
        ch_amrfinder_report = AMRFINDERPLUS_AMR.out.report
        ch_amrfinder_mutation_report = AMRFINDERPLUS_AMR.out.mutation_report
        ch_versions = ch_versions.mix(AMRFINDERPLUS_AMR.out.versions)
    } else {
        WRITE_SKIPPED_AMRFINDERPLUS(
            ch_amr_genomes,
            'amrfinderplus',
            'tsv',
            'no --amrfinderplus_db supplied and no automatic default setup is configured'
        )
        ch_amrfinder_report = WRITE_SKIPPED_AMRFINDERPLUS.out.report
        ch_versions = ch_versions.mix(WRITE_SKIPPED_AMRFINDERPLUS.out.versions)
    }

    ABRICATE_AMR(ch_amr_genomes, ch_abricate_db)
    if (has_rgi_db) {
        STAGE_RGI_DATABASE_ALIASES(
            params.rgi_preparecarddb_dir.toString().trim(),
            VALIDATE_BACTERIAL_AMR_DATABASES.out.ready
        )
        ch_versions = ch_versions.mix(STAGE_RGI_DATABASE_ALIASES.out.versions)
        ch_rgi_wildcard = params.rgi_include_wildcard
            ? STAGE_RGI_DATABASE_ALIASES.out.wildcard
            : channel.value([])
        RGI_AMR(
            ch_amr_genomes,
            STAGE_RGI_DATABASE_ALIASES.out.card,
            ch_rgi_wildcard
        )
        ch_rgi_tsv = RGI_AMR.out.tsv
        ch_rgi_json = RGI_AMR.out.json
        ch_versions = ch_versions.mix(RGI_AMR.out.versions)
    } else {
        WRITE_SKIPPED_RGI(
            ch_amr_genomes,
            'rgi',
            'txt',
            'no --rgi_preparecarddb_dir supplied and no automatic CARD setup is configured'
        )
        ch_rgi_tsv = WRITE_SKIPPED_RGI.out.report
        ch_versions = ch_versions.mix(WRITE_SKIPPED_RGI.out.versions)
    }

    if (has_resfinder_db) {
        ch_resfinder_db = VALIDATE_BACTERIAL_AMR_DATABASES.out.ready
            .map { _ready -> file(params.resfinder_db, checkIfExists: true) }
        ch_pointfinder_db = VALIDATE_BACTERIAL_AMR_DATABASES.out.ready
            .map { _ready -> params.pointfinder_db ? file(params.pointfinder_db, checkIfExists: true) : [] }
        ch_resfinder_input = ch_amr_genomes.map { meta, genome ->
            tuple(meta, [], genome)
        }
        RESFINDER_AMR(ch_resfinder_input, ch_pointfinder_db, ch_resfinder_db)
        ch_resfinder_json = RESFINDER_AMR.out.json
        ch_resfinder_results = RESFINDER_AMR.out.resfinder_results_tab
        ch_pointfinder_results = RESFINDER_AMR.out.pointfinder_results
        ch_pointfinder_prediction = RESFINDER_AMR.out.pointfinder_prediction
        ch_versions = ch_versions.mix(RESFINDER_AMR.out.versions)
    } else {
        WRITE_SKIPPED_RESFINDER(
            ch_amr_genomes,
            'resfinder',
            'json',
            'no --resfinder_db supplied and no automatic ResFinder database setup is configured'
        )
        ch_resfinder_json = WRITE_SKIPPED_RESFINDER.out.report
        ch_versions = ch_versions.mix(WRITE_SKIPPED_RESFINDER.out.versions)
    }

    if (has_deeparg_db) {
        ch_deeparg_db = VALIDATE_BACTERIAL_AMR_DATABASES.out.ready
            .map { _ready -> file(params.deeparg_db, checkIfExists: true) }
        PYRODIGAL_AMR(ch_amr_genomes, 'gff')
        GUNZIP_DEEPARG_ORFS(PYRODIGAL_AMR.out.fna)
        ch_deeparg_input = GUNZIP_DEEPARG_ORFS.out.gunzip.map { meta, orfs ->
            tuple(meta, orfs, 'LS')
        }
        DEEPARG_AMR(ch_deeparg_input, ch_deeparg_db)
        ch_deeparg_arg = DEEPARG_AMR.out.arg
        ch_deeparg_potential = DEEPARG_AMR.out.potential_arg
        ch_versions = ch_versions
            .mix(PYRODIGAL_AMR.out.versions_pyrodigal)
            .mix(GUNZIP_DEEPARG_ORFS.out.versions_gunzip)
            .mix(DEEPARG_AMR.out.versions_deeparg)
    } else {
        WRITE_SKIPPED_DEEPARG(
            ch_amr_genomes,
            'deeparg',
            'ARG',
            'no --deeparg_db supplied and automatic database downloads are disabled'
        )
        ch_deeparg_arg = WRITE_SKIPPED_DEEPARG.out.report
        ch_versions = ch_versions.mix(WRITE_SKIPPED_DEEPARG.out.versions)
    }

    ch_versions = ch_versions.mix(ABRICATE_AMR.out.versions_abricate)

    if (has_amrfinder_db) {
        HAMRONIZATION_AMRFINDERPLUS_AMR(
            ch_amrfinder_report,
            'amrfinderplus',
            'tsv',
            AMRFINDERPLUS_AMR.out.tool_version,
            AMRFINDERPLUS_AMR.out.db_version,
            true
        )
        ARGNORM_AMRFINDERPLUS_AMR(
            ch_amrfinder_report,
            'amrfinderplus',
            'ncbi'
        )
        ch_hamronization_reports = ch_hamronization_reports.mix(HAMRONIZATION_AMRFINDERPLUS_AMR.out.tsv)
        ch_argnorm_outputs = ch_argnorm_outputs.mix(ARGNORM_AMRFINDERPLUS_AMR.out.tsv)
        ch_versions = ch_versions
            .mix(HAMRONIZATION_AMRFINDERPLUS_AMR.out.versions)
            .mix(ARGNORM_AMRFINDERPLUS_AMR.out.versions)
    }

    HAMRONIZATION_ABRICATE_AMR(
        ABRICATE_AMR.out.report,
        'abricate',
        'tsv',
        '1.0.1',
        params.abricate_db_version,
        false
    )
    ch_hamronization_reports = ch_hamronization_reports.mix(HAMRONIZATION_ABRICATE_AMR.out.tsv)

    if (has_rgi_db) {
        HAMRONIZATION_RGI_AMR(
            ch_rgi_tsv,
            'rgi',
            'tsv',
            RGI_AMR.out.tool_version,
            RGI_AMR.out.db_version,
            true
        )
        ch_hamronization_reports = ch_hamronization_reports.mix(HAMRONIZATION_RGI_AMR.out.tsv)
        ch_versions = ch_versions.mix(HAMRONIZATION_RGI_AMR.out.versions)
    }

    if (has_deeparg_db) {
        ch_deeparg_ham_input = ch_deeparg_arg
            .map { meta, report ->
                def meta_new = meta.clone()
                meta_new.amr_result_id = "${meta.genome_id}.deeparg_arg"
                tuple(meta_new, report)
            }
            .mix(
                ch_deeparg_potential.map { meta, report ->
                    def meta_new = meta.clone()
                    meta_new.amr_result_id = "${meta.genome_id}.deeparg_potential"
                    tuple(meta_new, report)
                }
            )
        HAMRONIZATION_DEEPARG_AMR(
            ch_deeparg_ham_input,
            'deeparg',
            'tsv',
            '1.0.4',
            params.deeparg_db_version,
            true
        )
        ARGNORM_DEEPARG_AMR(
            ch_deeparg_ham_input,
            'deeparg',
            'deeparg'
        )
        ch_hamronization_reports = ch_hamronization_reports.mix(
            HAMRONIZATION_DEEPARG_AMR.out.tsv
        )
        ch_argnorm_outputs = ch_argnorm_outputs.mix(ARGNORM_DEEPARG_AMR.out.tsv)
        ch_versions = ch_versions
            .mix(HAMRONIZATION_DEEPARG_AMR.out.versions)
            .mix(ARGNORM_DEEPARG_AMR.out.versions)
    }

    ch_versions = ch_versions.mix(HAMRONIZATION_ABRICATE_AMR.out.versions)

    ch_argnorm_abricate = channel.empty()
    abricate_normalizable = params.abricate_db_id in [
        'ncbi',
        'resfinder',
        'argannot',
        'megares'
    ]
    if (abricate_normalizable) {
        ARGNORM_ABRICATE_AMR(
            ABRICATE_AMR.out.report,
            'abricate',
            params.abricate_db_id
        )
        ch_argnorm_abricate = ARGNORM_ABRICATE_AMR.out.tsv
        ch_versions = ch_versions.mix(ARGNORM_ABRICATE_AMR.out.versions)
    }

    ch_hamronization_report_files = ch_hamronization_reports
        .map { _meta, report -> report }
        .collect()
        .filter { reports -> reports && reports.size() > 0 }
    HAMRONIZATION_SUMMARIZE_AMR(ch_hamronization_report_files, 'tsv')
    ch_versions = ch_versions.mix(
        HAMRONIZATION_SUMMARIZE_AMR.out.versions
    )

    ch_argnorm_outputs = ch_argnorm_outputs
        .mix(ch_argnorm_abricate)

    ch_manifest_input = keyedAmrOutput(ch_amrfinder_report)
        .join(keyedAmrOutput(ABRICATE_AMR.out.report), by: 0)
        .map { key, meta, amrfinder, _abricateMeta, abricate ->
            tuple(key, meta, amrfinder, abricate)
        }
        .join(keyedAmrOutput(ch_rgi_tsv), by: 0)
        .map { key, meta, amrfinder, abricate, _rgiMeta, rgi ->
            tuple(key, meta, amrfinder, abricate, rgi)
        }
        .join(keyedAmrOutput(ch_resfinder_json), by: 0)
        .map { key, meta, amrfinder, abricate, rgi, _resfinderMeta, resfinder ->
            tuple(key, meta, amrfinder, abricate, rgi, resfinder)
        }
        .join(keyedAmrOutput(ch_deeparg_arg), by: 0)
        .map {
            _key,
            meta,
            amrfinder,
            abricate,
            rgi,
            resfinder,
            _deepargMeta,
            deeparg ->
                tuple(meta, amrfinder, abricate, rgi, resfinder, deeparg)
        }

    manifest_script = file(
        "${pipeline_root}/bin/py_scripts/assembly_based/write_amr_manifest.py",
        checkIfExists: true
    )
    WRITE_AMR_SAMPLE_MANIFEST(
        ch_manifest_input,
        abricate_normalizable,
        manifest_script,
        VALIDATE_BACTERIAL_AMR_DATABASES.out.validation.first()
    )
    ch_versions = ch_versions.mix(WRITE_AMR_SAMPLE_MANIFEST.out.versions)

    emit:
    database_validation = VALIDATE_BACTERIAL_AMR_DATABASES.out.validation
    amrfinder_report = ch_amrfinder_report
    amrfinder_mutations = ch_amrfinder_mutation_report
    abricate_report = ABRICATE_AMR.out.report
    rgi_tsv = ch_rgi_tsv
    rgi_json = ch_rgi_json
    resfinder_json = ch_resfinder_json
    resfinder_results = ch_resfinder_results
    pointfinder_results = ch_pointfinder_results
    pointfinder_prediction = ch_pointfinder_prediction
    deeparg_arg = ch_deeparg_arg
    deeparg_potential = ch_deeparg_potential
    hamronized = ch_hamronization_reports
    hamronization_summary = HAMRONIZATION_SUMMARIZE_AMR.out.tsv
    normalized = ch_argnorm_outputs
    manifests = WRITE_AMR_SAMPLE_MANIFEST.out.manifest
    versions = ch_versions.unique()
}
