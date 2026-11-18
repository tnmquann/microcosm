#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// Authors: @tnmquann (co-author)
// Maintainers: @tnmquann (co-maintainer)

include { GUNZIP as GUNZIP_SALMONELLA_TYPING } from '../../../modules/nf-core/gunzip/main'
include { SEQSERO2 as SEQSERO2_SALMONELLA } from '../../../modules/nf-core/seqsero2/main'
include { SISTR as SISTR_SALMONELLA } from '../../../modules/nf-core/sistr/main'
include { MYKROBE_PREDICT as MYKROBE_SALMONELLA } from '../../../modules/nf-core/mykrobe/predict/main'
include { GENOTYPHI_PARSE as GENOTYPHI_SALMONELLA } from '../../../modules/nf-core/genotyphi/parse/main'
include { PARSE_SALMONELLA_TYPING } from '../../../modules/local/finalize/parse_salmonella_typing/main'
include { PARATYPE_DOWNLOAD_LOCAL } from '../../../modules/local/paratype/download/main'
include { VALIDATE_PARATYPE_BUNDLE } from '../../../modules/local/input/validate_paratype_bundle/main'
include { PARATYPE_RUN } from '../../../modules/local/paratype/run/main'
include { stripGtdbSuffixes } from '../utils_bacterial_lineage/main'

def keyedSalmonellaOutput(output_channel) {
    return output_channel.map { meta, output ->
        tuple(meta.salmonella_typing_key, meta, output)
    }
}

// Flow of playground/modular/4_sal.sh: SISTR + SeqSero2 on every Salmonella
// genome; the SISTR serovar alone then routes Typhi -> Mykrobe -> GenoTyphi
// and Paratyphi A -> Paratype. Each tool has a --skip_<tool> switch.
workflow SALMONELLA_TYPING {
    take:
    bacterial_genomes_ch

    main:
    ch_versions = channel.empty()
    pipeline_root = params.pipeline_root ?: projectDir
    parser_script = file(
        "${pipeline_root}/bin/py_scripts/bacterial_analysis/parse_salmonella_typing.py",
        checkIfExists: true
    )

    ch_salmonella_genomes = bacterial_genomes_ch
        .filter { meta, _genome -> stripGtdbSuffixes(meta.lineage).contains('g__Salmonella') }
        .map { meta, genome ->
            def meta_new = meta.clone()
            meta_new.salmonella_typing_key = "${meta.sample_id}::${meta.genome_id}"
            tuple(meta_new, genome)
        }

    ch_salmonella_splits = ch_salmonella_genomes.branch { _meta, genome ->
        compressed: genome.getName().endsWith('.gz')
        uncompressed: true
    }
    GUNZIP_SALMONELLA_TYPING(ch_salmonella_splits.compressed)
    ch_typing_genomes = GUNZIP_SALMONELLA_TYPING.out.gunzip
        .mix(ch_salmonella_splits.uncompressed)
    ch_versions = ch_versions.mix(GUNZIP_SALMONELLA_TYPING.out.versions_gunzip)

    ch_sistr_tsv = channel.empty()
    ch_seqsero2_tsv = channel.empty()
    ch_evidence = channel.empty()
    ch_routing = channel.empty()
    ch_mykrobe_json = channel.empty()
    ch_genotyphi_tsv = channel.empty()
    ch_paratype_report = channel.empty()

    if (!params.skip_sistr) {
        SISTR_SALMONELLA(ch_typing_genomes)
        ch_sistr_tsv = SISTR_SALMONELLA.out.tsv
        ch_versions = ch_versions.mix(SISTR_SALMONELLA.out.versions_sistr)
    }
    if (!params.skip_seqsero2) {
        SEQSERO2_SALMONELLA(ch_typing_genomes)
        ch_seqsero2_tsv = SEQSERO2_SALMONELLA.out.tsv
        ch_versions = ch_versions.mix(SEQSERO2_SALMONELLA.out.versions_seqsero2)
    }

    // The SISTR serovar is the only routing input; without SISTR there is no follow-up.
    if (!params.skip_sistr) {
        // Route only after both tools finished for the same genome.
        ch_parser_input = params.skip_seqsero2
            ? ch_sistr_tsv.map { meta, sistr_tsv -> tuple(meta, [], sistr_tsv) }
            : keyedSalmonellaOutput(ch_seqsero2_tsv)
                .join(keyedSalmonellaOutput(ch_sistr_tsv), by: 0)
                .map { _key, meta, seqsero2_tsv, _sistr_meta, sistr_tsv ->
                    tuple(meta, seqsero2_tsv, sistr_tsv)
                }
        PARSE_SALMONELLA_TYPING(ch_parser_input, parser_script)
        ch_evidence = PARSE_SALMONELLA_TYPING.out.evidence
        ch_routing = PARSE_SALMONELLA_TYPING.out.routing
        ch_versions = ch_versions.mix(PARSE_SALMONELLA_TYPING.out.versions)

        ch_routes = ch_routing
            .map { _meta, routing -> routing }
            .splitCsv(header: true, sep: '\t', strip: false)

        ch_routed_genomes = keyedSalmonellaOutput(ch_typing_genomes)
            .join(ch_routes.map { row -> tuple(row.routing_key, row) }, by: 0)
            .branch { _key, _meta, _genome, row ->
                typhi: row.route_genotyphi == 'true'
                paratyphi_a: row.route_paratype == 'true'
            }

        // Typhi: Mykrobe (typhi panel) then the GenoTyphi parser.
        if (!params.skip_genotyphi) {
            MYKROBE_SALMONELLA(
                ch_routed_genomes.typhi.map { _key, meta, genome, _row -> tuple(meta, genome) },
                'typhi'
            )
            GENOTYPHI_SALMONELLA(MYKROBE_SALMONELLA.out.json)
            ch_mykrobe_json = MYKROBE_SALMONELLA.out.json
            ch_genotyphi_tsv = GENOTYPHI_SALMONELLA.out.tsv
            ch_versions = ch_versions
                .mix(MYKROBE_SALMONELLA.out.versions_mykrobe)
                .mix(GENOTYPHI_SALMONELLA.out.versions_genotyphi)
        }

        // Paratyphi A: Paratype; fetch it only when a genome needs it.
        if (!params.skip_paratype) {
            ch_paratype_input = ch_routed_genomes.paratyphi_a
                .map { _key, meta, genome, _row -> tuple(meta, genome) }
            if (params.paratype_bundle) {
                ch_paratype_bundle = channel.value(file(params.paratype_bundle, checkIfExists: true))
            } else {
                PARATYPE_DOWNLOAD_LOCAL(ch_paratype_input.first().map { _input -> params.paratype_revision })
                ch_paratype_bundle = PARATYPE_DOWNLOAD_LOCAL.out.bundle
                ch_versions = ch_versions.mix(PARATYPE_DOWNLOAD_LOCAL.out.versions)
            }
            VALIDATE_PARATYPE_BUNDLE(
                ch_paratype_bundle,
                file("${pipeline_root}/bin/py_scripts/bacterial_analysis/validate_paratype_bundle.py", checkIfExists: true)
            )
            PARATYPE_RUN(
                ch_paratype_input
                    .combine(VALIDATE_PARATYPE_BUNDLE.out.ready)
                    .map { meta, genome, _ready -> tuple(meta, genome) },
                ch_paratype_bundle
            )
            ch_paratype_report = PARATYPE_RUN.out.report
            ch_versions = ch_versions
                .mix(VALIDATE_PARATYPE_BUNDLE.out.versions)
                .mix(PARATYPE_RUN.out.versions)
        }
    }

    emit:
    seqsero2_tsv = ch_seqsero2_tsv
    sistr_tsv = ch_sistr_tsv
    mykrobe_json = ch_mykrobe_json
    evidence = ch_evidence
    summary = ch_evidence
    routing = ch_routing
    genotyphi_tsv = ch_genotyphi_tsv
    paratype_report = ch_paratype_report
    versions = ch_versions.unique()
}
