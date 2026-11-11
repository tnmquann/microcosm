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
include { VALIDATE_PARATYPE_BUNDLE } from '../../../modules/local/input/validate_paratype_bundle/main'
include { PARATYPE_RUN } from '../../../modules/local/paratype/run/main'

def keyedSalmonellaOutput(output_channel) {
    return output_channel.map { meta, output ->
        tuple(meta.salmonella_typing_key, meta, output)
    }
}

workflow SALMONELLA_TYPING {
    take:
    bacterial_genomes_ch

    main:
    ch_versions = channel.empty()
    ch_paratype_report = channel.empty()
    pipeline_root = params.pipeline_root ?: projectDir
    parser_script = file(
        "${pipeline_root}/bin/py_scripts/assembly_based/parse_salmonella_typing.py",
        checkIfExists: true
    )

    ch_salmonella_genomes = bacterial_genomes_ch
        .filter { meta, _genome ->
            meta.domain == 'd__Bacteria' && meta.genus == 'g__Salmonella'
        }
        .map { meta, genome ->
            def meta_new = meta.clone()
            meta_new.salmonella_typing_key = "${meta.sample_id}::${meta.genome_id}"
            tuple(meta_new, genome)
        }

    ch_salmonella_splits = ch_salmonella_genomes.branch { _meta, genome ->
        compressed: genome.getName().endsWith('.gz')
        uncompressed: !genome.getName().endsWith('.gz')
    }
    GUNZIP_SALMONELLA_TYPING(ch_salmonella_splits.compressed)
    ch_typing_genomes = GUNZIP_SALMONELLA_TYPING.out.gunzip
        .mix(ch_salmonella_splits.uncompressed)
    ch_versions = ch_versions.mix(GUNZIP_SALMONELLA_TYPING.out.versions_gunzip)

    SEQSERO2_SALMONELLA(ch_typing_genomes)
    SISTR_SALMONELLA(ch_typing_genomes)
    // Mykrobe's current supported Salmonella panel is Typhi.  Its structured
    // coverage evidence is parsed before any Typhi-only follow-up is routed.
    MYKROBE_SALMONELLA(ch_typing_genomes, 'typhi')

    ch_versions = ch_versions
        .mix(SEQSERO2_SALMONELLA.out.versions_seqsero2)
        .mix(SISTR_SALMONELLA.out.versions_sistr)
        .mix(MYKROBE_SALMONELLA.out.versions_mykrobe)

    ch_parser_input = keyedSalmonellaOutput(SEQSERO2_SALMONELLA.out.tsv)
        .join(keyedSalmonellaOutput(MYKROBE_SALMONELLA.out.json), by: 0)
        .join(keyedSalmonellaOutput(SISTR_SALMONELLA.out.tsv), by: 0)
        .map { _key, meta, seqsero2_tsv, _mykrobe_meta, mykrobe_json, _sistr_meta, sistr_tsv ->
            tuple(meta, seqsero2_tsv, mykrobe_json, sistr_tsv)
        }
    PARSE_SALMONELLA_TYPING(ch_parser_input, parser_script)
    ch_versions = ch_versions.mix(PARSE_SALMONELLA_TYPING.out.versions)

    ch_routing_rows = PARSE_SALMONELLA_TYPING.out.routing
        .map { _meta, routing -> routing }
        .splitCsv(header: true, sep: '\t', strip: false)

    ch_genotyphi_keys = ch_routing_rows
        .filter { row -> row.route_genotyphi == 'true' }
        .map { row -> tuple(row.routing_key, true) }
    GENOTYPHI_SALMONELLA(
        keyedSalmonellaOutput(MYKROBE_SALMONELLA.out.json)
            .join(ch_genotyphi_keys, by: 0)
            .map { _key, meta, mykrobe_json, _route -> tuple(meta, mykrobe_json) }
    )
    ch_versions = ch_versions.mix(GENOTYPHI_SALMONELLA.out.versions_genotyphi)

    ch_paratype_keys = ch_routing_rows
        .filter { row -> row.route_paratype == 'true' }
        .map { row -> tuple(row.routing_key, true) }
    ch_paratype_input = keyedSalmonellaOutput(ch_typing_genomes)
        .join(ch_paratype_keys, by: 0)
        .map { _key, meta, genome, _route -> tuple(meta, genome) }

    if (params.paratype_bundle) {
        paratype_bundle_path = file(params.paratype_bundle, checkIfExists: true)
        paratype_validator = file(
            "${pipeline_root}/bin/py_scripts/assembly_based/validate_paratype_bundle.py",
            checkIfExists: true
        )
        VALIDATE_PARATYPE_BUNDLE(paratype_bundle_path, paratype_validator)
        ch_versions = ch_versions.mix(VALIDATE_PARATYPE_BUNDLE.out.versions)
        ch_paratype_ready_input = ch_paratype_input
            .combine(VALIDATE_PARATYPE_BUNDLE.out.ready)
            .map { meta, genome, _ready -> tuple(meta, genome) }
        PARATYPE_RUN(ch_paratype_ready_input, channel.value(paratype_bundle_path))
        ch_paratype_report = PARATYPE_RUN.out.report
        ch_versions = ch_versions.mix(PARATYPE_RUN.out.versions)
    } else {
        // Do not silently omit an explicitly identified Paratyphi A genome.
        // This closure is evaluated only if the evidence-driven branch emits.
        ch_paratype_input.map { meta, _genome ->
            throw new IllegalArgumentException(
                "Paratyphi A evidence was detected for ${meta.sample_id}:${meta.genome_id}; " +
                "supply --paratype_bundle containing paratype.py and its reference/definition files."
            )
        }.view()
    }

    emit:
    seqsero2_tsv = SEQSERO2_SALMONELLA.out.tsv
    seqsero2_txt = SEQSERO2_SALMONELLA.out.txt
    seqsero2_log = SEQSERO2_SALMONELLA.out.log
    sistr_tsv = SISTR_SALMONELLA.out.tsv
    sistr_allele_fasta = SISTR_SALMONELLA.out.allele_fasta
    sistr_allele_json = SISTR_SALMONELLA.out.allele_json
    sistr_cgmlst_csv = SISTR_SALMONELLA.out.cgmlst_csv
    mykrobe_csv = MYKROBE_SALMONELLA.out.csv
    mykrobe_json = MYKROBE_SALMONELLA.out.json
    evidence = PARSE_SALMONELLA_TYPING.out.evidence
    summary = PARSE_SALMONELLA_TYPING.out.evidence
    routing = PARSE_SALMONELLA_TYPING.out.routing
    genotyphi_tsv = GENOTYPHI_SALMONELLA.out.tsv
    paratype_report = ch_paratype_report
    versions = ch_versions.unique()
}
