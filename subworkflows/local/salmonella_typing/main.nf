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

// Paratype commit used when --paratype_bundle is not supplied.
def paratypeRevision() {
    return '099db9b4de828e9c993fe76164f85731f6368e86'
}

def keyedSalmonellaOutput(output_channel) {
    return output_channel.map { meta, output ->
        tuple(meta.salmonella_typing_key, meta, output)
    }
}

// Flow of playground/modular/4_sal.sh: SISTR + SeqSero2 on every Salmonella
// genome; the SISTR serovar alone then routes Typhi -> Mykrobe -> GenoTyphi
// and Paratyphi A -> Paratype.
workflow SALMONELLA_TYPING {
    take:
    bacterial_genomes_ch

    main:
    ch_versions = channel.empty()
    pipeline_root = params.pipeline_root ?: projectDir
    parser_script = file(
        "${pipeline_root}/bin/py_scripts/assembly_based/parse_salmonella_typing.py",
        checkIfExists: true
    )

    ch_salmonella_genomes = bacterial_genomes_ch
        .filter { meta, _genome -> (meta.lineage ?: '').contains('g__Salmonella') }
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

    SISTR_SALMONELLA(ch_typing_genomes)
    SEQSERO2_SALMONELLA(ch_typing_genomes)
    ch_versions = ch_versions
        .mix(SEQSERO2_SALMONELLA.out.versions_seqsero2)
        .mix(SISTR_SALMONELLA.out.versions_sistr)

    // Route only after both tools finished for the same genome.
    ch_parser_input = keyedSalmonellaOutput(SEQSERO2_SALMONELLA.out.tsv)
        .join(keyedSalmonellaOutput(SISTR_SALMONELLA.out.tsv), by: 0)
        .map { _key, meta, seqsero2_tsv, _sistr_meta, sistr_tsv ->
            tuple(meta, seqsero2_tsv, sistr_tsv)
        }
    PARSE_SALMONELLA_TYPING(ch_parser_input, parser_script)
    ch_versions = ch_versions.mix(PARSE_SALMONELLA_TYPING.out.versions)

    ch_routes = PARSE_SALMONELLA_TYPING.out.routing
        .map { _meta, routing -> routing }
        .splitCsv(header: true, sep: '\t', strip: false)

    ch_routed_genomes = keyedSalmonellaOutput(ch_typing_genomes)
        .join(ch_routes.map { row -> tuple(row.routing_key, row) }, by: 0)
        .branch { _key, _meta, _genome, row ->
            typhi: row.route_genotyphi == 'true'
            paratyphi_a: row.route_paratype == 'true'
        }

    // Typhi: Mykrobe (typhi panel) then the GenoTyphi parser.
    MYKROBE_SALMONELLA(
        ch_routed_genomes.typhi.map { _key, meta, genome, _row -> tuple(meta, genome) },
        'typhi'
    )
    GENOTYPHI_SALMONELLA(MYKROBE_SALMONELLA.out.json)
    ch_versions = ch_versions
        .mix(MYKROBE_SALMONELLA.out.versions_mykrobe)
        .mix(GENOTYPHI_SALMONELLA.out.versions_genotyphi)

    // Paratyphi A: Paratype; fetch it only when a genome needs it.
    ch_paratype_input = ch_routed_genomes.paratyphi_a
        .map { _key, meta, genome, _row -> tuple(meta, genome) }
    if (params.paratype_bundle) {
        ch_paratype_bundle = channel.value(file(params.paratype_bundle, checkIfExists: true))
    } else {
        PARATYPE_DOWNLOAD_LOCAL(ch_paratype_input.first().map { _input -> paratypeRevision() })
        ch_paratype_bundle = PARATYPE_DOWNLOAD_LOCAL.out.bundle
        ch_versions = ch_versions.mix(PARATYPE_DOWNLOAD_LOCAL.out.versions)
    }
    VALIDATE_PARATYPE_BUNDLE(
        ch_paratype_bundle,
        file("${pipeline_root}/bin/py_scripts/assembly_based/validate_paratype_bundle.py", checkIfExists: true)
    )
    PARATYPE_RUN(
        ch_paratype_input
            .combine(VALIDATE_PARATYPE_BUNDLE.out.ready)
            .map { meta, genome, _ready -> tuple(meta, genome) },
        ch_paratype_bundle
    )
    ch_versions = ch_versions
        .mix(VALIDATE_PARATYPE_BUNDLE.out.versions)
        .mix(PARATYPE_RUN.out.versions)

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
    paratype_report = PARATYPE_RUN.out.report
    versions = ch_versions.unique()
}
