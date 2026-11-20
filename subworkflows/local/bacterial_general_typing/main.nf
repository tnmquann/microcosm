#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// Authors: @tnmquann (co-author)
// Maintainers: @tnmquann (co-maintainer)

include { GUNZIP as GUNZIP_GENERAL_TYPING } from '../../../modules/nf-core/gunzip/main'
include { MLST as MLST_TYPING } from '../../../modules/nf-core/mlst/main'
include { VIRULENCEFINDER_DOWNLOADDB_LOCAL } from '../../../modules/local/virulencefinder/downloaddb/main'
include { VIRULENCEFINDER_RUN as VIRULENCEFINDER_TYPING } from '../../../modules/local/virulencefinder/main'
include { PLASMIDFINDER as PLASMIDFINDER_TYPING } from '../../../modules/nf-core/plasmidfinder/main'
include { stepWanted } from '../utils_bacterial_lineage/main'

// Steps of playground/modular/2_general_bac.sh for every d__Bacteria genome:
// MLST, VirulenceFinder and PlasmidFinder (container-bundled database).
// Each tool can be turned off with --skip_mlst / --skip_virulencefinder / --skip_plasmidfinder.
workflow BACTERIAL_GENERAL_TYPING {
    take:
    bacterial_genomes_ch

    main:
    ch_versions = channel.empty()

    // --bacterial_decision_table: only genomes with at least one of these steps
    ch_genome_splits = bacterial_genomes_ch
        .filter { meta, _genome -> ['mlst', 'virulencefinder', 'plasmidfinder'].any { step -> stepWanted(meta, step) } }
        .branch { _meta, genome ->
            compressed: genome.getName().endsWith('.gz')
            uncompressed: true
        }
    GUNZIP_GENERAL_TYPING(ch_genome_splits.compressed)
    ch_typing_genomes = GUNZIP_GENERAL_TYPING.out.gunzip
        .mix(ch_genome_splits.uncompressed)
    ch_versions = ch_versions.mix(GUNZIP_GENERAL_TYPING.out.versions_gunzip)

    ch_mlst_tsv = channel.empty()
    ch_vf_json = channel.empty()
    ch_vf_legacy_json = channel.empty()
    ch_vf_tsv = channel.empty()
    ch_vf_txt = channel.empty()
    ch_pf_json = channel.empty()
    ch_pf_txt = channel.empty()
    ch_pf_tsv = channel.empty()

    if (!params.skip_mlst) {
        ch_mlst_genomes = ch_typing_genomes.filter { meta, _genome -> stepWanted(meta, 'mlst') }
        MLST_TYPING(ch_mlst_genomes)
        ch_mlst_tsv = MLST_TYPING.out.tsv
        ch_versions = ch_versions.mix(MLST_TYPING.out.versions_mlst)
    }

    if (!params.skip_virulencefinder) {
        ch_virulencefinder_genomes = ch_typing_genomes.filter { meta, _genome -> stepWanted(meta, 'virulencefinder') }
        if (params.virulencefinder_db) {
            ch_virulencefinder_db = channel.value(file(params.virulencefinder_db, checkIfExists: true))
        } else {
            VIRULENCEFINDER_DOWNLOADDB_LOCAL(
                ch_virulencefinder_genomes.first().map { _genome -> params.virulencefinder_db_revision }
            )
            ch_virulencefinder_db = VIRULENCEFINDER_DOWNLOADDB_LOCAL.out.db
            ch_versions = ch_versions.mix(VIRULENCEFINDER_DOWNLOADDB_LOCAL.out.versions)
        }
        VIRULENCEFINDER_TYPING(ch_virulencefinder_genomes, ch_virulencefinder_db)
        ch_vf_json = VIRULENCEFINDER_TYPING.out.json
        ch_vf_legacy_json = VIRULENCEFINDER_TYPING.out.legacy_json
        ch_vf_tsv = VIRULENCEFINDER_TYPING.out.tsv
        ch_vf_txt = VIRULENCEFINDER_TYPING.out.txt
        ch_versions = ch_versions.mix(VIRULENCEFINDER_TYPING.out.versions)
    }

    if (!params.skip_plasmidfinder) {
        ch_plasmidfinder_genomes = ch_typing_genomes.filter { meta, _genome -> stepWanted(meta, 'plasmidfinder') }
        PLASMIDFINDER_TYPING(ch_plasmidfinder_genomes)
        ch_pf_json = PLASMIDFINDER_TYPING.out.json
        ch_pf_txt = PLASMIDFINDER_TYPING.out.txt
        ch_pf_tsv = PLASMIDFINDER_TYPING.out.tsv
        ch_versions = ch_versions.mix(PLASMIDFINDER_TYPING.out.versions_plasmidfinder)
    }

    emit:
    mlst_tsv = ch_mlst_tsv
    virulencefinder_json = ch_vf_json
    virulencefinder_legacy_json = ch_vf_legacy_json
    virulencefinder_tsv = ch_vf_tsv
    virulencefinder_txt = ch_vf_txt
    plasmidfinder_json = ch_pf_json
    plasmidfinder_txt = ch_pf_txt
    plasmidfinder_tsv = ch_pf_tsv
    versions = ch_versions.unique()
}
