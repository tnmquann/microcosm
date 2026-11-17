#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// Authors: @tnmquann (co-author)
// Maintainers: @tnmquann (co-maintainer)

include { GUNZIP as GUNZIP_GENERAL_TYPING } from '../../../modules/nf-core/gunzip/main'
include { MLST as MLST_TYPING } from '../../../modules/nf-core/mlst/main'
include { VIRULENCEFINDER_DOWNLOADDB_LOCAL } from '../../../modules/local/virulencefinder/downloaddb/main'
include { VIRULENCEFINDER_RUN as VIRULENCEFINDER_TYPING } from '../../../modules/local/virulencefinder/main'
include { PLASMIDFINDER as PLASMIDFINDER_TYPING } from '../../../modules/nf-core/plasmidfinder/main'

// virulencefinder_db revision downloaded when --virulencefinder_db is not set
// (tag virulencefinder-3.2.1, matching the VirulenceFinder 3.2.1 container).
def virulencefinderDbRevision() {
    return '065b8d1'
}

// Steps of playground/modular/2_general_bac.sh for every d__Bacteria genome:
// MLST, VirulenceFinder and PlasmidFinder (container-bundled database).
workflow BACTERIAL_GENERAL_TYPING {
    take:
    bacterial_genomes_ch

    main:
    ch_versions = channel.empty()

    ch_genome_splits = bacterial_genomes_ch.branch { _meta, genome ->
        compressed: genome.getName().endsWith('.gz')
        uncompressed: true
    }
    GUNZIP_GENERAL_TYPING(ch_genome_splits.compressed)
    ch_typing_genomes = GUNZIP_GENERAL_TYPING.out.gunzip
        .mix(ch_genome_splits.uncompressed)
    ch_versions = ch_versions.mix(GUNZIP_GENERAL_TYPING.out.versions_gunzip)

    if (params.virulencefinder_db) {
        ch_virulencefinder_db = channel.value(file(params.virulencefinder_db, checkIfExists: true))
    } else {
        VIRULENCEFINDER_DOWNLOADDB_LOCAL(
            ch_typing_genomes.first().map { _genome -> virulencefinderDbRevision() }
        )
        ch_virulencefinder_db = VIRULENCEFINDER_DOWNLOADDB_LOCAL.out.db
        ch_versions = ch_versions.mix(VIRULENCEFINDER_DOWNLOADDB_LOCAL.out.versions)
    }

    MLST_TYPING(ch_typing_genomes)
    VIRULENCEFINDER_TYPING(ch_typing_genomes, ch_virulencefinder_db)
    PLASMIDFINDER_TYPING(ch_typing_genomes)
    ch_versions = ch_versions
        .mix(MLST_TYPING.out.versions_mlst)
        .mix(VIRULENCEFINDER_TYPING.out.versions)
        .mix(PLASMIDFINDER_TYPING.out.versions_plasmidfinder)

    emit:
    mlst_tsv = MLST_TYPING.out.tsv
    virulencefinder_json = VIRULENCEFINDER_TYPING.out.json
    virulencefinder_legacy_json = VIRULENCEFINDER_TYPING.out.legacy_json
    virulencefinder_tsv = VIRULENCEFINDER_TYPING.out.tsv
    virulencefinder_txt = VIRULENCEFINDER_TYPING.out.txt
    plasmidfinder_json = PLASMIDFINDER_TYPING.out.json
    plasmidfinder_txt = PLASMIDFINDER_TYPING.out.txt
    plasmidfinder_tsv = PLASMIDFINDER_TYPING.out.tsv
    versions = ch_versions.unique()
}
