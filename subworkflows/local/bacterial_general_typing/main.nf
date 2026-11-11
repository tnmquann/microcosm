#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// Authors: @tnmquann (co-author)
// Maintainers: @tnmquann (co-maintainer)

include { VALIDATE_BACTERIAL_TYPING_DATABASES } from '../../../modules/local/input/validate_bacterial_typing_databases/main'
include { GUNZIP as GUNZIP_GENERAL_TYPING } from '../../../modules/nf-core/gunzip/main'
include { VIRULENCEFINDER_RUN as VIRULENCEFINDER_TYPING } from '../../../modules/local/virulencefinder/main'
include { PLASMIDFINDER as PLASMIDFINDER_TYPING } from '../../../modules/nf-core/plasmidfinder/main'
include { MLST as MLST_TYPING } from '../../../modules/nf-core/mlst/main'
include { BLAST_MAKEBLASTDB as BLAST_MAKEBLASTDB_TYPING } from '../../../modules/nf-core/blast/makeblastdb/main'
include { BLAST_BLASTN as BLAST_BLASTN_TYPING } from '../../../modules/nf-core/blast/blastn/main'
include { FILTER_BACTERIAL_BLAST_HITS } from '../../../modules/local/finalize/filter_bacterial_blast_hits/main'
include { WRITE_GENERAL_TYPING_MANIFEST } from '../../../modules/local/finalize/write_general_typing_manifest/main'
include { SUMMARIZE_GENERAL_TYPING } from '../../../modules/local/finalize/summarize_general_typing/main'
include { WRITE_SKIPPED_BLASTN } from '../../../modules/local/finalize/write_skipped_blastn/main'
include { WRITE_SKIPPED_VIRULENCEFINDER } from '../../../modules/local/finalize/write_skipped_virulencefinder/main'

def keyedTypingOutput(outputChannel) {
    return outputChannel.map { meta, output ->
        tuple(meta.typing_key, meta, output)
    }
}

workflow BACTERIAL_GENERAL_TYPING {
    take:
    bacterial_genomes_ch

    main:
    ch_versions = channel.empty()
    ch_virulencefinder_json = channel.empty()
    ch_virulencefinder_legacy_json = channel.empty()
    ch_virulencefinder_tsv = channel.empty()
    ch_virulencefinder_txt = channel.empty()
    ch_virulencefinder_genome_seq = channel.empty()
    ch_virulencefinder_reference_seq = channel.empty()
    ch_blastn_output = channel.empty()
    ch_blastn_filtered = channel.empty()
    ch_blastn_summary = channel.empty()
    ch_blastn_manifest_output = channel.empty()
    pipeline_root = params.pipeline_root ?: projectDir
    virulencefinder_db_path = params.virulencefinder_db
        ? file(params.virulencefinder_db, checkIfExists: true).toString()
        : ''
    blast_reference_path = params.bacterial_blast_reference
        ? file(params.bacterial_blast_reference, checkIfExists: true).toString()
        : ''
    blast_database_path = params.bacterial_blast_db
        ? file(params.bacterial_blast_db, checkIfExists: true).toString()
        : ''
    blast_database_label = params.bacterial_blast_db_label?.toString()?.trim() ?: ''
    blast_min_identity = params.bacterial_blast_min_identity != null
        ? params.bacterial_blast_min_identity
        : ''
    blast_min_coverage = params.bacterial_blast_min_coverage != null
        ? params.bacterial_blast_min_coverage
        : ''
    validator_script = file(
        "${pipeline_root}/bin/py_scripts/assembly_based/validate_general_typing_databases.py",
        checkIfExists: true
    )
    validator_options = tuple(
        blast_database_label,
        blast_min_identity,
        blast_min_coverage,
        validator_script
    )
    manifest_script = file(
        "${pipeline_root}/bin/py_scripts/assembly_based/write_general_typing_manifest.py",
        checkIfExists: true
    )
    blast_filter_script = file(
        "${pipeline_root}/bin/py_scripts/assembly_based/filter_blasthits.py",
        checkIfExists: true
    )
    summary_script = file(
        "${pipeline_root}/bin/py_scripts/assembly_based/summarize_gene_hits.py",
        checkIfExists: true
    )

    ch_validate_trigger = bacterial_genomes_ch
        .map { _meta, _genome -> true }
        .ifEmpty(false)
        .first()
        .filter { value -> value }

    VALIDATE_BACTERIAL_TYPING_DATABASES(
        ch_validate_trigger,
        virulencefinder_db_path,
        blast_reference_path,
        blast_database_path,
        validator_options
    )
    ch_versions = ch_versions.mix(VALIDATE_BACTERIAL_TYPING_DATABASES.out.versions)

    ch_validated_genomes = bacterial_genomes_ch
        .combine(VALIDATE_BACTERIAL_TYPING_DATABASES.out.ready)
        .map { meta, genome, _ready ->
            def meta_new = meta.clone()
            meta_new.typing_key = "${meta.sample_id}::${meta.genome_id}"
            tuple(meta_new, genome)
        }

    ch_genome_splits = ch_validated_genomes.branch { _meta, genome ->
        compressed: genome.getName().endsWith('.gz')
        uncompressed: !genome.getName().endsWith('.gz')
    }
    GUNZIP_GENERAL_TYPING(ch_genome_splits.compressed)
    ch_typing_genomes = GUNZIP_GENERAL_TYPING.out.gunzip
        .mix(ch_genome_splits.uncompressed)
    ch_versions = ch_versions.mix(GUNZIP_GENERAL_TYPING.out.versions_gunzip)

    has_virulencefinder_db = params.virulencefinder_db?.toString()?.trim()
    if (has_virulencefinder_db) {
        ch_virulencefinder_db = VALIDATE_BACTERIAL_TYPING_DATABASES.out.ready
            .map { _ready -> file(params.virulencefinder_db, checkIfExists: true) }
        VIRULENCEFINDER_TYPING(ch_typing_genomes, ch_virulencefinder_db)
        ch_virulencefinder_json = VIRULENCEFINDER_TYPING.out.json
        ch_virulencefinder_legacy_json = VIRULENCEFINDER_TYPING.out.legacy_json
        ch_virulencefinder_tsv = VIRULENCEFINDER_TYPING.out.tsv
        ch_virulencefinder_txt = VIRULENCEFINDER_TYPING.out.txt
        ch_virulencefinder_genome_seq = VIRULENCEFINDER_TYPING.out.genome_seq
        ch_virulencefinder_reference_seq = VIRULENCEFINDER_TYPING.out.virulence_seq
        ch_versions = ch_versions.mix(VIRULENCEFINDER_TYPING.out.versions)
    } else {
        WRITE_SKIPPED_VIRULENCEFINDER(ch_typing_genomes)
        ch_virulencefinder_json = WRITE_SKIPPED_VIRULENCEFINDER.out.json
        ch_virulencefinder_legacy_json = WRITE_SKIPPED_VIRULENCEFINDER.out.legacy_json
        ch_virulencefinder_tsv = WRITE_SKIPPED_VIRULENCEFINDER.out.tsv
        ch_virulencefinder_txt = WRITE_SKIPPED_VIRULENCEFINDER.out.txt
        ch_virulencefinder_genome_seq = WRITE_SKIPPED_VIRULENCEFINDER.out.genome_seq
        ch_virulencefinder_reference_seq = WRITE_SKIPPED_VIRULENCEFINDER.out.virulence_seq
        ch_versions = ch_versions.mix(WRITE_SKIPPED_VIRULENCEFINDER.out.versions)
    }
    PLASMIDFINDER_TYPING(ch_typing_genomes)
    MLST_TYPING(ch_typing_genomes)

    if (params.bacterial_blast_reference) {
        ch_blast_reference = VALIDATE_BACTERIAL_TYPING_DATABASES.out.ready
            .map { _ready ->
                tuple(
                    [id: params.bacterial_blast_db_label],
                    file(params.bacterial_blast_reference, checkIfExists: true)
                )
            }
        BLAST_MAKEBLASTDB_TYPING(ch_blast_reference, [])
        ch_blast_db = BLAST_MAKEBLASTDB_TYPING.out.db
        blast_reference_type = 'user_supplied_reference_fasta'
        ch_versions = ch_versions.mix(BLAST_MAKEBLASTDB_TYPING.out.versions_makeblastdb)
        BLAST_BLASTN_TYPING(
            ch_typing_genomes,
            ch_blast_db,
            [],
            '',
            false
        )
        ch_blastn_output = BLAST_BLASTN_TYPING.out.txt
        ch_versions = ch_versions.mix(BLAST_BLASTN_TYPING.out.versions_blastn)
    } else if (params.bacterial_blast_db) {
        ch_blast_db = VALIDATE_BACTERIAL_TYPING_DATABASES.out.ready
            .map { _ready ->
                tuple(
                    [id: params.bacterial_blast_db_label],
                    file(params.bacterial_blast_db, checkIfExists: true)
                )
            }
        blast_reference_type = 'user_supplied_prebuilt_database'
        BLAST_BLASTN_TYPING(
            ch_typing_genomes,
            ch_blast_db,
            [],
            '',
            false
        )
        ch_blastn_output = BLAST_BLASTN_TYPING.out.txt
        ch_versions = ch_versions.mix(BLAST_BLASTN_TYPING.out.versions_blastn)
    } else {
        blast_reference_type = 'not_configured'
        WRITE_SKIPPED_BLASTN(ch_typing_genomes)
        ch_blastn_output = WRITE_SKIPPED_BLASTN.out.txt
        ch_blastn_manifest_output = ch_blastn_output.map { meta, _output ->
            tuple(meta, [])
        }
        ch_versions = ch_versions.mix(WRITE_SKIPPED_BLASTN.out.versions)
    }

    if (params.bacterial_blast_reference || params.bacterial_blast_db) {
        FILTER_BACTERIAL_BLAST_HITS(
            ch_blastn_output,
            params.bacterial_blast_min_identity,
            params.bacterial_blast_min_coverage,
            params.bacterial_blast_overlap_bp ?: 1,
            blast_filter_script
        )
        ch_blastn_filtered = FILTER_BACTERIAL_BLAST_HITS.out.filtered_hits
        ch_blastn_summary = FILTER_BACTERIAL_BLAST_HITS.out.summary
        ch_blastn_manifest_output = ch_blastn_filtered
        ch_versions = ch_versions.mix(FILTER_BACTERIAL_BLAST_HITS.out.versions)
    }

    ch_versions = ch_versions
        .mix(PLASMIDFINDER_TYPING.out.versions_plasmidfinder)
        .mix(MLST_TYPING.out.versions_mlst)

    ch_manifest_input = keyedTypingOutput(ch_virulencefinder_tsv)
        .join(keyedTypingOutput(PLASMIDFINDER_TYPING.out.tsv), by: 0)
        .map { key, meta, virulencefinder, _plasmidMeta, plasmidfinder ->
            tuple(key, meta, virulencefinder, plasmidfinder)
        }
        .join(keyedTypingOutput(MLST_TYPING.out.tsv), by: 0)
        .map { key, meta, virulencefinder, plasmidfinder, _mlstMeta, mlst ->
            tuple(key, meta, virulencefinder, plasmidfinder, mlst)
        }
        .join(keyedTypingOutput(ch_blastn_output), by: 0)
        .map { key, meta, virulencefinder, plasmidfinder, mlst, _blastMeta, blastn ->
            tuple(key, meta, virulencefinder, plasmidfinder, mlst, blastn)
        }
        .join(keyedTypingOutput(ch_blastn_manifest_output), by: 0)
        .map {
            _key,
            meta,
            virulencefinder,
            plasmidfinder,
            mlst,
            blastn,
            _filteredMeta,
            blastn_filtered ->
                tuple(
                    meta,
                    virulencefinder,
                    plasmidfinder,
                    mlst,
                    blastn,
                    blastn_filtered
                )
        }

    WRITE_GENERAL_TYPING_MANIFEST(
        ch_manifest_input,
        blast_reference_type,
        params.bacterial_blast_db_label?.toString()?.trim() ?: '',
        manifest_script
    )
    ch_versions = ch_versions.mix(WRITE_GENERAL_TYPING_MANIFEST.out.versions)

    ch_summary_input = ch_manifest_input.map {
        meta,
        virulencefinder,
        plasmidfinder,
        mlst,
        blastn,
        blastn_filtered ->
            tuple(
                meta,
                virulencefinder,
                plasmidfinder,
                mlst,
                blastn_filtered ? blastn_filtered : blastn
            )
    }
    SUMMARIZE_GENERAL_TYPING(ch_summary_input, summary_script)
    ch_versions = ch_versions.mix(SUMMARIZE_GENERAL_TYPING.out.versions)

    emit:
    database_validation = VALIDATE_BACTERIAL_TYPING_DATABASES.out.validation
    virulencefinder_json = ch_virulencefinder_json
    virulencefinder_legacy_json = ch_virulencefinder_legacy_json
    virulencefinder_tsv = ch_virulencefinder_tsv
    virulencefinder_txt = ch_virulencefinder_txt
    virulencefinder_genome_seq = ch_virulencefinder_genome_seq
    virulencefinder_reference_seq = ch_virulencefinder_reference_seq
    plasmidfinder_json = PLASMIDFINDER_TYPING.out.json
    plasmidfinder_txt = PLASMIDFINDER_TYPING.out.txt
    plasmidfinder_tsv = PLASMIDFINDER_TYPING.out.tsv
    plasmidfinder_genome_seq = PLASMIDFINDER_TYPING.out.genome_seq
    plasmidfinder_plasmid_seq = PLASMIDFINDER_TYPING.out.plasmid_seq
    mlst_tsv = MLST_TYPING.out.tsv
    blastn = ch_blastn_output
    blastn_filtered = ch_blastn_filtered
    blastn_summary = ch_blastn_summary
    evidence = SUMMARIZE_GENERAL_TYPING.out.evidence
    mlst_summary = SUMMARIZE_GENERAL_TYPING.out.mlst_summary
    manifests = WRITE_GENERAL_TYPING_MANIFEST.out.manifest
    versions = ch_versions.unique()
}
