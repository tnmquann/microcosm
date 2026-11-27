#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// Authors: @tnmquann
// Maintainers: @tnmquann

/*
 * Read-based profiling with sylph (--readbased_taxotools sylph).
 * Paired reads go straight to sylph profile (no MERGE_PAIREDENDSEQS).
 * --enable_sylphtax adds sylph-tax taxprof and a merged table in the
 * --sylph_datatype column; --readbased_postprocess phyloseq|all adds phyloseq tables.
 * run_rgi: RGI-BWT runs here when the microcosm (sourmash) branch is not selected.
 */

include { SYLPH_PROFILE } from '../../../modules/nf-core/sylph/profile/main'
include { SYLPHTAX_TAXPROF } from '../../../modules/nf-core/sylphtax/taxprof/main'
include { SYLPH_POSTPROCESS } from '../../../modules/local/sylph/postprocess/main'
include { RGI_PREPARECARDDB } from '../../../modules/local/rgi/preparecarddb/main'
include { RGI_BWT } from '../../../modules/local/rgi/bwt/main'
include { POSTPROCESS_READBASED } from '../../../modules/local/finalize/postprocess_readbased/main'
include { pathList } from '../utils_microcosm/main'

workflow READ_BASED_SYLPH {
    take:
    cleaned_reads_ch // channel: [ val(meta), [ reads ] ]
    output_prefix    // val: run prefix (outdir basename)
    run_rgi          // val: true when RGI-BWT is not run by the microcosm branch

    main:
    def versions_ch = channel.empty()
    def ch_rgi_results = channel.empty()
    def ch_sylphtax = channel.empty()
    def ch_merged = channel.empty()
    def ch_phyloseq = channel.empty()

    if (!params.sylph_database) {
        error "Missing required parameter for sylph profiling: --sylph_database"
    }

    SYLPH_PROFILE(cleaned_reads_ch, pathList(params.sylph_database, 'sylph_database'))

    if (params.enable_sylphtax) {
        if (!params.sylph_taxonomy) {
            error "--enable_sylphtax needs --sylph_taxonomy (sylph-tax metadata file(s), e.g. gtdb_r226_metadata.tsv.gz)"
        }
        def sylph_databases = pathList(params.sylph_database, 'sylph_database')
        def sylph_taxonomies = pathList(params.sylph_taxonomy, 'sylph_taxonomy')
        // sylph-tax labels a genome only from the metadata of its database; a missing file leaves that
        // database's genomes without lineage and sylph-tax does not fail, so check the pairing up front.
        if (sylph_taxonomies.size() != sylph_databases.size()) {
            error("--sylph_taxonomy has ${sylph_taxonomies.size()} file(s) but --sylph_database has ${sylph_databases.size()}: give one sylph-tax metadata file per database, in the same order (e.g. --sylph_database gtdb-r226-c200-dbv2.syl2db,fungi-refseq-2025-10-11-c200-dbv1.syldb --sylph_taxonomy gtdb_r226_metadata.tsv.gz,fungi_refseq_2025-10-11_metadata.tsv.gz)")
        }
        SYLPHTAX_TAXPROF(SYLPH_PROFILE.out.profile_out, sylph_taxonomies)
        ch_sylphtax = SYLPHTAX_TAXPROF.out.taxprof_output

        def make_phyloseq = params.readbased_postprocess in ['all', 'phyloseq']
        SYLPH_POSTPROCESS(
            ch_sylphtax
                .toSortedList { a, b -> a[0].id <=> b[0].id }
                .filter { records -> records }
                .map { records -> tuple(records.collect { record -> record[0].id }, records.collect { record -> record[1] }) },
            output_prefix,
            params.sylph_datatype,
            make_phyloseq,
            file("${projectDir}/bin/py_scripts/read_based/sylph_postprocess.py", checkIfExists: true)
        )
        ch_merged = SYLPH_POSTPROCESS.out.merged
        ch_phyloseq = SYLPH_POSTPROCESS.out.phyloseq
    }

    if (run_rgi && params.enable_rgi_bwt) {
        def ch_rgi_db
        if (params.rgi_preparecarddb_dir) {
            ch_rgi_db = channel.value(file(params.rgi_preparecarddb_dir))
        } else {
            RGI_PREPARECARDDB([])
            ch_rgi_db = RGI_PREPARECARDDB.out.db
            versions_ch = versions_ch.mix(RGI_PREPARECARDDB.out.versions)
        }
        RGI_BWT(cleaned_reads_ch.map { meta, reads -> [meta, reads[0], reads[1]] }, ch_rgi_db)
        versions_ch = versions_ch.mix(RGI_BWT.out.versions)
        ch_rgi_results = RGI_BWT.out.outdir

        // Collecting the RGI-BWT folders of all samples needs no taxonomy table, so it also runs without sourmash + YACHT
        if (params.readbased_postprocess in ['rgi_bwt', 'all']) {
            POSTPROCESS_READBASED(
                RGI_BWT.out.outdir
                    .map { _meta, rgi_dir -> rgi_dir }
                    .collect()
                    .map { rgi_dirs -> [[id: output_prefix], output_prefix, 'rgi_bwt', false, false, ' ', [], rgi_dirs] },
                params.microcosm_datatype,
                file("${projectDir}/bin/py_scripts/read_based/post_processing.py", checkIfExists: true)
            )
            versions_ch = versions_ch.mix(POSTPROCESS_READBASED.out.versions)
        }
    }

    emit:
    profile = SYLPH_PROFILE.out.profile_out
    sylphtax = ch_sylphtax
    merged = ch_merged
    phyloseq = ch_phyloseq
    rgi_results = ch_rgi_results
    versions = versions_ch
}
