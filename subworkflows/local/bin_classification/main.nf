#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// Import nf-core modules
include { SOURMASH_TAXANNOTATE }            from '../../../modules/nf-core/sourmash/taxannotate/main'

// Import local modules
include { CREATE_BINS_CSV }                         from '../../../modules/local/sourmash/createbinscsv/main'
include { SOURMASH_MANYSKETCH as SOURMASH_MANYSKETCH_BINS }     from '../../../modules/local/sourmash/manysketch/main'
include { EXTRACT_SOURMASH_SINGLESKETCHES }         from '../../../modules/local/finalize/extract_sourmash_singlesketches/main'
include { SOURMASH_FASTMULTIGATHER as SOURMASH_FASTMULTIGATHER_BINS } from '../../../modules/local/sourmash/fastmultigather/main'
include { COMBINE_FASTMULTIGATHER_RESULTS } from '../../../modules/local/finalize/combine_fastmultigather_results/main'
include { COMBINE_FASTMULTIGATHER_RESULTS as COMBINE_TAXANNOTATE_RESULTS} from '../../../modules/local/finalize/combine_fastmultigather_results/main'
include { COMBINE_FASTMULTIGATHER_RESULTS as COMBINE_TAXGENOME_RESULTS} from '../../../modules/local/finalize/combine_fastmultigather_results/main'
include { SOURMASH_TAXGENOME }              from '../../../modules/local/sourmash/taxgenome/main'
include { GTDBTK_CLASSIFYWF }               from '../../../modules/nf-core/gtdbtk/classifywf/main'
include { GTDBTK_TAXGENOME }                from '../../../modules/local/gtdbtk/taxgenome/main'
include { COMBINE_FASTMULTIGATHER_RESULTS as COMBINE_GTDBTK_RESULTS } from '../../../modules/local/finalize/combine_fastmultigather_results/main'
include { pathList; binClassifier } from '../utils_microcosm/main'

// GTDB-Tk --extension of a bin file (all genomes of one GTDB-Tk call must share it)
def gtdbtkExtension(genome) {
    def match = genome.name =~ /(?i)\.((?:fa|fasta|fna|fas)(?:\.gz)?)$/
    if (!match.find()) {
        error("GTDB-Tk: cannot tell the FASTA extension of ${genome.name} (expected .fa/.fasta/.fna/.fas, optionally .gz)")
    }
    return match.group(1)
}

workflow BIN_CLASSIFICATION {
    take:
    bins_ch  // channel: [ val(meta), path(bin) ] - individual bins from BINNING/BINNING_REFINEMENT

    main:
    ch_versions = channel.empty()
    ch_bins_csv = channel.empty()
    ch_sketches = channel.empty()
    ch_single_sketches = channel.empty()
    ch_gather_csv = channel.empty()
    ch_taxannotate = channel.empty()
    ch_combined_gather = channel.empty()
    ch_combined_taxannotate = channel.empty()
    ch_genome_classification = channel.empty()
    ch_combined_taxgenome = channel.empty()

    if (binClassifier() == 'gtdbtk') {
        // --bin_classifier gtdbtk (idea from nf-core/mag): GTDB-Tk classify_wf instead of sourmash.
        // One call per sample (all binners/refiners together; --genome_input rows in one call), split by
        // FASTA extension. Summaries are converted to the sourmash tax genome CSV layout, so the routing,
        // decision table and MultiQC read the same columns (query_name = bin name without extension).
        if (!params.gtdbtk_db) {
            error "Missing required parameter for --bin_classifier gtdbtk: --gtdbtk_db (GTDB-Tk reference data directory, e.g. release226)"
        }
        def ch_gtdbtk_input = bins_ch
            .map { meta, bin ->
                def batch = meta.source_route == 'standalone' ? 'genome_input' : meta.id
                def extension = gtdbtkExtension(bin)
                [[batch, extension], bin]
            }
            .groupTuple()
            .map { key, bins -> [[id: key[0], gtdbtk_extension: key[1]], bins] }

        GTDBTK_CLASSIFYWF(
            ch_gtdbtk_input,
            channel.value(['gtdb', file(params.gtdbtk_db, checkIfExists: true)]),
            params.gtdbtk_pplacer_scratch
        )

        GTDBTK_TAXGENOME(
            GTDBTK_CLASSIFYWF.out.summary,
            file("${projectDir}/bin/py_scripts/assembly_based/gtdbtk_to_taxgenome.py", checkIfExists: true)
        )
        ch_versions = ch_versions.mix(GTDBTK_TAXGENOME.out.versions)
        ch_genome_classification = GTDBTK_TAXGENOME.out.genome_classification

        COMBINE_GTDBTK_RESULTS(
            ch_genome_classification
                .map { meta, csv -> [meta.id, csv] }
                .groupTuple()
                .map { id, csvs -> [[id: id, file_type: 'taxgenome'], csvs.flatten()] }
        )
        ch_versions = ch_versions.mix(COMBINE_GTDBTK_RESULTS.out.versions)
        ch_combined_taxgenome = COMBINE_GTDBTK_RESULTS.out.combined_csv
    } else {
        // Validate required parameters
        if (!params.sourmash_database) {
            error "Missing required parameter for bin classification: --sourmash_database"
        }
        if (!params.sourmash_taxonomy_csv) {
            error "Missing required parameter for bin classification: --sourmash_taxonomy_csv (path to Sourmash taxonomy lineages CSV)"
        }

        // Step 1: Group bins by assembly/binner/binrefine/sample to create CSV manifests
        ch_bins_grouped = bins_ch
            .map { meta, bin ->
                def meta_clean = meta.clone()
                meta_clean.remove('bin_id')
                meta_clean.remove('domain')
                [meta_clean, bin]
            }
            .groupTuple()

        // Create CSV files for each group
        CREATE_BINS_CSV(ch_bins_grouped)
        ch_versions = ch_versions.mix(CREATE_BINS_CSV.out.versions)

        // Step 2: Run sourmash manysketch on CSV files
        SOURMASH_MANYSKETCH_BINS(CREATE_BINS_CSV.out.bins_csv)
        ch_versions = ch_versions.mix(SOURMASH_MANYSKETCH_BINS.out.versions)

        // Use multiMap to split the sketch outputs for parallel processing
        ch_sketch_split = SOURMASH_MANYSKETCH_BINS.out.sketch_zip_file
            .multiMap { meta, sketch_zip ->
                for_fastmultigather: [meta, sketch_zip]
                for_extract: [meta, sketch_zip]
            }

        // Step 3: Run fastmultigather
        SOURMASH_FASTMULTIGATHER_BINS(
            ch_sketch_split.for_fastmultigather,
            pathList(params.sourmash_database, 'sourmash_database')
        )
        ch_versions = ch_versions.mix(SOURMASH_FASTMULTIGATHER_BINS.out.versions)

        // Step 9: Extract single sketches
        sourmash_extract_zip_script = file("${projectDir}/bin/py_scripts/read_based/sourmash_extract_zip_with_ksize.py", checkIfExists: true)
        EXTRACT_SOURMASH_SINGLESKETCHES(
            ch_sketch_split.for_extract,
            params.sourmash_ksize,
            sourmash_extract_zip_script
        )
        ch_versions = ch_versions.mix(EXTRACT_SOURMASH_SINGLESKETCHES.out.versions)

        // Skip sourmash tax steps for groups with no gather hits.
        ch_gather_split = SOURMASH_FASTMULTIGATHER_BINS.out.gather_csv.branch { _meta, gather_csv ->
            nonempty: gather_csv.toFile().length() > 0
            empty: true
        }

        // Step 4: Run sourmash tax annotate
        // Input: gather_csv from step 3
        SOURMASH_TAXANNOTATE(
            ch_gather_split.nonempty,
            pathList(params.sourmash_taxonomy_csv, 'sourmash_taxonomy_csv')
        )
        ch_versions = ch_versions.mix(SOURMASH_TAXANNOTATE.out.versions_sourmash)

        // Step 5: Run sourmash tax genome
        // Input: TAXANNOTATE output (with-lineages.csv.gz)
        SOURMASH_TAXGENOME(
            SOURMASH_TAXANNOTATE.out.result,
            pathList(params.sourmash_taxonomy_csv, 'sourmash_taxonomy_csv')
        )
        ch_versions = ch_versions.mix(SOURMASH_TAXGENOME.out.versions)

        // Step 6: Combine fastmultigather results by sample ID
        ch_gather_for_combine = SOURMASH_FASTMULTIGATHER_BINS.out.gather_csv
            .map { meta, csv ->
                def meta_combine = [:]
                meta_combine.id = meta.id
                meta_combine.file_type = 'gather'
                [meta_combine.id, meta_combine, csv]
            }
            .groupTuple(by: 0)
            .map { _id, metas, csvs ->
                [metas[0], csvs.flatten()]
            }

        COMBINE_FASTMULTIGATHER_RESULTS(ch_gather_for_combine)
        ch_versions = ch_versions.mix(COMBINE_FASTMULTIGATHER_RESULTS.out.versions)

        // Step 7: Combine taxannotate results by sample ID
        // Note: Decompress .gz files before combining
        ch_taxannotate_for_combine = SOURMASH_TAXANNOTATE.out.result
            .map { meta, csv_gz ->
                def meta_combine = [:]
                meta_combine.id = meta.id
                meta_combine.file_type = 'taxannotate'
                [meta_combine.id, meta_combine, csv_gz]
            }
            .groupTuple(by: 0)
            .map { _id, metas, csvs ->
                [metas[0], csvs.flatten()]
            }

        COMBINE_TAXANNOTATE_RESULTS(ch_taxannotate_for_combine)
        ch_versions = ch_versions.mix(COMBINE_TAXANNOTATE_RESULTS.out.versions)

        // Step 8: Combine taxgenome results by sample ID
        ch_taxgenome_for_combine = SOURMASH_TAXGENOME.out.genome_classification
            .map { meta, csv ->
                def meta_combine = [:]
                meta_combine.id = meta.id
                meta_combine.file_type = 'taxgenome'
                [meta_combine.id, meta_combine, csv]
            }
            .groupTuple(by: 0)
            .map { _id, metas, csvs ->
                [metas[0], csvs.flatten()]
            }

        COMBINE_TAXGENOME_RESULTS(ch_taxgenome_for_combine)
        ch_versions = ch_versions.mix(COMBINE_TAXGENOME_RESULTS.out.versions)

        ch_bins_csv = CREATE_BINS_CSV.out.bins_csv
        ch_sketches = SOURMASH_MANYSKETCH_BINS.out.sketch_zip_file
        ch_single_sketches = EXTRACT_SOURMASH_SINGLESKETCHES.out.sig_zip_files
        ch_gather_csv = SOURMASH_FASTMULTIGATHER_BINS.out.gather_csv
        ch_taxannotate = SOURMASH_TAXANNOTATE.out.result
        ch_genome_classification = SOURMASH_TAXGENOME.out.genome_classification
        ch_combined_gather = COMBINE_FASTMULTIGATHER_RESULTS.out.combined_csv
        ch_combined_taxannotate = COMBINE_TAXANNOTATE_RESULTS.out.combined_csv
        ch_combined_taxgenome = COMBINE_TAXGENOME_RESULTS.out.combined_csv
    }

    emit:
    bins_csv             = ch_bins_csv
    sketches             = ch_sketches
    single_sketches      = ch_single_sketches
    gather_csv           = ch_gather_csv
    taxannotate          = ch_taxannotate
    genome_classification = ch_genome_classification  // [meta, CSV in the sourmash tax genome layout]
    combined_gather      = ch_combined_gather
    combined_taxannotate = ch_combined_taxannotate
    combined_taxgenome   = ch_combined_taxgenome
    versions             = ch_versions
}
