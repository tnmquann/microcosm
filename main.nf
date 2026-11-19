#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { METAFLOW } from './workflow/metaflow'

// Typed declarations for the non-string params (defaults live in conf/*.config, which win).
// With the strict syntax (NXF_SYNTAX_PARSER=v2) undeclared CLI values reach scripts as
// strings, so `--skip_binqc false` was the truthy string "false"; the types convert them.
params {
    annotation_bakta_activate_plot: Boolean = false
    annotation_bakta_complete: Boolean = false
    annotation_bakta_compliant: Boolean = false
    annotation_bakta_crispr: Boolean = false
    annotation_bakta_gap: Boolean = false
    annotation_bakta_ncrna: Boolean = false
    annotation_bakta_ncrnaregion: Boolean = false
    annotation_bakta_ori: Boolean = false
    annotation_bakta_pseudo: Boolean = false
    annotation_bakta_renamecontigheaders: Boolean = false
    annotation_bakta_rrna: Boolean = false
    annotation_bakta_singlemode: Boolean = false
    annotation_bakta_skipcds: Boolean = false
    annotation_bakta_skipsorf: Boolean = false
    annotation_bakta_tmrna: Boolean = false
    annotation_bakta_trna: Boolean = false
    annotation_prokka_addgenes: Boolean = false
    annotation_prokka_cdsrnaolap: Boolean = false
    annotation_prokka_compliant: Boolean = true
    annotation_prokka_rawproduct: Boolean = false
    annotation_prokka_retaincontigheaders: Boolean = false
    annotation_prokka_rnammer: Boolean = false
    annotation_prokka_singlemode: Boolean = false
    annotation_pyrodigal_closed: Boolean = false
    annotation_pyrodigal_forcenonsd: Boolean = false
    annotation_pyrodigal_singlemode: Boolean = false
    annotation_pyrodigal_usespecialstopcharacter: Boolean = false
    bin_domain_classification: Boolean = false
    busco_clean: Boolean = true
    cleanup: Boolean = false
    enable_bacterial_downstream: Boolean = false
    enable_copyintermediate: Boolean = true
    enable_copymergedseqs: Boolean = true
    enable_copysketch: Boolean = true
    enable_copysketch_singles: Boolean = true
    enable_readbase: Boolean = false
    enable_rgi_bwt: Boolean = false
    enable_singlesketch: Boolean = false
    enable_sylphtax: Boolean = false
    exclude_unbins_from_postbinning: Boolean = false
    fastp_save_merged: Boolean = false
    fastp_save_trimmed_fail: Boolean = false
    gtdbtk_pplacer_scratch: Boolean = true
    gunc_save_db: Boolean = false
    help: Boolean = false
    help_full: Boolean = false
    megahit_fix_cpu_1: Boolean = false
    prokka_bins_fast_mode: Boolean = false
    prokka_bins_with_compliance: Boolean = false
    rgi_include_other_models: Boolean = false
    rgi_include_wildcard: Boolean = true
    run_gunc: Boolean = false
    save_amrfinderplus_db: Boolean = true
    save_assembly_mapped_reads: Boolean = true
    save_bakta_db: Boolean = true
    save_bin_annotations: Boolean = true
    save_busco_db: Boolean = false
    save_checkm2_data: Boolean = false
    save_checkm_data: Boolean = false
    save_contig_annotations: Boolean = true
    save_deeparg_db: Boolean = true
    save_genomad_db: Boolean = true
    save_paratype_db: Boolean = true
    save_resfinder_db: Boolean = true
    save_rgi_db: Boolean = true
    save_virulencefinder_db: Boolean = true
    show_hidden: Boolean = false
    skip_abricate: Boolean = false
    skip_amrfinderplus: Boolean = false
    skip_argnorm: Boolean = false
    skip_bin_annotation: Boolean = false
    skip_bin_classification: Boolean = false
    skip_binning: Boolean = false
    skip_binning_bamabund: Boolean = false
    skip_binning_refinement: Boolean = false
    skip_binqc: Boolean = false
    skip_comebin: Boolean = false
    skip_concoct: Boolean = true
    skip_contig_analysis: Boolean = false
    skip_contig_annotation: Boolean = false
    skip_deeparg: Boolean = false
    skip_genomad_e2e: Boolean = false
    skip_genotyphi: Boolean = false
    skip_kleborate: Boolean = false
    skip_maxbin2: Boolean = false
    skip_megahit: Boolean = false
    skip_metabat2: Boolean = false
    skip_mlst: Boolean = false
    skip_multiqc: Boolean = false
    skip_paratype: Boolean = false
    skip_plasmidfinder: Boolean = false
    skip_preprocess: Boolean = false
    skip_quast: Boolean = false
    skip_resfinder: Boolean = false
    skip_rgi: Boolean = false
    skip_semibin: Boolean = false
    skip_seqsero2: Boolean = false
    skip_sistr: Boolean = false
    skip_skani_search: Boolean = true
    skip_spades: Boolean = false
    skip_tiara_contigqc: Boolean = false
    skip_vamb: Boolean = true
    skip_virulencefinder: Boolean = false
    skip_yacht: Boolean = false
    annotation_prokka_evalue: Float = 0.000001
    genomad_min_score: Float = 0.7
    refine_bins_dastool_threshold: Float = 0.5
    virulencefinder_min_coverage: Float = 0.6
    virulencefinder_min_identity: Float = 0.9
    annotation_bakta_mincontiglen: Integer = 1
    annotation_bakta_translationtable: Integer = 11
    annotation_prokka_coverage: Integer = 80
    annotation_prokka_gcode: Integer = 11
    annotation_prokka_mincontiglen: Integer = 1
    annotation_pyrodigal_transtable: Integer = 11
    checkm2_db_version: Integer = 14897628
    genomad_splits: Integer = 1
    gtdbtk_pplacer_cpus: Integer = 1
    max_cpus: Integer = 25
    max_unbinned_contigs: Integer = 100
    metabat_rng_seed: Integer = 1
    min_contig_size: Integer = 1500
    min_length_unbinned_contigs: Integer = 1000000
    refine_bins_binette_contamination_weight: Integer = 5
    refine_bins_binette_min_completeness: Integer = 40
    sourmash_ksize: Integer = 31
    spades_fix_cpus: Integer = -1
    tiara_min_length: Integer = 3000
    yacht_ksize: Integer = 31
}

workflow {
    if (!(workflow.preview && !params.input && !params.assembly_input && !params.genome_input)) {
        METAFLOW()
    }
}
