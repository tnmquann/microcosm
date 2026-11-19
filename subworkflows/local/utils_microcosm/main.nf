#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// Authors: @tnmquann
// Maintainers: @tnmquann

/*
 * Shared helpers for microcosm subworkflows.
 */

// Comma-separated path parameter -> list of checked files
// (e.g. --sourmash_database "db1.zip,db2.zip", --sylph_database "a.syldb,b.syldb").
def pathList(value, param_name) {
    def paths = (value ?: '').toString().tokenize(',').collect { path_text -> path_text.trim() }.findAll { path_text -> path_text }
    if (!paths) {
        error("--${param_name} is empty")
    }
    return paths.collect { path_text -> file(path_text, checkIfExists: true) }
}

// --bin_classifier: sourmash (default) or gtdbtk; both give GTDB lineages to the bacterial routing.
def binClassifier() {
    def classifier = (params.bin_classifier ?: 'sourmash').toString().trim().toLowerCase()
    if (!(classifier in ['sourmash', 'gtdbtk'])) {
        error("Unsupported --bin_classifier '${params.bin_classifier}'. Use sourmash or gtdbtk.")
    }
    return classifier
}

// Database inputs of the chosen bin classifier are set (--genome_input rows without lineage need them).
def binClassificationAvailable() {
    return binClassifier() == 'gtdbtk'
        ? params.gtdbtk_db as boolean
        : (params.sourmash_database && params.sourmash_taxonomy_csv) as boolean
}

// YACHT takes one database, so several --sourmash_database entries imply --skip_yacht.
def yachtSkipped() {
    def database_count = (params.sourmash_database ?: '').toString().tokenize(',').findAll { path_text -> path_text.trim() }.size()
    return params.skip_yacht || database_count > 1
}

// --readbased_taxotools: 'microcosm' (sourmash + YACHT), 'sylph', or both comma-separated.
def readbasedTaxotools() {
    def tools = (params.readbased_taxotools ?: 'microcosm').toString().tokenize(',').collect { tool -> tool.trim().toLowerCase() }.findAll { tool -> tool }.unique()
    def unknown = tools.findAll { tool -> !(tool in ['microcosm', 'sylph']) }
    if (!tools || unknown) {
        error("Unsupported --readbased_taxotools '${params.readbased_taxotools}'. Use microcosm, sylph, or microcosm,sylph.")
    }
    return tools
}

// --binqc_tool: one or more of busco, checkm, checkm2, comma-separated (e.g. "checkm2,checkm").
def binqcTools() {
    def tools = (params.binqc_tool ?: '').toString().tokenize(',').collect { tool -> tool.trim().toLowerCase() }.findAll { tool -> tool }.unique()
    def unknown = tools.findAll { tool -> !(tool in ['busco', 'checkm', 'checkm2']) }
    if (!tools || unknown) {
        error("Unsupported --binqc_tool '${params.binqc_tool}'. Use busco, checkm, checkm2, or several of them comma-separated.")
    }
    return tools
}
