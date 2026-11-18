#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// Authors: @tnmquann
// Maintainers: @tnmquann

/*
 * Lineage lookup tables for the bacterial downstream layer
 * (--kleborate_preset_map, --amrfinderplus_organism_map).
 * A table is a TSV with a lineage_pattern column plus value columns; rows are
 * tried in order and the first pattern found inside the GTDB lineage wins.
 * Both sides are matched with GTDB placeholder suffixes removed
 * (s__Enterococcus_B faecium -> s__Enterococcus faecium); meta.lineage itself is left unchanged.
 */

// GTDB rs226: every '_' inside a taxon name is a placeholder suffix of upper-case letters.
def stripGtdbSuffixes(lineage) {
    return (lineage ?: '').replaceAll('(?<=[A-Za-z0-9])_[A-Z]+(?=[\\s;]|$)', '')
}

def readLineageMap(path, required_columns, non_empty_columns) {
    def rows = file(path, checkIfExists: true).splitCsv(header: true, sep: '\t')
    rows.eachWithIndex { row, index ->
        def missing = (['lineage_pattern'] + required_columns).findAll { column -> !row.containsKey(column) }
        if (missing) {
            error("${path}: missing column(s) ${missing.join(', ')} (row ${index + 2})")
        }
        def empty = (['lineage_pattern'] + non_empty_columns).findAll { column -> !row[column]?.trim() }
        if (empty) {
            error("${path}: empty ${empty.join(', ')} (row ${index + 2})")
        }
    }
    return rows
}

def lookupLineage(rows, lineage) {
    def text = stripGtdbSuffixes(lineage)
    return rows.find { row -> text.contains(stripGtdbSuffixes(row.lineage_pattern.trim())) }
}
