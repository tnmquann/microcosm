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

/*
 * --bacterial_decision_table (optional): lineage_pattern at any rank (p__, g__, s__ ...), optional
 * label, then one y/n column per step. The most specific matching row decides its y/n steps for
 * that genome; a blank cell, an absent
 * column or a lineage with no matching row follows the default microcosm route.
 * --skip_<tool> still turns a tool off everywhere.
 */
def bacterialDecisionSteps() {
    return ['amrfinderplus', 'abricate', 'rgi', 'resfinder', 'deeparg', 'mlst', 'virulencefinder', 'plasmidfinder',
            'kleborate', 'sistr', 'seqsero2', 'genotyphi', 'paratype']
}

def readDecisionTable(path) {
    def steps = bacterialDecisionSteps()
    def rows = readLineageMap(path, [], [])
    def unknown = rows ? rows[0].keySet() - (['lineage_pattern', 'label'] + steps) : []
    if (unknown) {
        error("${path}: unknown column(s) ${unknown.join(', ')}; step columns are ${steps.join(', ')}")
    }
    rows.eachWithIndex { row, index ->
        if (patternRankDepth(row.lineage_pattern.trim()) < 0) {
            error("${path}: lineage_pattern '${row.lineage_pattern}' needs a GTDB rank prefix, e.g. g__Escherichia (row ${index + 2})")
        }
        steps.each { step ->
            if (!((row[step] ?: '').trim().toLowerCase() in ['', 'y', 'n', 'yes', 'no'])) {
                error("${path}: ${step} must be y or n, got '${row[step]}' (row ${index + 2})")
            }
        }
    }
    return rows
}

// Rank depth of a pattern (d__ = 0 ... s__ = 6) from its last rank prefix; -1 without one.
def patternRankDepth(pattern) {
    def prefixes = (pattern =~ /(?:^|;)\s*([dpcofgs])__/).collect { match -> match[1] }
    return prefixes ? 'dpcofgs'.indexOf(prefixes[-1]) : -1
}

// The most specific matching row wins (s__ over g__ over p__ ...; same rank: first row), so rows can
// mix ranks in any order. A pattern matches whole taxa (one or more, e.g. f__X;g__Y): g__Bacillus
// matches g__Bacillus_A (GTDB suffixes ignored as in lookupLineage) but not g__Bacillusx.
def lookupDecisionRow(rows, lineage) {
    def taxa = ';' + stripGtdbSuffixes(lineage).split(';').collect { taxon -> taxon.trim() }.join(';') + ';'
    def matches = rows.findAll { row ->
        def pattern = stripGtdbSuffixes(row.lineage_pattern.trim()).split(';').collect { taxon -> taxon.trim() }.join(';')
        taxa.contains(';' + pattern + ';')
    }
    return matches ? matches.max { row -> patternRankDepth(row.lineage_pattern.trim()) } : null
}

// [pattern: matched lineage_pattern or '', steps: [step: true|false] for the y/n cells of that row]
def decisionFor(rows, lineage) {
    def row = lookupDecisionRow(rows, lineage)
    def steps = row ? bacterialDecisionSteps()
        .findAll { step -> row[step]?.trim() }
        .collectEntries { step -> [(step): row[step].trim().toLowerCase().startsWith('y')] } : [:]
    return [pattern: row ? row.lineage_pattern.trim() : '', steps: steps]
}

// Decision for one genome and step; defaultRoute is what microcosm does without a table.
def stepWanted(meta, step, defaultRoute = true) {
    def decided = meta.decision_steps?.get(step)
    return decided == null ? defaultRoute : decided
}
