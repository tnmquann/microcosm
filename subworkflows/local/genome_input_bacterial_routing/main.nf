#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// Authors: @tnmquann
// Maintainers: @tnmquann

include { VALIDATE_GENOME_INPUT } from '../../../modules/local/input/validate_genome_input/main'
include { NORMALIZE_GENOME_CLASSIFICATIONS } from '../../../modules/local/finalize/normalize_genome_classifications/main'
include { PUBLISH_GENOME_ROUTING } from '../../../modules/local/finalize/publish_genome_routing/main'
// Intentional nested use: standalone --genome_input rows with blank lineage
// are classified here; assembly-produced bins are classified by workflow/microcosm.nf.
include { BIN_CLASSIFICATION } from '../bin_classification/main'
include { readDecisionTable; decisionFor } from '../utils_bacterial_lineage/main'
include { binClassificationAvailable; binClassifier } from '../utils_microcosm/main'

def routingTsvValue(value) {
    return value == null ? '' : value.toString().replaceAll(/[\t\r\n]/, ' ')
}

def routingTsvLine(meta, genome, status, reason) {
    return [
        meta.sample_id,
        meta.genome_id,
        meta.run_id,
        meta.group,
        meta.assembler,
        genome,
        meta.lineage,
        meta.domain,
        meta.phylum,
        meta['class'],
        meta.order,
        meta.family,
        meta.genus,
        meta.species,
        meta.lineage_source,
        meta.source_route,
        status,
        reason,
        meta.decision_pattern,
        meta.decision_steps?.collect { step, run -> "${step}=${run ? 'y' : 'n'}" }?.join(';'),
        meta.decision_params?.collect { column, value -> "${column}=${value}" }?.join('; ')
    ].collect { value -> routingTsvValue(value) }.join('\t') + '\n'
}

workflow GENOME_INPUT_AND_BACTERIAL_ROUTING {
    take:
    assembly_genomes_ch
    assembly_classifications_ch

    main:
    ch_versions = channel.empty()
    ch_provided_genomes = channel.empty()
    ch_classification_candidates = channel.empty()
    ch_classification_csvs = channel.empty()
    pipeline_root = params.pipeline_root ?: projectDir

    if (params.genome_input) {
        validator_script = file(
            "${pipeline_root}/bin/py_scripts/bacterial_analysis/validate_genome_input.py",
            checkIfExists: true
        )
        classification_available = binClassificationAvailable()

        VALIDATE_GENOME_INPUT(
            file(params.genome_input, checkIfExists: true),
            classification_available,
            validator_script
        )
        ch_versions = ch_versions.mix(VALIDATE_GENOME_INPUT.out.versions)

        ch_standalone_genomes = VALIDATE_GENOME_INPUT.out.validated_csv
            .splitCsv(header: true, sep: ',', strip: false)
            .map { row ->
                def meta = [
                    id: row.sample_id,
                    sample_id: row.sample_id,
                    genome_id: row.sample_id,
                    run_id: row.run_id ?: null,
                    group: row.group ?: null,
                    assembler: row.assembler ?: null,
                    lineage: row.lineage ?: '',
                    domain: row.domain ?: '',
                    phylum: row.phylum ?: '',
                    'class': row['class'] ?: '',
                    order: row.order ?: '',
                    family: row.family ?: '',
                    genus: row.genus ?: '',
                    species: row.species ?: '',
                    lineage_source: row.lineage_source,
                    source_route: row.source_route,
                    classification_key: row.classification_key,
                    input_row: row.row_number as Integer,
                    classification_status: row.lineage_source == 'provided' ? 'provided' : ''
                ]
                tuple(meta, file(row.fasta, checkIfExists: true))
            }

        ch_standalone_split = ch_standalone_genomes.branch { meta, _genome ->
            provided: meta.lineage_source == 'provided'
            classify: meta.lineage_source == 'classified'
        }
        ch_provided_genomes = ch_standalone_split.provided
        ch_classification_candidates = ch_standalone_split.classify

        if (classification_available) {
            BIN_CLASSIFICATION(ch_classification_candidates)
            ch_classification_csvs = BIN_CLASSIFICATION.out.genome_classification
            ch_versions = ch_versions.mix(BIN_CLASSIFICATION.out.versions)
        }
    } else {
        ch_classification_candidates = assembly_genomes_ch.map { meta, genome ->
            def meta_new = meta.clone()
            // Same name CREATE_BINS_CSV gives the sketch (x.fa.gz -> x), i.e. sourmash query_name.
            def bin_name = genome.name.replaceFirst('\\.gz$', '').replaceFirst('\\.(fa|fasta|fna|fas)$', '')
            meta_new.sample_id = meta.sample_id ?: meta.id
            meta_new.genome_id = meta.bin_id
                ? meta.bin_id.toString().replaceFirst('\\.gz$', '').replaceFirst('\\.(fa|fasta|fna|fas)$', '')
                : bin_name
            meta_new.run_id = meta.containsKey('run_id') ? meta.run_id : null
            meta_new.group = meta.containsKey('group') ? meta.group : null
            meta_new.assembler = meta.containsKey('assembler') ? meta.assembler : null
            meta_new.lineage = ''
            meta_new.domain = ''
            meta_new.phylum = ''
            meta_new['class'] = ''
            meta_new.order = ''
            meta_new.family = ''
            meta_new.genus = ''
            meta_new.species = ''
            meta_new.lineage_source = 'classified'
            meta_new.source_route = 'assembly_bins'
            meta_new.classification_key = bin_name
            meta_new.classification_status = ''
            tuple(meta_new, genome)
        }
        ch_classification_csvs = assembly_classifications_ch
    }

    empty_classification = file(
        "${pipeline_root}/assets/empty_genome_classification.csv",
        checkIfExists: true
    )
    normalizer_script = file(
        "${pipeline_root}/bin/py_scripts/bacterial_analysis/normalize_genome_classifications.py",
        checkIfExists: true
    )
    ch_classification_files = ch_classification_csvs
        .map { _meta, classification_csv -> classification_csv }
        .mix(channel.of(empty_classification))
        .collect()

    NORMALIZE_GENOME_CLASSIFICATIONS(ch_classification_files, normalizer_script)
    ch_versions = ch_versions.mix(NORMALIZE_GENOME_CLASSIFICATIONS.out.versions)

    // Unresolved rows are joined too, so their classifier detail reaches the routing report.
    ch_resolved_classification_rows = NORMALIZE_GENOME_CLASSIFICATIONS.out.resolved_tsv
        .mix(NORMALIZE_GENOME_CLASSIFICATIONS.out.unresolved_tsv)
        .splitCsv(header: true, sep: '\t', strip: false)
        .map { row -> tuple(row.query_name, row) }

    ch_candidate_records = ch_classification_candidates
        .map { meta, genome -> tuple(meta.classification_key, meta, genome) }
        .join(ch_resolved_classification_rows, by: 0, remainder: true)
        // Drop classification rows with no matching candidate genome.
        .filter { record -> record[1] != null }
        .branch { record ->
            resolved: record.size() == 4 && record[3]?.lineage
            unresolved: true
        }

    ch_classified_resolved = ch_candidate_records.resolved.map { record ->
        def meta = record[1]
        def genome = record[2]
        def lineage_row = record[3]
        def meta_new = meta.clone()
        meta_new.lineage = lineage_row.lineage
        meta_new.domain = lineage_row.domain
        meta_new.phylum = lineage_row.phylum
        meta_new['class'] = lineage_row['class']
        meta_new.order = lineage_row.order
        meta_new.family = lineage_row.family
        meta_new.genus = lineage_row.genus
        meta_new.species = lineage_row.species
        meta_new.lineage_source = 'classified'
        meta_new.classification_status = lineage_row.status
        tuple(meta_new, genome)
    }

    ch_unresolved_genomes = ch_candidate_records.unresolved.map { record ->
        def meta = record[1]
        def genome = record[2]
        def meta_new = meta.clone()
        meta_new.classification_status = 'unresolved'
        meta_new.unresolved_detail = record[3] == null ? 'no classification row' : (record[3].detail ?: record[3].status)
        tuple(meta_new, genome)
    }

    ch_resolved_genomes = ch_provided_genomes
        .mix(ch_classified_resolved)

    ch_routed_genomes = ch_resolved_genomes.branch { meta, _genome ->
        bacterial: meta.domain == 'd__Bacteria'
        excluded: true
    }

    // --bacterial_decision_table: per-lineage y/n steps (meta keys only added when a table is given,
    // so runs without one keep their -resume cache)
    decision_rows = params.bacterial_decision_table ? readDecisionTable(params.bacterial_decision_table) : null
    ch_bacterial_genomes = decision_rows == null ? ch_routed_genomes.bacterial : ch_routed_genomes.bacterial.map { meta, genome ->
        def decision = decisionFor(decision_rows, meta.lineage)
        def meta_new = meta.clone()
        meta_new.decision_pattern = decision.pattern
        meta_new.decision_steps = decision.steps
        meta_new.decision_params = decision.params
        tuple(meta_new, genome)
    }

    routing_header = [
        'sample_id',
        'genome_id',
        'run_id',
        'group',
        'assembler',
        'fasta',
        'lineage',
        'domain',
        'phylum',
        'class',
        'order',
        'family',
        'genus',
        'species',
        'lineage_source',
        'source_route',
        'status',
        'reason',
        'decision_pattern',
        'decision_steps',
        'decision_params'
    ].join('\t') + '\n'

    ch_bacterial_report = ch_bacterial_genomes
        .map { meta, genome ->
            routingTsvLine(
                meta,
                genome,
                meta.classification_status,
                'eligible_bacterial_genome'
            )
        }
        .ifEmpty('')
        .collectFile(
            name: 'bacterial_genomes.routing.tsv',
            seed: routing_header,
            sort: true
        )

    ch_excluded_report = ch_routed_genomes.excluded
        .map { meta, genome ->
            routingTsvLine(
                meta,
                genome,
                meta.classification_status,
                'resolved_domain_is_not_d__Bacteria'
            )
        }
        .ifEmpty('')
        .collectFile(
            name: 'excluded_genomes.routing.tsv',
            seed: routing_header,
            sort: true
        )

    ch_unresolved_report = ch_unresolved_genomes
        .map { meta, genome ->
            routingTsvLine(
                meta,
                genome,
                meta.classification_status,
                "no_resolved_lineage_from_BIN_CLASSIFICATION (${binClassifier()}: ${meta.unresolved_detail ?: 'no detail'})"
            )
        }
        .ifEmpty('')
        .collectFile(
            name: 'unresolved_genomes.routing.tsv',
            seed: routing_header,
            sort: true
        )

    PUBLISH_GENOME_ROUTING(
        ch_bacterial_report,
        ch_excluded_report,
        ch_unresolved_report
    )
    ch_versions = ch_versions.mix(PUBLISH_GENOME_ROUTING.out.versions)

    emit:
    bacterial_genomes = ch_bacterial_genomes
    excluded_genomes = ch_routed_genomes.excluded
    unresolved_genomes = ch_unresolved_genomes
    bacterial_report = PUBLISH_GENOME_ROUTING.out.bacterial
    excluded_report = PUBLISH_GENOME_ROUTING.out.excluded
    unresolved_report = PUBLISH_GENOME_ROUTING.out.unresolved
    versions = ch_versions
}
