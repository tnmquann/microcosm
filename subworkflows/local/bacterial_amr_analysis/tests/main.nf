#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { BACTERIAL_AMR_ANALYSIS } from '../main'

// -stub-run checks wiring (database setup stubs, no downloads); --genomes_dir
// runs the real tools on <name>.fasta. E. coli and K. pneumoniae get an
// AMRFinderPlus --organism; Salmonella CT18 gets PointFinder at genus level.
workflow {
    def common = [
        run_id: 'AMR_TEST', group: 'AMR_TEST', assembler: 'external',
        lineage_source: 'provided', source_route: 'standalone',
        classification_status: 'provided', domain: 'd__Bacteria'
    ]
    def ent = 'd__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae'
    def cases = [
        ['ecoli_K12', 'g__Escherichia', 's__Escherichia coli'],
        ['kpneumoniae_HS11286', 'g__Klebsiella', 's__Klebsiella pneumoniae'],
        ['styphi_CT18', 'g__Salmonella', 's__Salmonella enterica']
    ]
    ch_genomes = channel.fromList(cases).map { name, genus, species ->
        def genome = params.genomes_dir
            ? file("${params.genomes_dir}/${name}.fasta", checkIfExists: true)
            : file("${projectDir}/../../../../test_data/genome_input/bacterial.fa", checkIfExists: true)
        tuple(
            common + [id: name, sample_id: name, genome_id: name, genus: genus, species: species,
                      lineage: "${ent};${genus};${species}"],
            genome
        )
    }

    BACTERIAL_AMR_ANALYSIS(ch_genomes)
    BACTERIAL_AMR_ANALYSIS.out.hamronized.view { meta, tsv ->
        "AMR_HAMRONIZED ${meta.genome_id} ${tsv.name} rows=${Math.max(tsv.readLines().size() - 1, 0)}"
    }
    BACTERIAL_AMR_ANALYSIS.out.hamronization_summary.view { tsv -> "AMR_SUMMARY ${tsv.name} rows=${Math.max(tsv.readLines().size() - 1, 0)}" }
    BACTERIAL_AMR_ANALYSIS.out.normalized.view { meta, tsv -> "AMR_ARGNORM ${meta.genome_id} ${tsv.name}" }
    BACTERIAL_AMR_ANALYSIS.out.normalized_summary.view { _meta, tsv -> "AMR_ARGNORM_SUMMARY ${tsv.name}" }
}
