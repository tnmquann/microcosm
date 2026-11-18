#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// Authors: @tnmquann
// Maintainers: @tnmquann

/*
 * Run-level reporting:
 *  - ${outdir}/pipeline_info/microcosm_software_versions.yml (every tool of the run)
 *  - ${outdir}/MultiQC/<prefix>_multiqc_report.html (unless --skip_multiqc), styled by
 *    assets/multiqc_config.yml + assets/microcosm_logo.png, with a software-versions
 *    section and the run summary.
 * The MultiQC copy of the versions is written once every 'versions' topic producer has
 * finished, so MultiQC itself (patched to not publish on the topic) is added afterwards.
 */

include { MULTIQC } from '../../../modules/nf-core/multiqc/main'
include { versionRecords; versionsToYaml } from '../utils_software_versions/main'

workflow REPORTING {
    take:
    legacy_versions // channel: versions.yml files / version tuples from the subworkflows
    multiqc_files   // channel: files for MultiQC (QC reports, *_mqc.tsv, workflow summary)

    main:
    def records = legacy_versions
        .flatMap { item -> versionRecords(item) }
        .mix(channel.topic('versions').map { process_name, tool, version -> [process_name.toString(), tool.toString(), version.toString()] })
        .unique()

    def late_versions = channel.empty()
    def ch_report = channel.empty()
    if (!params.skip_multiqc) {
        def ch_mqc_versions = records
            .collect(flat: false)
            .ifEmpty([])
            .map { rows -> versionsToYaml(rows) }
            .collectFile(name: 'microcosm_mqc_versions.yml')

        def pipeline_root = params.pipeline_root ?: projectDir
        MULTIQC(
            multiqc_files
                .mix(ch_mqc_versions)
                .filter { mqc_file -> mqc_file }
                .collect()
                .map { mqc_files ->
                    [
                        [id: 'microcosm'],
                        mqc_files,
                        ["${pipeline_root}/assets/multiqc_config.yml", params.multiqc_config, params.multiqc_extra_config]
                            .findAll { config -> config }
                            .collect { config -> file(config, checkIfExists: true) },
                        file(params.multiqc_logo ?: "${pipeline_root}/assets/microcosm_logo.png", checkIfExists: true),
                        params.multiqc_replace_names ? file(params.multiqc_replace_names, checkIfExists: true) : [],
                        params.multiqc_sample_names ? file(params.multiqc_sample_names, checkIfExists: true) : []
                    ]
                }
        )
        late_versions = MULTIQC.out.versions_multiqc
        ch_report = MULTIQC.out.report
    }

    def ch_yaml = records
        .mix(late_versions.map { process_name, tool, version -> [process_name.toString(), tool.toString(), version.toString()] })
        .unique()
        .collect(flat: false)
        .ifEmpty([])
        .map { rows -> versionsToYaml(rows) }
        .collectFile(storeDir: "${params.outdir}/pipeline_info", name: 'microcosm_software_versions.yml')

    emit:
    versions_yaml = ch_yaml
    multiqc_report = ch_report
}
