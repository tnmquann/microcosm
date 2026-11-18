#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// Authors: @tnmquann
// Maintainers: @tnmquann

/*
 * Helpers for the run-wide software versions YAML (used by subworkflows/local/reporting):
 * legacy versions.yml files and (process, tool, version) tuples from the 'versions' topic.
 */

def versionRecordsFromYaml(yaml_file) {
    def data = new org.yaml.snakeyaml.Yaml().load(file(yaml_file).text) ?: [:]
    return data.collectMany { process_name, tools ->
        def tool_map = (tools instanceof Map) ? tools : [:]
        tool_map.collect { tool, version ->
            [process_name.toString(), tool.toString(), version.toString()]
        }
    }
}

// Items may be versions.yml paths, (process, tool, version) tuples, or lists of either.
def versionRecords(item) {
    if (item == null) {
        return []
    }
    if (item instanceof java.nio.file.Path) {
        return versionRecordsFromYaml(item)
    }
    if (item instanceof List) {
        def is_tuple = item.size() == 3 && item.every { value -> !(value instanceof java.nio.file.Path) && !(value instanceof List) }
        if (is_tuple) {
            return [item.collect { value -> value.toString() }]
        }
        return item.collectMany { value -> versionRecords(value) }
    }
    return []
}

def versionsToYaml(records) {
    def by_process = new TreeMap()
    records.each { process_name, tool, version ->
        def short_name = process_name.tokenize(':').last()
        by_process.computeIfAbsent(short_name) { _key -> new TreeMap() }.put(tool, version)
    }
    by_process.put('Workflow', [
        (workflow.manifest.name ?: 'microcosm'): (workflow.manifest.version ?: ''),
        Nextflow: workflow.nextflow.version.toString()
    ])
    def options = new org.yaml.snakeyaml.DumperOptions()
    options.setDefaultFlowStyle(org.yaml.snakeyaml.DumperOptions.FlowStyle.BLOCK)
    return new org.yaml.snakeyaml.Yaml(options).dump(by_process)
}
