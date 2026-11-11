process VALIDATE_PARATYPE_BUNDLE {
    tag 'Paratype bundle'
    label 'process_low'

    input:
    path paratype_bundle
    path validator_script

    output:
    path 'paratype_bundle', emit: bundle
    path 'paratype_bundle.ready', emit: ready
    path 'versions.yml', emit: versions

    script:
    """
    python ${validator_script} --bundle "${paratype_bundle}"
    touch paratype_bundle.ready

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        validate_paratype_bundle.py: \$(python ${validator_script} --version)
        python: \$(python --version 2>&1 | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    """
    python ${validator_script} --bundle "${paratype_bundle}"
    touch paratype_bundle.ready

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        validate_paratype_bundle.py: \$(python ${validator_script} --version)
        python: \$(python --version 2>&1 | sed 's/Python //')
    END_VERSIONS
    """
}
