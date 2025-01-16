#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// Import subworkflow
//include { PREPROCESS_READS } from './subworkflows/local/preprocess_reads/main'

// Define permitted workflow stages (similar to EBP pipeline)
//def workflow_permitted_stages = ['preprocess']  // We can expand this later

def workflow_permitted_stages = ['preprocess']

workflow {
    // Check input steps
    def workflow_steps = params.steps.tokenize(",")
    if (!workflow_steps.every { it in workflow_permitted_stages }) {
        error "Unrecognised workflow step in $params.steps ( $workflow_permitted_stages )"
    }

    // Print pipeline information
    log.info """
    HiFi READ PREPROCESSING PIPELINE
    ===============================
    input     : ${params.input}
    outdir    : ${params.outdir}
    """
    .stripIndent()

    // Setup sink channels
    ch_versions = Channel.empty()

    // Create channel for input YAML and validate it exists
    input_data = Channel.fromPath(params.input)
        .map { yaml_file ->
            def meta = [id: yaml_file.simpleName]
            [meta, yaml_file]
        }

    // Run preprocessing if included in steps
    if ('preprocess' in workflow_steps) {
        PREPROCESS_READS(input_data)
        ch_versions = ch_versions.mix(PREPROCESS_READS.out.versions)
    }
}

// Workflow completion handler
workflow.onComplete {
    def msg = """\
        Pipeline execution summary
        ---------------------------
        Completed at: ${workflow.complete}
        Duration    : ${workflow.duration}
        Success     : ${workflow.success}
        workDir     : ${workflow.workDir}
        exit status : ${workflow.exitStatus}
        """
        .stripIndent()

    log.info(msg)
}
