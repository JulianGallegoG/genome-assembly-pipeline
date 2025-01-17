#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// Import subworkflow
include { PREPARE_INPUT } from "$projectDir/subworkflows/local/prepare_input/main"
include { PREPROCESS_READS } from "$projectDir/subworkflows/local/preprocess_reads/main"
include { KMER_ANALYSIS   } from "$projectDir/subworkflows/local/kmer_analysis/main"

workflow {
    // Define constants
    def workflow_permitted_stages = [
        'prepare_input',
        'preprocess_reads',
        'kmer_analysis'
        ]

     // Check input
    def workflow_steps = params.steps.tokenize(",")
    if ( ! workflow_steps.every { it in workflow_permitted_stages } ) {
        error "Unrecognised workflow step in $params.steps ( $workflow_permitted_stages )"
    }

    // Print pipeline information
    log.info """
    Running FishEvoLab Assembly Workflow
    ===============================
    input     : ${params.input}
    outdir    : ${params.outdir}
    """
    .stripIndent()

    // Setup sink channels
    ch_versions = Channel.empty()

    // Read in data
    if ('prepare_input' in workflow_steps) {
        PREPARE_INPUT (params.input)
    }

    // Run preprocessing if included in steps
    if ('preprocess_reads' in workflow_steps) {
        PREPROCESS_READS(PREPARE_INPUT.out.hifi)
        ch_versions = ch_versions.mix(PREPROCESS_READS.out.versions)
    }

    // Run k-mer analysis if included in steps
    if ('kmer_analysis' in workflow_steps) {
        KMER_ANALYSIS(PREPARE_INPUT.out.hifi)
        ch_versions = ch_versions.mix(KMER_ANALYSIS.out.versions)
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
