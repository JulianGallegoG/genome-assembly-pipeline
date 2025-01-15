#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// Pipeline parameter declarations
params.reads = null        // Path to input reads (e.g., '*_R{1,2}.fastq.gz')
params.outdir = 'results' // Output directory path

// Load modules and subworkflows
include { FASTQC } from './modules/nf-core/fastqc/main'

// Log pipeline info
log.info """
         Genome Assembly Pipeline
         ===================================
         input reads    : ${params.reads}
         output dir     : ${params.outdir}
         """
         .stripIndent()

// Main workflow
workflow {
    // Input validation
    if (!params.reads) {
        error "Please provide input reads with --reads"
    }

    // Create input channel for reads
    Channel
        .fromFilePairs(params.reads, checkIfExists: true)
        .ifEmpty { error "Cannot find any reads matching: ${params.reads}" }
        .set { input_reads }

    // Quality control
    FASTQC(input_reads)
}

// Completion handler
workflow.onComplete {
    log.info """
    Pipeline execution summary
    ---------------------------
    Completed at : ${workflow.complete}
    Duration    : ${workflow.duration}
    Success     : ${workflow.success}
    workDir     : ${workflow.workDir}
    exit status : ${workflow.exitStatus}
    """
}
