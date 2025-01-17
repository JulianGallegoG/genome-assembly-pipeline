#!/usr/bin/env nextflow

// Import required modules
// BAM processing module - will be used in production
include { CUTADAPT      } from "$projectDir/modules/nf-core/cutadapt/main"
include { FASTQC        } from "$projectDir/modules/nf-core/fastqc/main"

workflow PREPROCESS_READS {
    take:
    hifi_reads
    main:
    // Initialize versions channel
    ch_versions = Channel.empty()

    // Run adapter trimming with Cutadapt
    CUTADAPT(hifi_reads)
    ch_versions = ch_versions.mix(CUTADAPT.out.versions)

    /*
    // Perform quality control analysis with FastQC
    FASTQC(CUTADAPT.out.reads)
    ch_versions = ch_versions.mix(FASTQC.out.versions)
    */

    emit:
    reads_raw   = hifi_reads               // channel: [ meta, fastq ]
    reads_clean = CUTADAPT.out.reads      // channel: [ meta, trimmed_fastq ]
    //fastqc_html = FASTQC.out.html         // channel: [ meta, html ]
    //fastqc_zip  = FASTQC.out.zip          // channel: [ meta, zip ]
    versions    = ch_versions             // channel: [ versions.yml ]
}
