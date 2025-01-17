#!/usr/bin/env nextflow

// Import required modules
include { MERYL_COUNT    } from "$projectDir/modules/nf-core/meryl/count/main"
include { MERYL_UNIONSUM } from "$projectDir/modules/nf-core/meryl/unionsum/main"
include { MERYL_HISTOGRAM } from "$projectDir/modules/nf-core/meryl/histogram/main"
include { GENOMESCOPE2   } from "$projectDir/modules/nf-core/genomescope2/main"

workflow KMER_ANALYSIS {
    take:
    hifi_reads

    main:
    // Create k-mer database and histogram with meryl
    MERYL_COUNT(hifi_reads, params.meryl.kmer_size)
    MERYL_UNIONSUM(MERYL_COUNT.out.meryl_db.transpose().groupTuple(), params.meryl.kmer_size)
    MERYL_HISTOGRAM(MERYL_UNIONSUM.out.meryl_db, params.meryl.kmer_size)

    // Run GenomeScope2 analysis
    GENOMESCOPE2(MERYL_HISTOGRAM.out.hist)

    MERYL_COUNT.out.versions.first().mix(
        MERYL_UNIONSUM.out.versions.first(),
        MERYL_HISTOGRAM.out.versions.first(),
        GENOMESCOPE2.out.versions.first()
    ).set { versions }

    emit:
    meryldb    = MERYL_COUNT.out.meryl_db      // channel: [ meta, meryldb ]
    histogram  = MERYL_HISTOGRAM.out.hist // channel: [ meta, hist ]
    summary    = GENOMESCOPE2.out.summary     // channel: [ meta, summary ]
    versions
}
