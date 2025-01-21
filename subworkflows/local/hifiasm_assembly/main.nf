#!/usr/bin/env nextflow

// Import required modules with distinct names for each GFASTATS instance
include { HIFIASM                        } from "$projectDir/modules/nf-core/hifiasm/main"
include { GFASTATS as GFASTATS_PRIMARY   } from "$projectDir/modules/nf-core/gfastats/main"
include { GFASTATS as GFASTATS_PATERNAL  } from "$projectDir/modules/nf-core/gfastats/main"
include { GFASTATS as GFASTATS_MATERNAL  } from "$projectDir/modules/nf-core/gfastats/main"

workflow HIFIASM_ASSEMBLY {
    take:
    hifi_reads   // channel: [ meta, reads ] from prepare_input
    hic_reads    // channel: [ meta, [read1, read2] ] from prepare_input

    main:
    // Group HiFi reads by sample for assembly
    ch_hifi_grouped = hifi_reads
        .groupTuple()
        .map { meta, reads -> [meta, reads[0]] }

    // Process Hi-C reads to match HIFIASM module expectations
    ch_hic_processed = hic_reads
        .map { meta, reads ->
            def new_meta = meta.findAll { key, value -> key != 'pair_id' }
            [new_meta, reads[0], reads[1]]
        }
        .groupTuple(by: [0])
        .map { meta, r1s, r2s -> [meta, r1s[0], r2s[0]] }

    // Run HIFIASM with Hi-C integration
    HIFIASM(
        ch_hifi_grouped,                   // HiFi reads
        Channel.value([[:], [], []]),      // Empty trio binning channel
        ch_hic_processed                   // Hi-C reads
    )

    // Convert primary contigs to FASTA
    GFASTATS_PRIMARY(
        HIFIASM.out.primary_contigs,
        "fasta",              // Output format
        [],                   // No genome size needed
        "",                   // No target
        [],                   // No AGP file
        [],                   // No include bed
        [],                   // No exclude bed
        []                    // No SAK instructions
    )

    // Convert paternal contigs to FASTA
    GFASTATS_PATERNAL(
        HIFIASM.out.paternal_contigs,
        "fasta",
        [], "", [], [], [], []
    )

    // Convert maternal contigs to FASTA
    GFASTATS_MATERNAL(
        HIFIASM.out.maternal_contigs,
        "fasta",
        [], "", [], [], [], []
    )

    // Collect all version information
    versions = HIFIASM.out.versions.first()
        .mix(GFASTATS_PRIMARY.out.versions.first())
        .mix(GFASTATS_PATERNAL.out.versions.first())
        .mix(GFASTATS_MATERNAL.out.versions.first())

    emit:
    // GFA format outputs from HIFIASM
    gfa_primary   = HIFIASM.out.primary_contigs     // channel: [ meta, gfa ]
    gfa_paternal  = HIFIASM.out.paternal_contigs    // channel: [ meta, gfa ]
    gfa_maternal  = HIFIASM.out.maternal_contigs    // channel: [ meta, gfa ]

    // FASTA format outputs from GFASTATS
    fasta_primary   = GFASTATS_PRIMARY.out.assembly    // channel: [ meta, fasta ]
    fasta_paternal  = GFASTATS_PATERNAL.out.assembly   // channel: [ meta, fasta ]
    fasta_maternal  = GFASTATS_MATERNAL.out.assembly   // channel: [ meta, fasta ]

    // Additional outputs
    logs     = HIFIASM.out.log      // channel: [ meta, log ]
    versions                        // channel: [ versions.yml ]
}
