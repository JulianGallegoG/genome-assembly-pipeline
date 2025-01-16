#!/usr/bin/env nextflow

import org.yaml.snakeyaml.Yaml

// Import required modules
// BAM processing module - will be used in production
include { SAMTOOLS_BAM2FQ } from "$projectDir/modules/nf-core/samtools/bam2fq/main"
include { CUTADAPT      } from "$projectDir/modules/nf-core/cutadapt/main"
include { FASTQC        } from "$projectDir/modules/nf-core/fastqc/main"

// Function to read and parse YAML file
def parseYamlFile(yamlPath) {
    def yaml = new Yaml()
    return yaml.load(new File(yamlPath).text)
}

workflow PREPROCESS_READS {
    take:
    input_data    // channel: [ meta, yaml_file ]

    main:
    ch_versions = Channel.empty()

    // Parse YAML and create read channels
    // The map operator transforms each input tuple into HiFi read entries
    reads_ch = input_data
        .map { meta, yaml_file ->
            // Parse the YAML file to extract sample information and read paths
            def data = parseYamlFile(yaml_file.toString())

            // Transform the HiFi entries into a list of meta+read pairs
            def hifiReads = data.hifi.collect { entry ->
                // Create a metadata map that includes sample information
                def read_meta = meta + [
                    sample: data.sample.name.replaceAll(/\s+/, '_'),
                    single_end: true
                ]
                // Return a tuple of metadata and read file path
                [read_meta, file(entry.reads)]
            }
            return hifiReads
        }
        .flatten()
        .buffer(size: 2)

    // Sort reads into appropriate channels based on file type
    reads_ch
        .branch { meta, reads ->
            fasta: reads.toString() =~ /\.fa(sta)?(.gz)?$/
                return [ meta, reads ]
            fastq: reads.toString() =~ /\.f(ast)?q(.gz)?$/
                return [ meta, reads ]
            fail: true
                return [ meta, reads ]
        }
        .set { ch_input_reads }

    // Provide clear error messages for unrecognized file types
    ch_input_reads.fail
        .view { meta, reads ->
            error "ERROR: Unrecognized file format for ${meta.sample}: ${reads}"
        }

    // Combine FASTA and FASTQ inputs into a single channel for processing
    ch_input_reads.fastq
        .mix(ch_input_reads.fasta)
        .set { ch_fastq }

    // Log which files are being processed
    ch_fastq
        .view { meta, reads ->
            log.info "Processing HiFi reads for ${meta.sample}: ${reads}"
            return [meta, reads]
        }

    // Run adapter trimming with Cutadapt
    CUTADAPT(
        ch_fastq
    )
    ch_versions = ch_versions.mix(CUTADAPT.out.versions)

    // Perform quality control analysis with FastQC
    FASTQC(
        CUTADAPT.out.reads
    )
    ch_versions = ch_versions.mix(FASTQC.out.versions)

    emit:
    reads_raw   = ch_fastq                // channel: [ meta, fastq ]
    reads_clean = CUTADAPT.out.reads      // channel: [ meta, trimmed_fastq ]
    fastqc_html = FASTQC.out.html         // channel: [ meta, html ]
    fastqc_zip  = FASTQC.out.zip          // channel: [ meta, zip ]
    versions    = ch_versions             // channel: [ versions.yml ]
}
