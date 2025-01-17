#!/usr/bin/env nextflow

// Import required YAML parsing functionality
import org.yaml.snakeyaml.Yaml

workflow PREPARE_INPUT {
    take:
    infile // The input YAML file path

    main:
    // Read in YAML and create initial channel
    ch_input = Channel.fromPath(infile, checkIfExists: true)
        .map { file -> readYAML(file) }
        .map { yaml -> yaml + [id: yaml.sample.name.replace(" ","_")] }

    // Process input data and create channels for HiFi and Hi-C
    ch_input
        .multiMap { data ->
            hic_ch: (data.hic ? [
                data.subMap('id', 'sample') + [single_end: false],
                data.hic.collect { [
                    file(it.read1, checkIfExists: true),
                    file(it.read2, checkIfExists: true)
                ]}
            ] : [])
            hifi_ch: (data.hifi ? [
                data.subMap('id', 'sample') + [single_end: true],
                data.hifi.collect { file(it.reads, checkIfExists: true) }
            ] : [])
        }
        .set { input }

    // Process Hi-C data
    input.hic_ch
        .filter { !it.isEmpty() }
        .flatMap { meta, hic_pairs ->
            hic_pairs.withIndex().collect { pair, index ->
                [meta + [pair_id: index], pair]
            }
        }
        .set { hic_fastx_ch }

    // Process HiFi data
    input.hifi_ch
        .filter { !it.isEmpty() }
        .transpose()   // Transform to [ [ id: 'sample_name'], file('/path/to/read') ]
        .set { hifi_fastx_ch }

    emit:
    hic = hic_fastx_ch.dump(tag: 'Input: Hi-C')
    hifi = hifi_fastx_ch.dump(tag: 'Input: PacBio HiFi')
}

// Helper function to read YAML files
def readYAML(yamlfile) {
    return new Yaml().load(new FileReader(yamlfile.toString()))
}
