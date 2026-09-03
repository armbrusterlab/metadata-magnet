#!/usr/bin/env nextflow
include { getMetadata } from './modules/getMetadata.nf'
include { checkTaxidCoverage } from './modules/checkCoverage.nf'

params {
    // always use params file to pass arguments; otherwise ints, floats, and Booleans may be interpreted as strings
    sra_query: String
    genome_length: Float
    taxid: Integer
    coverage_threshold: Integer
}

workflow {
    main:


    publish:

}


output {
    // SAVE METADATA CSV FILES TO metadata/ DIR BECAUSE THAT'S THE ASSUMPTION IN NEXTFLOW.CONFIG

}