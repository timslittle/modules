#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { TRINITY } from '../../../../modules/nf-core/trinity/main.nf'

workflow test_trinity_paired_end {
    
    input = [
        [ id:'test', single_end:false ], // meta map
        [ file(params.test_data['homo_sapiens']['illumina']['fastq']['test_rnaseq_1_fastq_gz'], checkIfExists: true),
                file(params.test_data['homo_sapiens']['illumina']['fastq']['test_rnaseq_2_fastq_gz'], checkIfExists: true) ]
    ]

    TRINITY ( input )
}