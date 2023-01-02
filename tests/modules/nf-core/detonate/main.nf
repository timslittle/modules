#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { DETONATE } from '../../../../modules/nf-core/detonate/main.nf'

workflow test_detonate_pairedEnd {
    
    rnaseq = [
        [ id:'test', single_end:false ], // meta map
        [
            file(params.test_data['sarscov2']['illumina']['test_1_fastq_gz'], checkIfExists: true),
            file(params.test_data['sarscov2']['illumina']['test_2_fastq_gz'], checkIfExists: true)
        ]
    ]

    fasta = [
        [ id:'test', single_end:false ], // meta map
        [
            file(params.test_data['sarscov2']['genome']['transcriptome_fasta'], checkIfExists: true)
        ]
    ]

    DETONATE ( rnaseq, fasta )
}
