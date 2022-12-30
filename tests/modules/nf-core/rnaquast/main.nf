#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { RNAQUAST } from '../../../../modules/nf-core/rnaquast/main.nf'

workflow test_rnaquast {
    
    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_paired_end_bam'], checkIfExists: true)
    ]

    RNAQUAST ( input )
}
