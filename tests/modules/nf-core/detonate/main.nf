#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { DETONATE } from '../../../../modules/nf-core/detonate/main.nf'

workflow test_detonate {
    
    input = [
        [ id:'test', single_end:false ], // meta map
        file(params.test_data['sarscov2']['illumina']['test_paired_end_bam'], checkIfExists: true)
    ]

    DETONATE ( input )
}
