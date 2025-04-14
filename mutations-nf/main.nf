#!/usr/bin/env nextflow
    
nextflow.enable.dsl = 2

include { INPUT_PREPARATION1; INPUT_PREPARATION2; FIND_MUTATIONS } from './modules/consensus-mutations'
include { VIRAL_POPULATION } from './modules/viral-population'


// Workflow
Channel
    .fromPath( params.file )
    .splitCsv( header: true, sep: ',' )
    .map { row -> 
        tuple(file(row.SAMPLE_DIR),
        row.prot,
        file(row.VP1consensus),
        file(row.EVreference),
        val(row.genotype))
    }
    .set { inputs }

workflow MUTATIONS_CONSENSUS { 

    fstStep = INPUT_PREPARATION1(inputs)
    scdStep = INPUT_PREPARATION2(fstStep)
    FIND_MUTATIONS(inputs, scdStep)

}

workflow VIRAL_POPULATION_wf() {

    fstStep = INPUT_PREPARATION1(inputs)
    VIRAL_POPULATION(inputs, fstStep)

}

workflow {
    if (params.viral_population == "yes"){
        VIRAL_POPULATION_wf()
    } else {
        MUTATIONS_CONSENSUS()
    }
}
