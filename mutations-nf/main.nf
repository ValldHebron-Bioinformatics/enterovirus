#!/usr/bin/env nextflow
    
nextflow.enable.dsl = 2

include { INPUT_PREPARATION1; INPUT_PREPARATION2; FIND_MUTATIONS } from './modules/consensus-mutations'
include { VIRAL_POPULATION } from './modules/viral-population'


// Workflow
Channel
    .fromPath( params.file )
    .splitCsv( header: false, sep: ';' )
    .map { row -> tuple( row[0], row[1], row[2], row[3] ) }
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

