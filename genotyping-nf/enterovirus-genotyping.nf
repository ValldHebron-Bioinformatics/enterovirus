#!/usr/bin/env nextflow
    
nextflow.enable.dsl = 2

include { CREATEDIR; GETFASTQS; QUALCONTROL; FILTHOST; TRIMPRIMERSR; TRIMPRIMERSL } from './modules/quality-control'
include { SPADES                                                                  } from './modules/assembly'
include { BLASTN; GETBLASTNMATCH; GETCDS; DIAMOND; GENOTYPEVP1                    } from './modules/genotyping'

//Checking user-defined parameters  
if (params.protocol != "complete" && params.protocol != "partial") {
    exit 1, "Unknown protocol. Choose any of [ complete | partial ]"
}   

// Workflow
Channel
    .fromPath( params.file )
    .splitCsv( header: false, sep: ';' )
    .map { row -> tuple( row[0], row[1], row[2] ) }
    .set { sample_run_ch }


workflow { 
    if (params.input == "fastq") {
        dir_ch = CREATEDIR(sample_run_ch)
        fastq_ch = GETFASTQS(dir_ch)
        fastq_qc_ch = QUALCONTROL(fastq_ch)
        filt_host_ch = FILTHOST(fastq_qc_ch)
        if (params.primers) {
            primers = Channel.of(params.primers)
            trimr_primers_ch = TRIMPRIMERSR(filt_host_ch.combine(primers))
            triml_primers_ch = TRIMPRIMERSL(trimr_primers_ch)
            spades_input_ch = triml_primers_ch
        }
        else {
            spades_input_ch = filt_host_ch
        }
        seqsFasta = SPADES(spades_input_ch)
    }
    else {
        // Asumes fasta
        seqsFasta = params.file
        
    }
    blastin_ch = Channel.of([params.user, seqsFasta])
    blastn_ch = BLASTN(blastin_ch)
    get_match_ch = GETBLASTNMATCH(blastn_ch)
    get_cds_ch = GETCDS(get_match_ch)
    prot_match_ch = DIAMOND(get_cds_ch)
    vp1_ch = GENOTYPEVP1(prot_match_ch)
}

