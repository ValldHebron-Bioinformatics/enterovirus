#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { CREATEDIR; QUALCONTROL; FILTHOST; TRIMPRIMERSR; TRIMPRIMERSL } from './modules/quality-control'
include { SPADES                                                                  } from './modules/assembly'
include { BLASTN; GETBLASTNMATCH; GETCDS; DIAMOND; GENOTYPEVP1                    } from './modules/genotyping'

// Checking user-defined parameters
if (params.protocol != "complete" && params.protocol != "partial") {
    exit 1, "Unknown protocol '${params.protocol}'. Choose any of [ complete | partial ]"
}

// Check user is defined
if (!params.user) {
    exit 1, "User is not defined. Please set the user parameter."
}

// Check if the input file is provided
if (!params.file) {
    exit 1, "Input file is not provided. Please set the file parameter."
}

// Check if the input file exists
if (!file(params.file).exists()) {
    exit 1, "Input file '${params.file}' does not exist."
}

// Check if the input file is a valid CSV
if (!params.file.endsWith('.csv')) {
    exit 1, "Input file '${params.file}' is not a valid CSV file. Please provide a CSV file."
}

// Check if the file type is defined
if (!params.fileType) {
    exit 1, "Input file type is not defined. Please set the fileType parameter."
}

// Check input file type is either fastq or fasta
if (params.fileType != 'fastq' && params.fileType != 'fasta') {
    exit 1, "Input file type '${params.fileType}' is not supported. Please choose either 'fastq' or 'fasta'."
}

// Check threads parameter
if (!params.threads) {
    exit 1, "Threads parameter is not defined. Please set the threads parameter."
}

// Workflow
Channel
    .fromPath(params.file)
    .splitCsv(header: false, sep: ';')
    .map { row -> tuple(row[0], row[1], row[2]) }
    .set { sample_run_ch }

workflow {
    // Create directories for all samples
    dir_ch = CREATEDIR(sample_run_ch)

    // Process FASTQ files
    if (params.fileType == 'fastq')
    {
        fastq_qc_ch = QUALCONTROL(dir_ch)
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
        processed_samples_ch = SPADES(spades_input_ch)

    }
    // Process FASTA files
    else if (params.fileType == 'fasta') {
        processed_samples_ch = dir_ch.map { tuple(it[0], it[1], it[2]) }
    }

    // Run final genotyping steps
    blastn_ch = BLASTN(processed_samples_ch)
    get_match_ch = GETBLASTNMATCH(blastn_ch)
    get_cds_ch = GETCDS(get_match_ch)
    prot_match_ch = DIAMOND(get_cds_ch)
    vp1_ch = GENOTYPEVP1(prot_match_ch)
}
