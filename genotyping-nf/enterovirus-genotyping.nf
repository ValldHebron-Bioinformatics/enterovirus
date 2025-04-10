#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { CREATEDIR; GETFASTQS; QUALCONTROL; FILTHOST; TRIMPRIMERSR; TRIMPRIMERSL } from './modules/quality-control'
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

// Workflow
Channel
    .fromPath(params.file)
    .splitCsv(header: false, sep: ';')
    .map { row -> tuple(row[0], row[1], row[2]) }
    .set { sample_run_ch }

workflow {
    // Create directories for all samples
    dir_ch = CREATEDIR(sample_run_ch.map { row -> row[0] })

    // Determine the file type and process samples
    sample_run_ch.map { sample_id, file1, file2 ->
        def extension = file1.tokenize('.').last()
        tuple(sample_id, file1, file2, extension)
    }.set { sample_run_with_ext_ch }

    sample_run_with_ext_ch.map { sample_id, file1, file2, extension ->
        if (extension == "fastq") {
            def seqsFasta = processFastq([sample_id, file1, file2])
            return tuple(sample_id, seqsFasta, extension)
        } else if (extension == "fasta") {
            def seqsFasta = Channel.of(file1)
            return tuple(sample_id, seqsFasta, extension)
        } else {
            exit 1, "Unsupported file extension '${extension}'. Only 'fastq' and 'fasta' are supported."
        }
    }.set { processed_samples_ch }

    // Ensure BLASTN waits for seqsFasta to be generated
    blastin_ch = processed_samples_ch.map { seqsFasta, extension -> [seqsFasta, extension, dir_ch.outputDir] }
    blastn_ch = BLASTN(blastin_ch)
    get_match_ch = GETBLASTNMATCH(blastn_ch)
    get_cds_ch = GETCDS(get_match_ch)
    prot_match_ch = DIAMOND(get_cds_ch)
    vp1_ch = GENOTYPEVP1(prot_match_ch)
}

// Helper function to process FASTQ files
def processFastq(input) {
    dir_ch = CREATEDIR(sample_run_ch)
    fastq_ch = GETFASTQS(dir_ch)
    fastq_qc_ch = QUALCONTROL(fastq_ch)
    filt_host_ch = FILTHOST(fastq_qc_ch)

    if (params.primers) {
        primers = Channel.of(params.primers)
        trimr_primers_ch = TRIMPRIMERSR(filt_host_ch.combine(primers))
        triml_primers_ch = TRIMPRIMERSL(trimr_primers_ch)
        spades_input_ch = triml_primers_ch
    } else {
        spades_input_ch = filt_host_ch
    }

    // Pass the required inputs to the SPADES module
    return SPADES(spades_input_ch.map { tuple(it[0], it[1], it[2], it[0]) })
}
