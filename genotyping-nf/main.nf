#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { CREATEDIR; QUALCONTROL; FILTHOST; TRIMPRIMERSR; TRIMPRIMERSL; GETEVREADS } from './modules/quality-control'
include { SPADES; NREPLACER; ASSEMBLYMETRICS                                       } from './modules/assembly'
include { BLASTN; GETBLASTNMATCH; GETCDS; DIAMOND; GENOTYPEVP1                     } from './modules/genotyping'

workflow CHECK_INPUTS {
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

}

// Workflow
workflow {
    CHECK_INPUTS()

    Channel
        .fromPath(params.file)
        .splitCsv(header: false, sep: ',')
        .map { row -> tuple(row[0], row[1], row[2]) }
        .set { sample_run_ch }

    // Create directories for all samples
    dir_ch = CREATEDIR(sample_run_ch)

    // Process FASTQ files
    if (params.fileType == 'fastq') {
        fastq_qc_ch = QUALCONTROL(dir_ch)
        if (params.primers) {
            primers = Channel.of(params.primers)
            trimr_primers_ch = TRIMPRIMERSR(fastq_qc_ch.combine(primers))
            filthost_input_ch = TRIMPRIMERSL(trimr_primers_ch)
            spades_input_ch = FILTHOST(filthost_input_ch)
        }
        else {
            spades_input_ch = FILTHOST(fastq_qc_ch)
        }

        processed_samples_ch = SPADES(spades_input_ch)

    // Process FASTA files
    } else if (params.fileType == 'fasta') {
        processed_samples_ch = dir_ch.map { tuple(it[0], it[1], it[3]) }

    // Error message for unsupported file types
    } else {
        exit 1, "Input file type is wrongly defined. Please set the fileType parameter bewteen fasta or fastq."
    }

    // Run final genotyping steps
    blastn_ch = BLASTN(processed_samples_ch)
    get_match_ch = GETBLASTNMATCH(blastn_ch)
    get_cds_ch = GETCDS(get_match_ch)
    prot_match_ch = DIAMOND(get_cds_ch)
    vp1_ch = GENOTYPEVP1(prot_match_ch)
    if (params.fileType == "fastq") {
    	ev_ch = GETEVREADS(spades_input_ch, get_match_ch)
        nreplaced_ch = NREPLACER(spades_input_ch, get_match_ch)
        // assembly_ch = 
        ASSEMBLYMETRICS(nreplaced_ch)
    } 
    //else { // TO DO: Integrarlo con el formato de entrada FASTA.
    // se le debe dar de input el directorio donde este el análisis, ya que va a buscar
    // el archivo ${dirSample}/results/ev-match.fasta
        // assembly_ch = ASSEMBLYMETRICS([params.dirFasta])
    //}
}
