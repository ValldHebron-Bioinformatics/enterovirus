#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process SPADES {
    /**
    De novo assembly with SPAdes.
    This process dynamically adjusts SPAdes parameters based on the protocol (complete or partial).
    */

    errorStrategy 'ignore'

    input:
    tuple val(sampleId), val(fastq1), val(fastq2), val(outputDir)

    output:
    tuple val(sample_id), env('seqsFasta'), val(outputDir)

    script:
    // Determine SPAdes mode based on the protocol
    def spadesMode = params.protocol == 'complete' ? '--rnaviral' : '--isolate'

    """
    #!/bin/bash
    spades.py $spadesMode --threads $params.threads -1 $fastq1 -2 $fastq2 -o ${outputDir}/assembly/spades

    # Check for output files and set the seqsFasta variable
    if [ -f ${outputDir}/assembly/spades/scaffolds.fasta ]; then
        seqsFasta=${outputDir}/assembly/spades/scaffolds.fasta
    elif [ -f ${outputDir}/assembly/spades/contigs.fasta ]; then
        seqsFasta=${outputDir}/assembly/spades/contigs.fasta
    else
        echo "No SPAdes contigs or scaffolds assembled." >> ${outputDir}/errors.log
        exit 1
    fi
    """
}
