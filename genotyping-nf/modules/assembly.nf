#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process SPADES {
    /**
    De novo assembly with SPAdes.
    This process dynamically adjusts SPAdes parameters based on the protocol (complete or partial).
    */

    errorStrategy 'ignore'

    input:
    tuple val(sample), val(fastq1), val(fastq2), val(userDir)

    output:
    path "${userDir}/${sample}/assembly/spades/scaffolds.fasta", optional: true
    path "${userDir}/${sample}/assembly/spades/contigs.fasta", optional: true
    path "${userDir}/${sample}/errors.log", emit: errors

    script:
    // Determine SPAdes mode based on the protocol
    def spadesMode = params.protocol == 'complete' ? '--rnaviral' : '--isolate'

    """
    #!/bin/bash
    spades.py $spadesMode --threads ${params.threads ?: 24} -1 $fastq1 -2 $fastq2 -o ${userDir}/assembly/spades

    # Check for output files and set the seqsFasta variable
    if [ -f ${userDir}/${sample}/assembly/spades/scaffolds.fasta ]; then
        seqsFasta=${userDir}/${sample}/assembly/spades/scaffolds.fasta
    elif [ -f ${userDir}/assembly/spades/contigs.fasta ]; then
        seqsFasta=${userDir}/${sample}/assembly/spades/contigs.fasta
    else
        echo "No SPAdes contigs or scaffolds assembled." >> ${userDir}/${sample}/errors.log
        exit 1
    fi
    """
}
