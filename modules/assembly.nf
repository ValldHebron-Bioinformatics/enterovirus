#!/usr/bin/env nextflow
    
nextflow.enable.dsl = 2


process SPADES {
    /**
    Ensamblado de novo con SPAdes
    */
    input:
    tuple val(sample), val(fastq1), val(fastq2), val(dirSample)

    output:
    tuple env('dirFASTA'), env('seqsFasta')

    script:
    if (params.protocol == 'complete')
        """
        #!/bin/bash
        spades.py --rnaviral --threads 24 -1 $fastq1 -2 $fastq2 -o ${dirSample}/assembly/spades
        if [ -f ${dirSample}/assembly/spades/scaffolds.fasta ]; then
            seqsFasta=${dirSample}/assembly/spades/scaffolds.fasta
        elif [ -f ${dirSample}/assembly/spades/contigs.fasta ]; then
            seqsFasta=${dirSample}/assembly/spades/contigs.fasta
        else
            echo "No spades contigs or scaffold files found." > ${dirSample}/assembly/spades/blastn.log
        fi
        dirFASTA=${dirSample}/assembly
        """
    else if (params.protocol == 'partial')
        """
        #!/bin/bash
        spades.py --isolate --threads 24 -1 $fastq1 -2 $fastq2 -o ${dirSample}/assembly/spades
        if [ -f ${dirSample}/assembly/spades/scaffolds.fasta ]; then
            seqsFasta=${dirSample}/assembly/spades/scaffolds.fasta
        elif [ -f ${dirSample}/assembly/spades/contigs.fasta ]; then
            seqsFasta=${dirSample}/assembly/spades/contigs.fasta
        else
            echo "No spades contigs or scaffold files found." > ${dirSample}/assembly/spades/blastn.log
        fi
        dirFASTA=${dirSample}/assembly
        """
}
