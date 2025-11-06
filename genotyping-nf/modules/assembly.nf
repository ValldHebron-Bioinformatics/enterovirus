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
    tuple val(sampleId), env('seqsFasta'), val(outputDir)

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


process NREPLACER {
    /**
    Reemplazo por N en las posiciones con menos de 20x
    */
    errorStrategy 'retry'
    maxRetries 2 

    input:
    tuple val(sampleId), val(fastq1), val(fastq2), val(outputDir)
    val(dirFASTA)

    output:
    val(outputDir)

    script:
    """
    #!/bin/bash
    minimap2 -ax sr ${outputDir}/results/ev-match.fasta $fastq1 $fastq2 > ${outputDir}/tmp/reads-mapped-consensus.sam
    samtools view -bS ${outputDir}/tmp/reads-mapped-consensus.sam | samtools sort - -o ${outputDir}/tmp/reads-mapped-consensus.bam
    echo -e "ref;pos;depth" > ${outputDir}/tmp/depth-consensus.tsv
    samtools depth -a ${outputDir}/tmp/reads-mapped-consensus.bam >> ${outputDir}/tmp/depth-consensus.tsv
    sed -i -z 's/\t/;/g' ${outputDir}/tmp/depth-consensus.tsv
    if [[ \$(wc -l < ${outputDir}/tmp/depth-consensus.tsv) -lt 3 ]]; then
        sleep 5
        exit 1  
    fi
    python3 $params.programs.Nreplacer --fasta ${outputDir}/results/ev-match.fasta --coverage ${outputDir}/tmp/depth-consensus.tsv
    EVmatch=${outputDir}/results/ev-match.fasta
    """
}

process ASSEMBLYMETRICS {
    /**
    obtener metricas de ensamblado
    */
    
    errorStrategy 'ignore'

    errorStrategy 'ignore'

    input:
    val(outputDir)

    script:
    """
    #!/bin/bash
    touch ${outputDir}/results/assembly-metrics.csv
    python3 $params.programs.qualitymetrics --fasta ${outputDir}/results/ev-match.fasta --csv ${outputDir}/results/assembly-metrics.csv --coverage ${outputDir}/tmp/depth-consensus.tsv
    """
}
