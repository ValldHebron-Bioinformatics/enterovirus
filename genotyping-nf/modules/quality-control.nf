#!/usr/bin/env nextflow
    
nextflow.enable.dsl = 2

process CREATEDIR {
    /**
    Create directories and set variables
    */

    input:
    tuple val(sample_id), val(file1), val(file2)

    output:
    tuple val(sample_id), val(file1), val(file2), env('outputDir'), env('extension')

    script:
    """
    #!/bin/bash
    outputDir=$params.workdir/$params.user/$sample_id/
    mkdir -p \$outputDir
    if [[ $file1 == *.fastq.gz ]]; then
        extension="fastq.gz"
    elif [[ $file1 == *.fastq ]]; then
        extension="fastq"
    elif [[ $file1 == *.fasta ]]; then
        extension="fasta"
    else
        echo "Unsupported file extension for $file1"
        exit 1
    fi
    """
}

process QUALCONTROL {
    /**
    Filtrado con trimmomatic
    */
    
    input:
    tuple val(sample_id), val(fastq1), val(fastq2), val(outputDir), val(extension)

    output:
    tuple val(sample_id), path('*_paired-trim_1.fastq.gz'), path('*_paired-trim_2.fastq.gz'), val(outputDir), val(extension)
    
    script:
    def output = "${sample_id}_paired-trim_1.fastq.gz ${sample_id}_unpaired-trim_1.fastq.gz ${sample_id}_paired-trim_2.fastq.gz ${sample_id}_unpaired-trim_2.fastq.gz"
    """
    #!/bin/bash
    trimmomatic PE -threads $params.threads $fastq1 $fastq2 $output LEADING:30 TRAILING:30 SLIDINGWINDOW:10:30
    """
}

process FILTHOST {
    /**
    Filtrado de host reads
    */
    
    input:
    tuple val(sample_id), path(fastq1), path(fastq2), val(outputDir), val(extension)

    output:
    tuple val(sample_id), path('*_host_removed_R1.fastq.gz'), path('*_host_removed_R2.fastq.gz'), val(outputDir), val(extension)
    """
    #!/bin/bash  
    bowtie2 -p $params.threads -x $params.references.refHuman -1 $fastq1 -2 $fastq2 --un-conc-gz ${sample_id}_host_removed > ${sample_id}_host_mapped_and_unmapped.sam
    mv ${sample_id}_host_removed.1 ${sample_id}_host_removed_R1.fastq.gz; mv ${sample_id}_host_removed.2 ${sample_id}_host_removed_R2.fastq.gz
    """
}

process TRIMPRIMERSR {
    /**
    Filtrado de primers
    */
    
    input:
    tuple val(sample), path(fastq1), path(fastq2), val(dirSample), path(primers)
    
    output:
    tuple val(sample), path('*R1_clean_r.fastq.gz'), path('*R2_clean_r.fastq.gz'), val(dirSample), path(primers)
    """
    #!/bin/bash  
    $params.programs.bbduk in1=$fastq1 in2=$fastq2 out1=${sample}_R1_clean_r.fastq.gz out2=${sample}_R2_clean_r.fastq.gz ref=$primers k=15 ktrim=r restrictright=30
    """
}

process TRIMPRIMERSL {
    /**
    Filtrado de primers
    */
    publishDir "$dirSample/fastq", pattern: '*_clean.fastq.gz', mode: 'copy'
    
    input:
    tuple val(sample), path(fastq1), path(fastq2), val(dirSample), path(primers)
    
    output:
    tuple val(sample), path('*R1_clean.fastq.gz'), path('*R2_clean.fastq.gz'), val(dirSample)
    """
    #!/bin/bash  
    $params.programs.bbduk in1=${sample}_R1_clean_r.fastq.gz in2=${sample}_R2_clean_r.fastq.gz out1=${sample}_R1_clean.fastq.gz out2=${sample}_R2_clean.fastq.gz ref=$primers k=15 ktrim=l restrictleft=30
    """
}


