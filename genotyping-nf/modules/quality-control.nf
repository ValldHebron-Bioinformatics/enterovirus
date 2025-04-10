#!/usr/bin/env nextflow
    
nextflow.enable.dsl = 2

process CREATEDIR {
    /**
    Create directories in samples folder from tmp files
    */

    input:
    val sample
    
    output:
    val outputDir
    
    script:
    """
    #!/bin/bash
    outputDir=$params.workdir/$params.user/$sample/
    mkdir -p $outputDir
    echo $outputDir
    """    
}

process GETFASTQS { 
    /**
    Get fastqs
    */
    
    input:
    tuple val(sample), val(fastq1), val(fastq2), val(dirSample)
    
    output:
    tuple val(sample), env('FASTQ1'), env('FASTQ2'), val(dirSample)
    
    script:
    """
    #!/bin/bash
    cp $params.rawfastqDir/$fastq1 $dirSample/fastq/$fastq1
    cp $params.rawfastqDir/$fastq2 $dirSample/fastq/$fastq2
    FASTQ1=$dirSample/fastq/$fastq1; FASTQ2=$dirSample/fastq/$fastq2
    """
}

process QUALCONTROL {
    /**
    Filtrado con trimmomatic
    */
    
    input:
    tuple val(sample), val(fastq1), val(fastq2), val(dirSample)
    
    output:
    tuple val(sample), path('*_paired-trim_1.fastq.gz'), path('*_paired-trim_2.fastq.gz'), val(dirSample)
    
    script:
    def output = "${sample}_paired-trim_1.fastq.gz ${sample}_unpaired-trim_1.fastq.gz ${sample}_paired-trim_2.fastq.gz ${sample}_unpaired-trim_2.fastq.gz"
    """
    #!/bin/bash
    trimmomatic PE -threads 24 $fastq1 $fastq2 $output LEADING:30 TRAILING:30 SLIDINGWINDOW:10:30
    """
}

process FILTHOST {
    /**
    Filtrado de host reads
    */
    
    input:
    tuple val(sample), path(fastq1), path(fastq2), val(dirSample)
    
    output:
    tuple val(sample), path('*_host_removed_R1.fastq.gz'), path('*_host_removed_R2.fastq.gz'), val(dirSample)
    """
    #!/bin/bash  
    bowtie2 -p 24 -x $params.references.refHuman -1 $fastq1 -2 $fastq2 --un-conc-gz ${sample}_host_removed > ${sample}_host_mapped_and_unmapped.sam
    mv ${sample}_host_removed.1 ${sample}_host_removed_R1.fastq.gz; mv ${sample}_host_removed.2 ${sample}_host_removed_R2.fastq.gz
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


