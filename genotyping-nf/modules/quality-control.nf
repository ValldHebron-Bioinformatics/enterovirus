#!/usr/bin/env nextflow
    
nextflow.enable.dsl = 2

process CREATEDIR {
    /**
    Create directories and set variables
    */

    input:
    tuple val(sampleId), val(file1), val(file2)

    output:
    tuple val(sampleId), val(file1), val(file2), env('outputDir')

    script:
    """
    #!/bin/bash
    outputDir=$params.workdir/$params.user/$sampleId/
    fastqDir=\$outputDir/fastq
    mkdir -p \$outputDir

    if [[ $file1 == *.fastq.gz ]]; then
        extension="fastq"
        if [[ $params.fileType == "fasta" ]]; then
            echo "File $file1 is in FASTQ format, but the input file type is set to FASTA."
            exit 1
        fi
    elif [[ $file1 == *.fastq ]]; then
        extension="fastq"
        if [[ $params.fileType == "fasta" ]]; then
            echo "File $file1 is in FASTQ format, but the input file type is set to FASTA."
            exit 1
        fi
    elif [[ $file1 == *.fasta ]]; then
        extension="fasta"
        if [[ $params.fileType == "fastq" ]]; then
            echo "File $file1 is in FASTA format, but the input file type is set to FASTQ."
            exit 1
        fi
    else
        echo "Unsupported file extension for $file1"
        exit 1
    fi

    # Check file1 path is absolute
    if [[ $file1 != /* ]]; then
        echo "File $file1 is not an absolute path."
        exit 1
    fi

    if [ ! -f $file1 ]; then
        echo "File $file1 does not exist."
        exit 1
    fi

    # Create fastq directory if extension is fastq
    if [[ \$extension == "fastq" ]]; then
        mkdir -p \$fastqDir
        ln -s $file1 \$fastqDir/
    fi

    # Check if second file is not '-'
    if [ "$file2" != "-" ]; then
        # Check file2 path is absolute
        if [[ $file2 != /* ]]; then
            echo "File $file2 is not an absolute path."
            exit 1
        fi

        if [ ! -f $file2 ]; then
            echo "File $file2 does not exist."
            exit 1
        fi
        ln -s $file2 \$fastqDir/

    else
        if [[ \$extension == "fastq" ]]; then
            echo "File $file2 is not provided, but the input file type is set to FASTQ."
            exit 1
        fi
    fi
    """
}

process QUALCONTROL {
    /**
    Filtrado con trimmomatic
    */
    
    input:
    tuple val(sampleId), val(fastq1), val(fastq2), val(outputDir)

    output:
    tuple val(sampleId), path('*_paired-trim_1.fastq.gz'), path('*_paired-trim_2.fastq.gz'), val(outputDir)
    
    script:
    def output = "${sampleId}_paired-trim_1.fastq.gz ${sampleId}_unpaired-trim_1.fastq.gz ${sampleId}_paired-trim_2.fastq.gz ${sampleId}_unpaired-trim_2.fastq.gz"
    """
    #!/bin/bash
    trimmomatic PE -threads $params.threads $fastq1 $fastq2 $output LEADING:30 TRAILING:30 SLIDINGWINDOW:10:30
    """
}

process FILTHOST {
    /**
    Filtrado de host reads
    */
    publishDir "$outputDir/fastq", pattern: '*_host_removed_R*.fastq.gz', mode: 'copy'
    
    input:
    tuple val(sampleId), path(fastq1), path(fastq2), val(outputDir)

    output:
    tuple val(sampleId), path('*_host_removed_R1.fastq.gz'), path('*_host_removed_R2.fastq.gz'), val(outputDir)
    """
    #!/bin/bash  
    bowtie2 -p $params.threads -x $params.references.refHuman -1 $fastq1 -2 $fastq2 --un-conc-gz ${sampleId}_host_removed > ${sampleId}_host_mapped_and_unmapped.sam
    mv ${sampleId}_host_removed.1 ${sampleId}_host_removed_R1.fastq.gz; mv ${sampleId}_host_removed.2 ${sampleId}_host_removed_R2.fastq.gz
    """
}

process TRIMPRIMERSR {
    /**
    Filtrado de primers
    */
    
    input:
    tuple val(sampleId), path(fastq1), path(fastq2), val(outputDir), path(primers)
    
    output:
    tuple val(sampleId), path('*R1_clean_r.fastq.gz'), path('*R2_clean_r.fastq.gz'), val(outputDir), path(primers)
    """
    #!/bin/bash  
    $params.programs.bbduk in1=$fastq1 in2=$fastq2 out1=${sampleId}_R1_clean_r.fastq.gz out2=${sampleId}_R2_clean_r.fastq.gz ref=$primers k=15 ktrim=r restrictright=30
    """
}

process TRIMPRIMERSL {
    /**
    Filtrado de primers
    */
    
    input:
    tuple val(sampleId), path(fastq1), path(fastq2), val(outputDir), path(primers)
    
    output:
    tuple val(sampleId), path('*R1_clean.fastq.gz'), path('*R2_clean.fastq.gz'), val(outputDir)
    """
    #!/bin/bash  
    $params.programs.bbduk in1=${sampleId}_R1_clean_r.fastq.gz in2=${sampleId}_R2_clean_r.fastq.gz out1=${sampleId}_R1_clean.fastq.gz out2=${sampleId}_R2_clean.fastq.gz ref=$primers k=15 ktrim=l restrictleft=30
    """
}
