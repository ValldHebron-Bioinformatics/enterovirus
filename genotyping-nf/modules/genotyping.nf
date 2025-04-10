#!/usr/bin/env nextflow
    
nextflow.enable.dsl = 2

process BLASTN {
    /**
    Mapeo contra base de datos
    */
    //publishDir "$dirSample/assembly", pattern: 'out-blastn.txt', mode: 'copy'
    errorStrategy 'ignore'

    input:
    tuple val(sample_id), val(seqsFasta), val(outputDir), val(extension)

    output:
    tuple val(outputDir), val(seqsFasta), val(extension), env('blastnout')

    script:
    if (extension == 'fasta')
        """
        #!/bin/bash
        blastn -query $seqsFasta -db $params.references.EVdb -out $outputDir/out-blastn.txt -outfmt "6 qacc sacc score evalue qstart qend sstart send"
        blastn -task dc-megablast -query $seqsFasta -db $params.references.VP1db -out $outputDir/out-blastn2.txt -outfmt "6 qacc sacc score evalue qstart qend sstart send"
        if [ -f $outputDir/out-blastn.txt ] && [ -f $outputDir/out-blastn2.txt ]; then
            cat $outputDir/out-blastn2.txt >> $outputDir/out-blastn.txt
            blastnout=$outputDir/out-blastn.txt
        elif [ -f $outputDir/out-blastn.txt ] && [ ! -f $outputDir/out-blastn2.txt ]; then
            blastnout=$outputDir/out-blastn.txt
        elif [ ! -f $outputDir/out-blastn.txt ] && [ -f $outputDir/out-blastn2.txt ]; then
            blastnout=$outputDir/out-blastn2.txt
        fi
        if [ \$(cat \$blastnout | wc -l) -eq 0 ]; then
            echo "No enterovirus sequences found." >> ${outputDir}/errors.log
            exit 1
        fi
        """
    else if (extension == 'fastq')
        if (params.protocol == 'complete')
            """
            #!/bin/bash
            blastn -query $seqsFasta -db $params.references.EVdb -out $outputDir/out-blastn.txt -outfmt "6 qacc sacc score evalue qstart qend sstart send"
            if [ -f $outputDir/out-blastn.txt ]; then
                blastnout=$outputDir/out-blastn.txt
            fi
            if [ \$(cat \$blastnout | wc -l) -eq 0 ]; then
                echo "No enterovirus sequences found." >> ${outputDir}/../errors.log
                exit 1
            fi
            """
        else if (params.protocol == 'partial')
            """
            blastn -task dc-megablast -query $seqsFasta -db $params.references.VP1db -out $outputDir/out-blastn.txt -outfmt "6 qacc sacc score evalue qstart qend sstart send"
            if [ -f $outputDir/out-blastn.txt ]; then
                blastnout=$outputDir/out-blastn.txt
            fi
            if [ \$(cat \$blastnout | wc -l) -eq 0 ]; then
                echo "No enterovirus sequences found." >> ${outputDir}/../errors.log
                exit 1
            fi
            """
}

process GETBLASTNMATCH {
    /**
    Recuperar los match de blastn
    */

    input:
    tuple val(outputDir), val(seqsFasta), val(extension), val(blastnout)

    output:
    val(outputDir)

    script:
    """
    #!/bin/bash
    python3 $params.programs.generateFastas --blast $blastnout --scaffolds $seqsFasta --out-dir $outputDir --protocol $params.protocol --input $extension --refs $params.references.speciesType
    mkdir -p $outputDir/results
    cp $outputDir/ev-match.fasta $outputDir/results/ev-match.fasta
    cp $outputDir/species-assignment.csv $outputDir/results/species-assignment.csv
    """
}

process GETCDS {

    input:
    val(outputDir)

    output:
    tuple val(outputDir), env('cdsFile')

    script:
    """
    #!/bin/bash
    python3 $params.programs.translateSeq --seq $outputDir/ev-match.fasta --fastaNucl $outputDir/ev-match_cds-nucl.fasta --fastaProt $outputDir/ev-match_cds-aa.fasta
    cdsFile=$outputDir/ev-match_cds-aa.fasta
    """
}

process DIAMOND {

    input:
    tuple val(outputDir), val(cdsFile)

    output:
    val(outputDir)

    script:
    """
    #!/bin/bash
    $params.programs.diamond blastp --query $cdsFile --db $params.references.VP1dbDiamond --out $outputDir/out-diamond.txt --outfmt 6 qseqid sseqid bitscore evalue qstart qend sstart send
    """
}

process GENOTYPEVP1 {

    input:
    val(outputDir)

    output:
    val(outputDir)

    script:
    """
    #!/bin/bash
    python3 $params.programs.getVP1 --dir $outputDir --diamond $outputDir/out-diamond.txt --refs $params.references.speciesType --pwd $params.workdir
    awk -v outdir=$outputDir '/^>/ {out = outdir "/" substr(\$1, 2) ".fasta"; print > out} !/^>/ {print >> out}' $outputDir/VP1_nucl.fasta
    cp $outputDir/species-assignment.csv $outputDir/results/genotype-assignment.csv
    """
}
