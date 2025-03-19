#!/usr/bin/env nextflow
    
nextflow.enable.dsl = 2

process BLASTN {
    /**
    Mapeo contra base de datos
    */
    //publishDir "$dirSample/assembly", pattern: 'out-blastn.txt', mode: 'copy'
    errorStrategy 'ignore'

    input:
    tuple val(dirFASTA), val(seqsFasta)

    output:
    tuple val(dirFASTA), val(seqsFasta), env('blastnout')

    script:
    if (params.input == 'fasta')
        """
        #!/bin/bash
        blastn -query $seqsFasta -db $params.references.EVdb -out $dirFASTA/out-blastn.txt -outfmt "6 qacc sacc score evalue qstart qend sstart send"
        blastn -task dc-megablast -query $seqsFasta -db $params.references.VP1db -out $dirFASTA/out-blastn2.txt -outfmt "6 qacc sacc score evalue qstart qend sstart send"
        if [ -f $dirFASTA/out-blastn.txt ] && [ -f $dirFASTA/out-blastn2.txt ]; then
            cat $dirFASTA/out-blastn2.txt >> $dirFASTA/out-blastn.txt
            blastnout=$dirFASTA/out-blastn.txt
        elif [ -f $dirFASTA/out-blastn.txt ] && [ ! -f $dirFASTA/out-blastn2.txt ]; then
            blastnout=$dirFASTA/out-blastn.txt
        elif [ ! -f $dirFASTA/out-blastn.txt ] && [ -f $dirFASTA/out-blastn2.txt ]; then
            blastnout=$dirFASTA/out-blastn2.txt
        fi
        if [ \$(cat \$blastnout | wc -l) -eq 0 ]; then
            echo "No enterovirus sequences found." >> ${dirFASTA}/errors.log
            exit 1
        fi
        """
    else if (params.input == 'fastq')
        if (params.protocol == 'complete')
            """
            #!/bin/bash
            blastn -query $seqsFasta -db $params.references.EVdb -out $dirFASTA/out-blastn.txt -outfmt "6 qacc sacc score evalue qstart qend sstart send"
            if [ -f $dirFASTA/out-blastn.txt ]; then
                blastnout=$dirFASTA/out-blastn.txt
            fi
            if [ \$(cat \$blastnout | wc -l) -eq 0 ]; then
                echo "No enterovirus sequences found." >> ${dirFASTA}/../errors.log
                exit 1
            fi
            """
        else if (params.protocol == 'partial')
            """
            blastn -task dc-megablast -query $seqsFasta -db $params.references.VP1db -out $dirFASTA/out-blastn.txt -outfmt "6 qacc sacc score evalue qstart qend sstart send"
            if [ -f $dirFASTA/out-blastn.txt ]; then
                blastnout=$dirFASTA/out-blastn.txt
            fi
            if [ \$(cat \$blastnout | wc -l) -eq 0 ]; then
                echo "No enterovirus sequences found." >> ${dirFASTA}/../errors.log
                exit 1
            fi
            """
}

process GETBLASTNMATCH {
    /**
    Recuperar los match de blastn
    */

    input:
    tuple val(dirFASTA), val(seqsFasta), val(blastnout)

    output:
    val(dirFASTA)

    script:
    """
    #!/bin/bash
    python3 $params.programs.generateFastas --blast $blastnout --scaffolds $seqsFasta --out-dir $dirFASTA --protocol $params.protocol --input $params.input --refs $params.references.speciesType
    if [[ $params.input == "fastq" ]]; then
        cp $dirFASTA/ev-match.fasta $dirFASTA/../results/ev-match.fasta
    else
        cp $dirFASTA/ev-match.fasta $dirFASTA/results/ev-match.fasta
    fi
    """
}

process GETCDS {

    input:
    val(dirFASTA)

    output:
    tuple val(dirFASTA), env('cdsFile')

    script:
    """
    #!/bin/bash
    python3 $params.programs.translateSeq --seq $dirFASTA/ev-match.fasta --fastaNucl $dirFASTA/ev-match_cds-nucl.fasta --fastaProt $dirFASTA/ev-match_cds-aa.fasta
    cdsFile=$dirFASTA/ev-match_cds-aa.fasta
    """
}

process DIAMOND {

    input:
    tuple val(dirFASTA), val(cdsFile)

    output:
    val(dirFASTA)

    script:
    """
    #!/bin/bash
    $params.programs.diamond blastp --query $cdsFile --db $params.references.VP1dbDiamond --out $dirFASTA/out-diamond.txt --outfmt 6 qseqid sseqid bitscore evalue qstart qend sstart send
    """
}

process GENOTYPEVP1 {

    input:
    val(dirFASTA)

    output:
    val(dirFASTA)

    script:
    """
    #!/bin/bash
    python3 $params.programs.getVP1 --dir $dirFASTA --diamond $dirFASTA/out-diamond.txt --refs $params.references.speciesType --pwd $params.workdir
    awk -v outdir=$dirFASTA '/^>/ {out = outdir "/" substr(\$1, 2) ".fasta"; print > out} !/^>/ {print >> out}' $dirFASTA/VP1_nucl.fasta
    if [[ $params.input == "fastq" ]]; then
        cp $dirFASTA/species-assignment.csv $dirFASTA/../results/genotype-assignment.csv
    else
        cp $dirFASTA/species-assignment.csv $dirFASTA/results/genotype-assignment.csv
    fi
    """
}
