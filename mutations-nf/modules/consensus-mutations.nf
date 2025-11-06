#!/usr/bin/env nextflow
nextflow.enable.dsl = 2
process INPUT_PREPARATION1 {
    errorStrategy 'terminate'
    input:
    tuple path(out_path), val(prot), path(VP1cons), path(EVref), val(genotype)

    output:
    tuple file("ref_*.fasta"), file("cons_VP1.fasta")

    script:
    """
    #!/bin/bash
    gt=\$(echo $genotype | sed -e 's/-/_/g' -e 's/(/_/g' -e 's/)//g')
    echo ${EVref}
    ref=$EVref
    grep "\$gt" ${EVref} | tr -d '>' > name.txt
    
    seqtk subseq ${EVref} name.txt > "ref_${genotype}.fasta"; rm name.txt

    cp ${VP1cons} cons_VP1.fasta
    cat cons_VP1.fasta
    cat ref_${genotype}.fasta
    """
}


process INPUT_PREPARATION2 {
    errorStrategy 'terminate'
    input:
    tuple file(VP1ref), file(VP1cons)

    output:
    tuple file("alignment.fasta"), file("alignment.mafft"), file("ref.fasta"), file("cons.fasta")

    script:
    """
    cat ${VP1ref} > alignment.fasta
    grep ">" ${VP1ref} | tr -d '>' > refname.txt
    cat ${VP1cons} >> alignment.fasta
    grep ">" ${VP1cons} | tr -d '>' > consensusname.txt
    mafft alignment.fasta > alignment.mafft
    seqtk subseq alignment.mafft refname.txt > ref.fasta; rm refname.txt
    seqtk subseq alignment.mafft consensusname.txt > cons.fasta; rm consensusname.txt
    """
}


process FIND_MUTATIONS {
    errorStrategy 'terminate'
    publishDir "${out_path}/mutations", mode: 'copy', pattern: "${reference}"
    input:
    tuple path(out_path), val(prot), path(VP1cons), path(EVref), val(genotype)
    tuple file(aln_fasta), file(aln_mafft), file(reference), file(consensus)

    output:
    tuple file("${out_path}/results/mutations_${prot}.csv"), file("${out_path}/results/Annotated_mutations.csv"), optional: true

    script:
    """
    RESULTS_DIR=${out_path}/results
    echo ref
    cat ${reference}
    echo cons
    cat ${consensus}

    python3 $params.programs.muts --ref_seq ${reference} --consensus_seq ${consensus} --sample_name ${out_path.baseName}_${genotype} --out_csv \$RESULTS_DIR/mutations_${prot}_${genotype}.csv --prot_name ${prot}
    python3 $params.programs.annotator --out_dir \$RESULTS_DIR --annotate ${params.project_data}/metadata/annotate.csv --sample_name ${out_path.baseName}_${genotype} --muts_file \$RESULTS_DIR/mutations_${prot}_${genotype}.csv 
    """
}
