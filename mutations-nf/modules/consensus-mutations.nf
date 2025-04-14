#!/usr/bin/env nextflow

process INPUT_PREPARATION1 {
    errorStrategy 'terminate'
    input:
    tuple path(out_path), val(prot), path(VP1cons), path(EVref)

    output:
    tuple file("ref_*.fasta"), file("cons_VP1.fasta")

    script:
    """
    DIR_SAMPLE=${out_path}
    if [[ -s \${DIR_SAMPLE}/results/species-assignment.csv ]]; then
        gt=\$(python3 -c "import pandas as pd;df = pd.read_csv('\${DIR_SAMPLE}/results/species-assignment.csv', sep=',');print(df.genotype.to_list())")
        if [[ \$gt == "[nan]" ]]; then
            gt=\$(grep ">" ${VP1cons} | tr -d '>' | cut -d\$'_' -f2 | tr -d '-')
            # gt=\$(echo ${VP1cons} | cut -d\$'_' -f1 | tr -d '-')
        fi
    else
        gt=\$(grep ">" ${VP1cons} | tr -d '>' | cut -d\$'_' -f2 | tr -d '-')
        # gt=\$(echo ${VP1cons} | cut -d\$'_' -f1 | tr -d '-')
    fi

    grep "\$gt" ${EVref} | tr -d '>' > name.txt

    # Check if the reference file is empty
    #if [ ! -s name.txt ]; then
    #    echo "Reference file is empty. Exiting."
    #    exit 1
    #fi
    
    seqtk subseq ${EVref} name.txt > "ref_\${gt}.fasta"; rm name.txt

    cp ${VP1cons} cons_VP1.fasta
    cat cons_VP1.fasta
    cat ref_\${gt}.fasta
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
    tuple path(out_path), val(prot), path(VP1cons), path(EVref)
    tuple file(aln_fasta), file(aln_mafft), file(reference), file(consensus)

    output:
    tuple file("${out_path}/results/mutations_${prot}.csv"), file("${out_path}/results/Annotated_mutations.csv"), optional: true

    script:
    """
    RESULTS_DIR=${out_path}/results
    cat ${reference}
    cat ${consensus}
    python3 $params.programs.muts --ref_seq ${reference} --consensus_seq ${consensus} --sample_name ${out_path.baseName} --out_csv \$RESULTS_DIR/mutations_${prot}.csv --prot_name ${prot}
    python3 $params.programs.annotator --out_dir \$RESULTS_DIR --annotate ${params.project_data}/metadata/annotate.csv --muts_file \$RESULTS_DIR/mutations_${prot}.csv
    """
}
