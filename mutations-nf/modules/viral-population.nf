#!/usr/bin/env nextflow

process VIRAL_POPULATION {
    errorStrategy 'terminate'

    input:
    tuple path(out_path), val(prot), path(VP1cons), path(EVref)
    tuple file(ref_seq), file(cons)

    output:
    stdout

    script:
    """
    DIR_SAMPLE=${out_path}
    FASTQ=${out_path}/fastq
    SAMPLE=${out_path.baseName}
    if [[ -d \$FASTQ ]]; then
        if [[ -s \$FASTQ/"\$SAMPLE"_R1_clean.fastq.gz ]]; then
            R1=\$FASTQ/"\$SAMPLE"_R1_clean.fastq.gz
            R2=\$FASTQ/"\$SAMPLE"_R2_clean.fastq.gz
        elif [[ -s \$FASTQ/"\$SAMPLE"_host_removed_R1.fastq.gz ]]; then
            R1=\$FASTQ/"\$SAMPLE"_host_removed_R1.fastq.gz
            R2=\$FASTQ/"\$SAMPLE"_host_removed_R2.fastq.gz
        else
            echo No fastqs found
            exit 1
        fi
    else
        echo No fastqs found
        exit 1
    fi

    nextflow run $params.programs.mMf --out_path $params.samples_dir/\$SAMPLE --r1 \$R1 --r2 \$R2 --syn_muts "no" --ref_seq ${ref_seq} --annotate ${params.project_data}/metadata/annotate.tsv
    MUT_FILE=$params.samples_dir/\$SAMPLE/mutations/${out_path.baseName}_mutations.csv
    ANNOT_FILE=$params.samples_dir/\$SAMPLE/mutations/Annotated_mutations.csv
    if (grep -q "Annotated_mutation" \$MUT_FILE); then
        cp \$MUT_FILE \$ANNOT_FILE
        cut -d\$';' -f1,2,3,4,5,6,7,8 \$ANNOT_FILE > \$MUT_FILE
        cp \$ANNOT_FILE $params.samples_dir/\$SAMPLE/results/Annotated_mutations.csv
    fi
    cp \$MUT_FILE $params.samples_dir/\$SAMPLE/results/mutations_${prot}.csv
    """
}
