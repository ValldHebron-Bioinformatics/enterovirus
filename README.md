# Enterovirus repository

## Prerrequisitos

Ensure that the following programs/packages are installed on your system before running the pipeline:

- nextflow v24.04.3
- python3 y modulos:
    - sys
    - re
    - os
    - SeqIO, Seq, CodonTable 
    - pandas
- Trimmomatic v0.39
- bowtie2 v2.5.1
- bbduk.sh (https://github.com/BioInfoTools/BBMap/tree/master)
- blastn v2.14.0+
- diamond v2.1.10.164
- seqkit v2.5.1

## Uso

Para utilizar los datasets de prueba, se debe crear un directorio donde se guardarán los outputs. Este directorio de destino debe escribirse (full path) en el archivo nextflow.config, en la linea de params.workdir

Para ejecutar el test de uso de multifasta, ejecutar la línea de código:

```
nextflow run enterovirus-genotyping.nf --input fasta --user /fullPath/test/demo-user --fastaFile /fullPath/test/demo-user/multifasta.fasta
```

Si en lugar de adjuntar un archivo se copian las secuencias en la web, se deberán guardar en un archivo, que será el input para el parámetro `--fastaFile`
