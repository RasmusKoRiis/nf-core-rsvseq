process CHOPPER {
    tag "$meta.id"
    label 'process_medium'
    //errorStrategy 'ignore'

    container 'quay.io/biocontainers/chopper:0.9.0--hdcf5f25_0'

    input:
    tuple val(meta), path(fastq)

    output:
    tuple val(meta), path("*fastq") , emit: chopperfastq

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    #fastp -i $fastq -o ${meta.id}_filtered.fastq
    chopper -q 11 -l 700 --maxlength 1500 -i $fastq > ${meta.id}_filtered.fastq
    #cutadapt -u -7 -u 7 -o ${meta.id}_filtered.fastq $fastq
    """
}
