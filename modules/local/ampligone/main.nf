process AMPLIGONE {
    tag "$meta.id"
    label 'process_medium'
    //errorStrategy 'ignore'

    container 'quay.io/biocontainers/ampligone:1.3.1--pyhdfd78af_0'

    input:
    tuple val(meta), path(fastq)
    path(primerdir)

    output:
    tuple val(meta), path("*fastq") , emit: primertrimmedfastq

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    # Check if meta.id contains "A" or "B"
    echo ${meta.id}
    if [[ ${meta.id} == *A* ]]; then
        var="a"
    elif [[ ${meta.id} == *B* ]]; then
        var="b"
    else
        echo "meta.id does not contain 'A' or 'B'. Exiting."
        exit 1
    fi

    ampligone \
        --input $fastq\
        --output ${meta.id}.fastq \
        --reference "$primerdir/\$var/reference.fasta" \
        --primers "$primerdir/\$var/primer.fasta" \
        --export-primers ${meta.id}_removed_coordinates.bed \
        --amplicon-type end-to-end
    """
}