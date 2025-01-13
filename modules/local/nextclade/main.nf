process NEXTCLADE {
    tag "${meta.id}"
    label 'process_single'
    //errorStrategy 'ignore'
    
    container 'docker.io/nextstrain/nextclade:latest'
    containerOptions = "-v ${baseDir}/bin:/project-bin" // Mount the bin directory
    //container logic as needed

    input:
    tuple val(meta), path(fasta)
   

    output:
    tuple val(meta), path("*nextclade.csv"), emit: nextclade_csv, optional: true
    path("*nextclade.csv"), emit: report_nextclade_csv, optional: true


    when:
    task.ext.when == null || task.ext.when


    script:
    """

    # Check if meta.id contains "A" or "B"
    if [[ ${meta.id} == *A* ]]; then
        var="nextstrain/rsv/a/EPI_ISL_412866"
    elif [[ ${meta.id} == *B* ]]; then
        var="nextstrain/rsv/b/EPI_ISL_1653999"
    else
        echo "meta.id does not contain 'A' or 'B'. Exiting."
        exit 1
    fi

    nextclade dataset get --name "\$var" --output-dir "${meta.id}_whuan_nextclade_dataset/"

    nextclade run \
        --input-dataset ${meta.id}_whuan_nextclade_dataset/ \
        --output-all=${meta.id}_whuan_nextclade_output/ \
        $fasta

    if compgen -G "${meta.id}_whuan_nextclade_output/*" > /dev/null; then
        for file in ${meta.id}_whuan_nextclade_output/*; do
            cat "\$file"
            basename=\$(basename \$file)
            if [[ "\$file" == *.csv ]]; then
                mv "\$file" ./${meta.id}_\$basename
            else
                mv "\$file" ./${meta.id}_\$basename
            fi
        done
    fi

    

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        : \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//' ))
    END_VERSIONS

    """
}
