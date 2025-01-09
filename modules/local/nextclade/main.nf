process NEXTCLADE {
    tag "${meta.id}"
    label 'process_single'
    errorStrategy 'ignore'
    
    container 'docker.io/nextstrain/nextclade:latest'
    containerOptions = "-v ${baseDir}/bin:/project-bin" // Mount the bin directory
    //container logic as needed

    input:
    tuple val(meta), path(fasta)
   

    output:
    tuple val(meta), path("*nextclade.csv"), emit: nextclade_csv, optional: true


    when:
    task.ext.when == null || task.ext.when


    script:
    """
    
    nextclade dataset get --name 'nextstrain/sars-cov-2/wuhan-hu-1/orfs' --output-dir "${meta.id}_whuan_nextclade_dataset/"

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
