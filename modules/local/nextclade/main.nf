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
    set -euo pipefail

    sample_id="${meta.id}"
    sample_subtype="${meta.subtype ?: ''}"

    dataset_candidates=()
    label_candidates=()
    if [[ "\$sample_subtype" == "RSVA" ]]; then
        dataset_candidates=("nextstrain/rsv/a/EPI_ISL_412866" "nextstrain/rsv/b/EPI_ISL_1653999")
        label_candidates=("RSVA" "RSVB")
    elif [[ "\$sample_subtype" == "RSVB" ]]; then
        dataset_candidates=("nextstrain/rsv/b/EPI_ISL_1653999" "nextstrain/rsv/a/EPI_ISL_412866")
        label_candidates=("RSVB" "RSVA")
    elif [[ "$sample_id" == *A* ]]; then
        dataset_candidates=("nextstrain/rsv/a/EPI_ISL_412866")
        label_candidates=("RSVA")
    elif [[ "$sample_id" == *B* ]]; then
        dataset_candidates=("nextstrain/rsv/b/EPI_ISL_1653999")
        label_candidates=("RSVB")
    else
        dataset_candidates=("nextstrain/rsv/a/EPI_ISL_412866" "nextstrain/rsv/b/EPI_ISL_1653999")
        label_candidates=("RSVA" "RSVB")
    fi

    success=0
    last_err=""
    for idx in "\${!dataset_candidates[@]}"; do
        dataset_name="\${dataset_candidates[\$idx]}"
        subtype_label="\${label_candidates[\$idx]}"
        ds_dir="\${sample_id}_nextclade_dataset_\${subtype_label}"
        out_dir="\${sample_id}_nextclade_output_\${subtype_label}"

        rm -rf "\$ds_dir" "\$out_dir"

        if ! nextclade dataset get --name "\$dataset_name" --output-dir "\$ds_dir"; then
            last_err="Failed dataset download for \$dataset_name"
            continue
        fi

        if ! nextclade run \
            --input-dataset "\$ds_dir" \
            --output-all="\$out_dir" \
            $fasta; then
            last_err="Nextclade run failed for \$dataset_name"
            continue
        fi

        if compgen -G "\$out_dir/*" > /dev/null; then
            for file in "\$out_dir"/*; do
                base=\$(basename "\$file")
                mv "\$file" "./\${sample_id}_\${base}"
            done
            success=1
            break
        else
            last_err="No output files produced for \$dataset_name"
        fi
    done

    if [[ "\$success" -ne 1 ]]; then
        echo "\${last_err:-All Nextclade attempts failed}"
        exit 1
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        : \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//' ))
    END_VERSIONS

    """
}
