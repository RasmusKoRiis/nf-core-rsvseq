process AMPLIGONE {
    tag { meta.id }
    label 'process_medium'
    //errorStrategy 'ignore'

    container 'quay.io/biocontainers/ampligone:1.3.1--pyhdfd78af_0'

    input:
    tuple val(meta), path(fastq)
    path(primerdir)
    path(primer_schemes_dir)
    val(primer_scheme)

    output:
    tuple val(meta), path("*.fastq"), emit: primertrimmedfastq

    when:
    task.ext.when == null || task.ext.when

    shell:
    '''
    set -euo pipefail

    sample_id="!{meta.id}"
    sample_subtype="!{meta.subtype ?: ''}"
    primer_scheme_value="!{primer_scheme}"
    primerdir="!{primerdir}"
    primer_schemes_dir="!{primer_schemes_dir}"

    scheme_candidates=()
    legacy_candidates=()
    if [[ "$sample_subtype" == "RSVA" ]]; then
        scheme_candidates=("RSVA" "RSVB")
        legacy_candidates=("a" "b")
    elif [[ "$sample_subtype" == "RSVB" ]]; then
        scheme_candidates=("RSVB" "RSVA")
        legacy_candidates=("b" "a")
    elif [[ "$sample_id" == *A* ]]; then
        scheme_candidates=("RSVA")
        legacy_candidates=("a")
    elif [[ "$sample_id" == *B* ]]; then
        scheme_candidates=("RSVB")
        legacy_candidates=("b")
    else
        # Sample IDs like RSV105 do not encode subtype; try both.
        scheme_candidates=("RSVA" "RSVB")
        legacy_candidates=("a" "b")
    fi

    if [[ "$primer_scheme_value" == "null" ]]; then
        primer_scheme_value=""
    fi

    run_ampligone() {
        local ref="$1"
        local pri="$2"
        local mode="$3"
        local log_file="${sample_id}.ampligone_${mode}.log"
        if ampligone \
            --input "!{fastq}" \
            --output "${sample_id}.fastq" \
            --reference "$ref" \
            --primers "$pri" \
            --export-primers "${sample_id}_removed_coordinates.bed" \
            --amplicon-type end-to-end >"$log_file" 2>&1; then
            return 0
        fi
        return 1
    }

    attempted=0
    last_error=""

    if [[ -n "$primer_scheme_value" ]]; then
        for scheme_subdir in "${scheme_candidates[@]}"; do
            scheme_dir="${primer_schemes_dir}/${scheme_subdir}/${primer_scheme_value}"
            [[ -d "$scheme_dir" ]] || continue

            reference_file=""
            if [[ -f "$scheme_dir/reference.fasta" ]]; then
                reference_file="$scheme_dir/reference.fasta"
            elif [[ -f "$scheme_dir/${scheme_subdir}.reference.fasta" ]]; then
                reference_file="$scheme_dir/${scheme_subdir}.reference.fasta"
            fi
            [[ -n "$reference_file" ]] || continue

            primer_file=""
            if [[ -f "$scheme_dir/primer.fasta" ]]; then
                primer_file="$scheme_dir/primer.fasta"
            elif [[ -f "$scheme_dir/${scheme_subdir}.primer.fasta" ]]; then
                primer_file="$scheme_dir/${scheme_subdir}.primer.fasta"
            else
                primer_bed=""
                if [[ -f "$scheme_dir/primer.bed" ]]; then
                    primer_bed="$scheme_dir/primer.bed"
                elif [[ -f "$scheme_dir/${scheme_subdir}.primer.bed" ]]; then
                    primer_bed="$scheme_dir/${scheme_subdir}.primer.bed"
                fi
                [[ -n "$primer_bed" ]] || continue

                primer_file="${sample_id}_${scheme_subdir}_generated_primers.fasta"
                bed_to_primer_fasta.py "$reference_file" "$primer_bed" "$primer_file"
                [[ -s "$primer_file" ]] || continue
            fi

            attempted=$((attempted + 1))
            if run_ampligone "$reference_file" "$primer_file" "$scheme_subdir"; then
                exit 0
            fi
            last_error="ampligone failed for scheme ${scheme_subdir} (see ${sample_id}.ampligone_${scheme_subdir}.log)"
        done
    else
        for legacy_subdir in "${legacy_candidates[@]}"; do
            reference_file="${primerdir}/${legacy_subdir}/reference.fasta"
            primer_file="${primerdir}/${legacy_subdir}/primer.fasta"
            [[ -f "$reference_file" && -f "$primer_file" ]] || continue

            attempted=$((attempted + 1))
            if run_ampligone "$reference_file" "$primer_file" "$legacy_subdir"; then
                exit 0
            fi
            last_error="ampligone failed for legacy ${legacy_subdir} (see ${sample_id}.ampligone_${legacy_subdir}.log)"
        done
    fi

    if [[ "$attempted" -eq 0 ]]; then
        echo "No valid reference/primer combination found for sample $sample_id"
        exit 1
    fi

    echo "${last_error:-All candidate ampligone runs failed for $sample_id}"
    exit 1
    '''
}
