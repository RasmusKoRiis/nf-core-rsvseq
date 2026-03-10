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
    primer_scheme_value="!{primer_scheme}"
    primerdir="!{primerdir}"
    primer_schemes_dir="!{primer_schemes_dir}"

    if [[ "$sample_id" == *A* ]]; then
        legacy_subdir="a"
        scheme_subdir="RSVA"
    elif [[ "$sample_id" == *B* ]]; then
        legacy_subdir="b"
        scheme_subdir="RSVB"
    else
        echo "meta.id does not contain 'A' or 'B': $sample_id"
        exit 1
    fi

    if [[ "$primer_scheme_value" == "null" ]]; then
        primer_scheme_value=""
    fi

    if [[ -n "$primer_scheme_value" ]]; then
        scheme_dir="${primer_schemes_dir}/${scheme_subdir}/${primer_scheme_value}"
        if [[ ! -d "$scheme_dir" ]]; then
            echo "Primer scheme directory not found: $scheme_dir"
            exit 1
        fi

        if [[ -f "$scheme_dir/reference.fasta" ]]; then
            reference_file="$scheme_dir/reference.fasta"
        elif [[ -f "$scheme_dir/${scheme_subdir}.reference.fasta" ]]; then
            reference_file="$scheme_dir/${scheme_subdir}.reference.fasta"
        else
            echo "Reference FASTA not found in $scheme_dir"
            exit 1
        fi

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

            if [[ -z "$primer_bed" ]]; then
                echo "No primer FASTA or BED found in $scheme_dir"
                exit 1
            fi

            primer_file="${sample_id}_generated_primers.fasta"
            bed_to_primer_fasta.py "$reference_file" "$primer_bed" "$primer_file"

            if [[ ! -s "$primer_file" ]]; then
                echo "Failed to generate primer FASTA from $primer_bed"
                exit 1
            fi
        fi
    else
        reference_file="${primerdir}/${legacy_subdir}/reference.fasta"
        primer_file="${primerdir}/${legacy_subdir}/primer.fasta"
    fi

    if [[ ! -f "$reference_file" ]]; then
        echo "Reference FASTA missing: $reference_file"
        exit 1
    fi

    if [[ ! -f "$primer_file" ]]; then
        echo "Primer FASTA missing: $primer_file"
        exit 1
    fi

    ampligone \
        --input "!{fastq}" \
        --output "${sample_id}.fastq" \
        --reference "$reference_file" \
        --primers "$primer_file" \
        --export-primers "${sample_id}_removed_coordinates.bed" \
        --amplicon-type end-to-end
    '''
}
