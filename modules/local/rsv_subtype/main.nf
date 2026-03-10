process RSV_SUBTYPE {
    tag { meta.id }
    label 'process_single'
    errorStrategy 'ignore'

    container 'quay.io/biocontainers/ampligone:1.3.1--pyhdfd78af_0'

    input:
    tuple val(meta), path(fastq)
    path(primerdir)
    path(primer_schemes_dir)
    val(primer_scheme)

    output:
    tuple val(meta), path(fastq), env(SUBTYPE), emit: typed_fastq

    when:
    task.ext.when == null || task.ext.when

    shell:
    '''
    set -euo pipefail

    subtype="$(rsv_subtype_guess.py \
        --fastq "!{fastq}" \
        --primerdir "!{primerdir}" \
        --primer-schemes-dir "!{primer_schemes_dir}" \
        --primer-scheme "!{primer_scheme}" || echo UNKNOWN)"

    case "$subtype" in
        RSVA|RSVB) ;;
        *) subtype="UNKNOWN" ;;
    esac

    export SUBTYPE="$subtype"
    '''
}
