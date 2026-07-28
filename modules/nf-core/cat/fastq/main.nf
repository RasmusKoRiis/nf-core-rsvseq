process CAT_FASTQ {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
        'nf-core/ubuntu:20.04' }"

    input:
    tuple val(meta), path(reads, stageAs: "input*/*")

    output:
    tuple val(meta), path("*.merged.fastq.gz"), emit: reads
    tuple val(meta), path("*.skipped_fastq.log"), emit: skipped_fastq_log
    path "versions.yml"                       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def readList = reads instanceof List ? reads.collect{ it.toString() } : [reads.toString()]
    def quote = { v -> "'${v.replace("'", "'\\''")}'" }
    if (meta.single_end) {
        if (readList.size >= 1) {
            """
            valid_reads=()
            : > ${prefix}.skipped_fastq.log

            for read in ${readList.collect{ quote(it) }.join(' ')}; do
                if gzip -t "\$read" >> ${prefix}.skipped_fastq.log 2>&1; then
                    valid_reads+=("\$read")
                else
                    echo "Skipping corrupt FASTQ for ${prefix}: \$read -> \$(readlink -f "\$read")" | tee -a ${prefix}.skipped_fastq.log >&2
                fi
            done

            if [[ "\${#valid_reads[@]}" -eq 0 ]]; then
                echo "No valid FASTQ files remain for ${prefix} after gzip integrity filtering." >&2
                cat ${prefix}.skipped_fastq.log >&2
                exit 1
            fi

            cat "\${valid_reads[@]}" > ${prefix}.merged.fastq.gz
            gzip -t ${prefix}.merged.fastq.gz

            cat <<-END_VERSIONS > versions.yml
            "${task.process}":
                cat: \$(echo \$(cat --version 2>&1) | sed 's/^.*coreutils) //; s/ .*\$//')
                gzip: \$(echo \$(gzip --version 2>&1) | sed 's/^gzip //; s/ .*\$//')
            END_VERSIONS
            """
        }
    } else {
        if (readList.size >= 2) {
            def read1 = []
            def read2 = []
            readList.eachWithIndex{ v, ix -> ( ix & 1 ? read2 : read1 ) << v }
            """
            read1=(${read1.collect{ quote(it) }.join(' ')})
            read2=(${read2.collect{ quote(it) }.join(' ')})
            valid_read1=()
            valid_read2=()
            : > ${prefix}.skipped_fastq.log

            for i in "\${!read1[@]}"; do
                r1="\${read1[\$i]}"
                r2="\${read2[\$i]}"
                r1_log="${prefix}.gzip_r1_\${i}.log"
                r2_log="${prefix}.gzip_r2_\${i}.log"
                r1_ok=1
                r2_ok=1
                gzip -t "\$r1" > "\$r1_log" 2>&1 && r1_ok=0
                gzip -t "\$r2" > "\$r2_log" 2>&1 && r2_ok=0

                if [[ "\$r1_ok" -eq 0 && "\$r2_ok" -eq 0 ]]; then
                    valid_read1+=("\$r1")
                    valid_read2+=("\$r2")
                else
                    echo "Skipping corrupt FASTQ pair for ${prefix}: \$r1 -> \$(readlink -f "\$r1") / \$r2 -> \$(readlink -f "\$r2")" | tee -a ${prefix}.skipped_fastq.log >&2
                    cat "\$r1_log" "\$r2_log" >> ${prefix}.skipped_fastq.log
                fi
                rm -f "\$r1_log" "\$r2_log"
            done

            if [[ "\${#valid_read1[@]}" -eq 0 ]]; then
                echo "No valid FASTQ pairs remain for ${prefix} after gzip integrity filtering." >&2
                cat ${prefix}.skipped_fastq.log >&2
                exit 1
            fi

            cat "\${valid_read1[@]}" > ${prefix}_1.merged.fastq.gz
            cat "\${valid_read2[@]}" > ${prefix}_2.merged.fastq.gz
            gzip -t ${prefix}_1.merged.fastq.gz
            gzip -t ${prefix}_2.merged.fastq.gz

            cat <<-END_VERSIONS > versions.yml
            "${task.process}":
                cat: \$(echo \$(cat --version 2>&1) | sed 's/^.*coreutils) //; s/ .*\$//')
                gzip: \$(echo \$(gzip --version 2>&1) | sed 's/^gzip //; s/ .*\$//')
            END_VERSIONS
            """
        }
    }

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def readList = reads instanceof List ? reads.collect{ it.toString() } : [reads.toString()]
    if (meta.single_end) {
        if (readList.size >= 1) {
            """
            touch ${prefix}.merged.fastq.gz
            touch ${prefix}.skipped_fastq.log

            cat <<-END_VERSIONS > versions.yml
            "${task.process}":
                cat: \$(echo \$(cat --version 2>&1) | sed 's/^.*coreutils) //; s/ .*\$//')
            END_VERSIONS
            """
        }
    } else {
        if (readList.size >= 2) {
            """
            touch ${prefix}_1.merged.fastq.gz
            touch ${prefix}_2.merged.fastq.gz
            touch ${prefix}.skipped_fastq.log

            cat <<-END_VERSIONS > versions.yml
            "${task.process}":
                cat: \$(echo \$(cat --version 2>&1) | sed 's/^.*coreutils) //; s/ .*\$//')
            END_VERSIONS
            """
        }
    }

}
