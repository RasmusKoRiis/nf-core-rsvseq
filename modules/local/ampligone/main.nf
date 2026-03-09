process AMPLIGONE {
    tag "$meta.id"
    label 'process_medium'
    errorStrategy 'ignore'

    container 'quay.io/biocontainers/ampligone:1.3.1--pyhdfd78af_0'

    input:
    tuple val(meta), path(fastq)
    path(primerdir)
    path(primer_schemes_dir)
    val(primer_scheme)

    output:
    tuple val(meta), path("*fastq") , emit: primertrimmedfastq

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    echo ${meta.id}
    if [[ ${meta.id} == *A* ]]; then
        legacy_subdir="a"
        scheme_subdir="RSVA"
    elif [[ ${meta.id} == *B* ]]; then
        legacy_subdir="b"
        scheme_subdir="RSVB"
    else
        echo "meta.id does not contain 'A' or 'B'. Exiting."
        exit 1
    fi

    primer_scheme_value="${primer_scheme:-}"
    if [[ "\$primer_scheme_value" == "null" ]]; then
        primer_scheme_value=""
    fi

    if [[ -n "\$primer_scheme_value" ]]; then
        scheme_dir="${primer_schemes_dir}/\$scheme_subdir/\$primer_scheme_value"
        if [[ ! -d "\$scheme_dir" ]]; then
            echo "Primer scheme directory not found: \$scheme_dir"
            exit 1
        fi

        if [[ -f "\$scheme_dir/reference.fasta" ]]; then
            reference_file="\$scheme_dir/reference.fasta"
        elif [[ -f "\$scheme_dir/\$scheme_subdir.reference.fasta" ]]; then
            reference_file="\$scheme_dir/\$scheme_subdir.reference.fasta"
        else
            echo "Reference FASTA not found in \$scheme_dir"
            exit 1
        fi

        if [[ -f "\$scheme_dir/primer.fasta" ]]; then
            primer_file="\$scheme_dir/primer.fasta"
        elif [[ -f "\$scheme_dir/\$scheme_subdir.primer.fasta" ]]; then
            primer_file="\$scheme_dir/\$scheme_subdir.primer.fasta"
        else
            primer_bed=""
            if [[ -f "\$scheme_dir/primer.bed" ]]; then
                primer_bed="\$scheme_dir/primer.bed"
            elif [[ -f "\$scheme_dir/\$scheme_subdir.primer.bed" ]]; then
                primer_bed="\$scheme_dir/\$scheme_subdir.primer.bed"
            fi

            if [[ -z "\$primer_bed" ]]; then
                echo "No primer FASTA or BED found in \$scheme_dir"
                exit 1
            fi

            primer_file="${meta.id}_generated_primers.fasta"
            python - "\$reference_file" "\$primer_bed" "\$primer_file" <<'PY'
import sys

ref_path, bed_path, out_path = sys.argv[1], sys.argv[2], sys.argv[3]

def revcomp(seq):
    table = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    return seq.translate(table)[::-1]

refs = {}
current = None
chunks = []
with open(ref_path, "r", encoding="utf-8") as fh:
    for line in fh:
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if current is not None:
                refs[current] = "".join(chunks)
            current = line[1:].split()[0]
            chunks = []
        else:
            chunks.append(line)
    if current is not None:
        refs[current] = "".join(chunks)

written = 0
with open(bed_path, "r", encoding="utf-8") as bed, open(out_path, "w", encoding="utf-8") as out:
    for idx, raw in enumerate(bed, start=1):
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.split("\t")
        if len(parts) < 3:
            continue

        chrom = parts[0]
        start = int(parts[1])
        end = int(parts[2])
        name = parts[3] if len(parts) >= 4 and parts[3] else f"primer_{idx}"
        strand = parts[5] if len(parts) >= 6 else "+"

        if chrom not in refs:
            raise SystemExit(f"BED contig '{chrom}' not found in reference FASTA")
        if start < 0 or end <= start or end > len(refs[chrom]):
            raise SystemExit(f"Invalid BED interval: {chrom}:{start}-{end}")
        if (end - start) > 120:
            raise SystemExit(
                f"BED interval is too long for a primer ({chrom}:{start}-{end}). "
                "Use a true primer BED or provide primer.fasta."
            )

        seq = refs[chrom][start:end]
        if strand == "-":
            seq = revcomp(seq)

        out.write(f">{name}\\n{seq}\\n")
        written += 1

if written == 0:
    raise SystemExit("No primers were written from BED")
PY

            if [[ ! -s "\$primer_file" ]]; then
                echo "Failed to generate primer FASTA from \$primer_bed"
                exit 1
            fi
        fi
    else
        reference_file="${primerdir}/\$legacy_subdir/reference.fasta"
        primer_file="${primerdir}/\$legacy_subdir/primer.fasta"
    fi

    if [[ ! -f "\$reference_file" ]]; then
        echo "Reference FASTA missing: \$reference_file"
        exit 1
    fi

    if [[ ! -f "\$primer_file" ]]; then
        echo "Primer FASTA missing: \$primer_file"
        exit 1
    fi

    ampligone \
        --input "$fastq" \
        --output ${meta.id}.fastq \
        --reference "\$reference_file" \
        --primers "\$primer_file" \
        --export-primers ${meta.id}_removed_coordinates.bed \
        --amplicon-type end-to-end
    """
}
