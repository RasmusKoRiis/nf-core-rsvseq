#!/usr/bin/env python3
import sys


def revcomp(seq):
    table = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    return seq.translate(table)[::-1]


def read_fasta(path):
    refs = {}
    current = None
    chunks = []
    with open(path, "r", encoding="utf-8") as fh:
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
    return refs


def main():
    if len(sys.argv) != 4:
        raise SystemExit("Usage: bed_to_primer_fasta.py <reference.fasta> <primer.bed> <out.fasta>")

    ref_path, bed_path, out_path = sys.argv[1], sys.argv[2], sys.argv[3]
    refs = read_fasta(ref_path)

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

            out.write(f">{name}\n{seq}\n")
            written += 1

    if written == 0:
        raise SystemExit("No primers were written from BED")


if __name__ == "__main__":
    main()
