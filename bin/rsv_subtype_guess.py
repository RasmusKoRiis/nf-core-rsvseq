#!/usr/bin/env python3
import argparse
import gzip
from pathlib import Path


def parse_args():
    p = argparse.ArgumentParser(description="Guess RSV subtype (A/B) from FASTQ reads.")
    p.add_argument("--fastq", required=True, help="Input FASTQ(.gz)")
    p.add_argument("--primerdir", required=True, help="Legacy primer directory containing a/ and b/")
    p.add_argument("--primer-schemes-dir", required=True, help="Primer schemes root containing RSVA/ and RSVB/")
    p.add_argument("--primer-scheme", default="", help="Primer scheme name")
    p.add_argument("--max-reads", type=int, default=4000, help="Max reads to inspect")
    p.add_argument("--kmer", type=int, default=17, help="K-mer size")
    return p.parse_args()


def read_fasta(path: Path) -> str:
    seqs = []
    with path.open("r", encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith(">"):
                continue
            seqs.append(line.upper())
    return "".join(seqs)


def kmer_set(seq: str, k: int):
    out = set()
    if len(seq) < k:
        return out
    for i in range(0, len(seq) - k + 1):
        kmer = seq[i : i + k]
        if "N" in kmer:
            continue
        out.add(kmer)
    return out


def revcomp(seq: str) -> str:
    table = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    return seq.translate(table)[::-1]


def detect_ref_paths(primerdir: Path, primer_schemes_dir: Path, primer_scheme: str):
    refs = {}
    if primer_scheme:
        for subtype in ("RSVA", "RSVB"):
            base = primer_schemes_dir / subtype / primer_scheme
            candidates = [base / "reference.fasta", base / f"{subtype}.reference.fasta"]
            ref = next((c for c in candidates if c.exists()), None)
            if ref:
                refs[subtype] = ref

    if "RSVA" not in refs:
        c = primerdir / "a" / "reference.fasta"
        if c.exists():
            refs["RSVA"] = c
    if "RSVB" not in refs:
        c = primerdir / "b" / "reference.fasta"
        if c.exists():
            refs["RSVB"] = c

    return refs


def score_read(seq: str, kmers_a, kmers_b, k: int):
    if len(seq) < k:
        return 0, 0
    s_a = 0
    s_b = 0
    rc = revcomp(seq)
    for source in (seq, rc):
        for i in range(0, len(source) - k + 1):
            kmer = source[i : i + k]
            if "N" in kmer:
                continue
            if kmer in kmers_a:
                s_a += 1
            if kmer in kmers_b:
                s_b += 1
    return s_a, s_b


def fastq_iter(path: Path):
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt", encoding="utf-8", errors="ignore") as fh:
        while True:
            h = fh.readline()
            if not h:
                return
            seq = fh.readline().strip()
            fh.readline()
            fh.readline()
            if not seq:
                continue
            yield seq.upper()


def main():
    args = parse_args()
    primer_scheme = args.primer_scheme if args.primer_scheme not in ("", "null", "None") else ""
    refs = detect_ref_paths(Path(args.primerdir), Path(args.primer_schemes_dir), primer_scheme)
    if "RSVA" not in refs or "RSVB" not in refs:
        print("UNKNOWN")
        return

    seq_a = read_fasta(refs["RSVA"])
    seq_b = read_fasta(refs["RSVB"])
    if not seq_a or not seq_b:
        print("UNKNOWN")
        return

    kmers_a = kmer_set(seq_a, args.kmer)
    kmers_b = kmer_set(seq_b, args.kmer)
    if not kmers_a or not kmers_b:
        print("UNKNOWN")
        return

    total_a = 0
    total_b = 0
    seen = 0
    for read in fastq_iter(Path(args.fastq)):
        sa, sb = score_read(read, kmers_a, kmers_b, args.kmer)
        total_a += sa
        total_b += sb
        seen += 1
        if seen >= args.max_reads:
            break

    if total_a == 0 and total_b == 0:
        print("UNKNOWN")
        return

    bigger = max(total_a, total_b)
    smaller = min(total_a, total_b)
    ratio = bigger / max(1, smaller)

    # Conservative assignment to avoid confident misclassification.
    if ratio < 1.15 or (bigger - smaller) < 25:
        print("UNKNOWN")
    else:
        print("RSVA" if total_a > total_b else "RSVB")


if __name__ == "__main__":
    main()
