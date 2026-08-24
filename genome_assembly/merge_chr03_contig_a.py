#!/usr/bin/env python3
"""
Merge CbWHA1_contig_a onto CbWHA1_chr03 in the final assembly FASTA.

Background: contig_a:104761-128850 (24,089 bp) aligns to
chr03:376709-400798 at >99.9% identity and is unique elsewhere in the
genome (see the minimap2 commands in chr03_contig_a_merge_commands.txt
for how this was established) -- chr03's own copy of this window is
kept unmodified. contig_a[:104761] (telomere + protein-coding genes +
contig_a's own copy of the miniexon_SL array) is reverse-complemented
and appended after a 100 bp N junction marker, since it has no
counterpart in the current chr03 assembly and chr03 itself lacks a
telomere at this end.

Usage:
    python3 merge_chr03_contig_a.py \\
        CbWHA1_assembly.final.fasta \\
        CbWHA1_assembly.chr03merged.v2.fasta
"""
import sys

GAP_N = 100
KEEP_LEN = 104761  # contig_a[0:104761] kept; contig_a[104761:128850] is the
                    # validated single-copy anchor, redundant with chr03's
                    # own unmodified copy of the same window


def revcomp(s):
    comp = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    return s.translate(comp)[::-1]


def load_fasta(path):
    seqs, order = {}, []
    name, buf = None, []
    with open(path) as f:
        for line in f:
            if line.startswith(">"):
                if name:
                    seqs[name] = "".join(buf)
                name = line[1:].split()[0]
                order.append(name)
                buf = []
            else:
                buf.append(line.strip())
        if name:
            seqs[name] = "".join(buf)
    return seqs, order


def main():
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    in_path, out_path = sys.argv[1:3]

    seqs, order = load_fasta(in_path)
    chr03 = seqs["CbWHA1_chr03"]
    contig_a = seqs["CbWHA1_contig_a"]

    merged_chr03 = chr03 + ("N" * GAP_N) + revcomp(contig_a[:KEEP_LEN])

    seqs["CbWHA1_chr03"] = merged_chr03
    del seqs["CbWHA1_contig_a"]
    order.remove("CbWHA1_contig_a")

    with open(out_path, "w") as f:
        for n in order:
            f.write(f">{n}\n")
            s = seqs[n]
            for i in range(0, len(s), 60):
                f.write(s[i:i + 60] + "\n")

    print(f"chr03 (unchanged): {len(chr03)} bp")
    print(f"gap: {GAP_N} bp")
    print(f"appended contig_a block (revcomp): {KEEP_LEN} bp")
    print(f"merged chr03 total: {len(merged_chr03)} bp")
    print(f"sequences in output: {len(order)}")
    print(f"Wrote: {out_path}")


if __name__ == "__main__":
    main()
