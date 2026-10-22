#!/usr/bin/env python3
"""
Codon usage analysis for C. bombi (CbWHA1) nuclear-encoded CDS sequences.

For each gene, the longest annotated isoform's CDS is used (avoids
double-counting near-identical isoforms from the same locus). Codon
counts, RSCU (Relative Synonymous Codon Usage), and per-codon-position GC
content (GC1/GC2/GC3) are computed using the standard genetic code
(kinetoplastid nuclear genomes use the standard code; only the
mitochondrial/maxicircle genome uses an alternative code, which is
excluded here since only nuclear chromosomes are analyzed).
"""
import os
import re
from collections import defaultdict

import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio.Data import CodonTable
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# Input/output directories, relative to the working directory by default.
# Override with CODON_INPUT_DIR / CODON_OUTPUT_DIR environment variables.
UP = os.environ.get("CODON_INPUT_DIR", ".")
OUT = os.environ.get("CODON_OUTPUT_DIR", ".")
os.makedirs(OUT, exist_ok=True)
FASTA = f"{UP}/CbWHA1_assembly.final_final.fasta"
GFF = f"{UP}/CbWHA1_full_annotation.final_numbered.TATE_reclassified.GAPDHcorrected.gff3"

STANDARD_TABLE = CodonTable.unambiguous_dna_by_id[1]
STOP_CODONS = set(STANDARD_TABLE.stop_codons)
CODON_TO_AA = dict(STANDARD_TABLE.forward_table)
for stop in STOP_CODONS:
    CODON_TO_AA[stop] = "*"

AA_TO_CODONS = defaultdict(list)
for codon, aa in CODON_TO_AA.items():
    AA_TO_CODONS[aa].append(codon)


def load_genome(fasta_path):
    return {rec.id: rec.seq for rec in SeqIO.parse(fasta_path, "fasta")}


def gene_id_from_mrna(mrna_attr):
    return mrna_attr.split("^")[0]


def load_cds_blocks(gff_path):
    """mRNA_id -> [(chrom, start, end, strand), ...], mRNA_id -> gene_id"""
    blocks = defaultdict(list)
    mrna_gene = {}
    with open(gff_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            if f[2] == "mRNA":
                m_id = re.search(r"ID=([^;]+)", f[8])
                m_parent = re.search(r"Parent=([^;]+)", f[8])
                if m_id and m_parent:
                    mrna_gene[m_id.group(1)] = gene_id_from_mrna(m_parent.group(1))
            elif f[2] == "CDS":
                m_parent = re.search(r"Parent=([^;]+)", f[8])
                if m_parent:
                    blocks[m_parent.group(1)].append((f[0], int(f[3]), int(f[4]), f[6]))
    return blocks, mrna_gene


def build_cds_seq(genome, block_list):
    strand = block_list[0][3]
    chrom = block_list[0][0]
    if chrom not in genome:
        return None
    blocks_sorted = sorted(block_list, key=lambda b: b[1], reverse=(strand == "-"))
    parts = []
    for c, start, end, s in blocks_sorted:
        frag = genome[c][start - 1:end]
        if s == "-":
            frag = frag.reverse_complement()
        parts.append(str(frag))
    return "".join(parts)


def main():
    print("Loading genome and CDS blocks...")
    genome = load_genome(FASTA)
    cds_blocks, mrna_gene = load_cds_blocks(GFF)

    # pick the longest-CDS isoform per gene
    gene_best = {}  # gene_id -> (length, mrna_id, seq)
    for mrna_id, blocks in cds_blocks.items():
        gene_id = mrna_gene.get(mrna_id)
        if gene_id is None:
            continue
        seq = build_cds_seq(genome, blocks)
        if seq is None:
            continue
        total_len = sum(e - s + 1 for _, s, e, _ in blocks)
        if gene_id not in gene_best or total_len > gene_best[gene_id][0]:
            gene_best[gene_id] = (total_len, mrna_id, seq)

    print(f"{len(gene_best)} genes with a usable CDS (longest isoform each)")

    codon_counts = defaultdict(int)
    gc1_sum = gc2_sum = gc3_sum = n_codons_total = 0
    n_used_genes = 0
    n_skipped = 0
    per_gene_gc = []

    for gene_id, (length, mrna_id, seq) in gene_best.items():
        seq = seq.upper()
        # trim to a multiple of 3, drop terminal stop codon if present
        usable_len = (len(seq) // 3) * 3
        seq = seq[:usable_len]
        if len(seq) < 30 or set(seq) - set("ACGT"):
            n_skipped += 1
            continue
        codons = [seq[i:i + 3] for i in range(0, len(seq), 3)]
        if codons and codons[-1] in STOP_CODONS:
            codons = codons[:-1]
        if not codons:
            n_skipped += 1
            continue
        n_used_genes += 1
        g1 = g2 = g3 = n_c = 0
        for c in codons:
            if c in CODON_TO_AA:
                codon_counts[c] += 1
                gc1_sum += c[0] in "GC"
                gc2_sum += c[1] in "GC"
                gc3_sum += c[2] in "GC"
                n_codons_total += 1
                g1 += c[0] in "GC"
                g2 += c[1] in "GC"
                g3 += c[2] in "GC"
                n_c += 1
        if n_c > 0:
            per_gene_gc.append({
                "gene_id": gene_id, "n_codons": n_c,
                "GC1": 100 * g1 / n_c, "GC2": 100 * g2 / n_c, "GC3": 100 * g3 / n_c,
            })

    print(f"Used {n_used_genes} genes ({n_skipped} skipped: too short or ambiguous bases), "
          f"{n_codons_total} codons total")
    print(f"GC1={100*gc1_sum/n_codons_total:.1f}%, GC2={100*gc2_sum/n_codons_total:.1f}%, "
          f"GC3={100*gc3_sum/n_codons_total:.1f}%")

    gc_df = pd.DataFrame(per_gene_gc)
    gc_df["GC12"] = (gc_df["GC1"] + gc_df["GC2"]) / 2
    gc_df.to_csv(f"{OUT}/codon_usage_GC123_per_gene.tsv", sep="\t", index=False)

    # --- Figure: GC1/GC2/GC3 skew ---
    plt.rcParams.update({"font.size": 13})
    fig2, axg = plt.subplots(1, 2, figsize=(11, 5))

    data_gc = [gc_df["GC1"], gc_df["GC2"], gc_df["GC3"]]
    colors_gc = ["#377eb8", "#4daf4a", "#e6550d"]
    vp = axg[0].violinplot(data_gc, showmedians=True, widths=0.7)
    for body, c in zip(vp["bodies"], colors_gc):
        body.set_facecolor(c)
        body.set_alpha(0.8)
        body.set_edgecolor("black")
    vp["cmedians"].set_color("black")
    axg[0].axhline(50, color="grey", linestyle="--", linewidth=0.8, label="50% (no skew)")
    axg[0].set_xticks([1, 2, 3])
    axg[0].set_xticklabels(["GC1", "GC2", "GC3"])
    axg[0].set_ylabel("%GC at codon position")
    axg[0].legend(fontsize=9, loc="lower right")
    means = [d.mean() for d in data_gc]
    for xi, m in zip([1, 2, 3], means):
        axg[0].text(xi, 102, f"mean={m:.1f}%", ha="center", fontsize=9)
    axg[0].set_ylim(0, 108)

    # neutrality plot: GC12 (mean of GC1,GC2) vs GC3 per gene
    axg[1].scatter(gc_df["GC3"], gc_df["GC12"], s=4, alpha=0.25, color="#377eb8")
    from scipy.stats import pearsonr
    r, p = pearsonr(gc_df["GC3"], gc_df["GC12"])
    lims = [0, 100]
    axg[1].plot(lims, lims, color="grey", linestyle="--", linewidth=1, label="y=x (mutation-driven expectation)")
    z = np.polyfit(gc_df["GC3"], gc_df["GC12"], 1)
    xs = np.linspace(gc_df["GC3"].min(), gc_df["GC3"].max(), 100)
    axg[1].plot(xs, np.polyval(z, xs), color="black", linewidth=1.5,
                label=f"observed fit (slope={z[0]:.2f}, r={r:.2f})")
    axg[1].set_xlim(0, 100)
    axg[1].set_ylim(0, 100)
    axg[1].set_xlabel("GC3 (%)")
    axg[1].set_ylabel("GC12 (mean of GC1, GC2) (%)")
    axg[1].legend(fontsize=8, loc="upper left")

    plt.tight_layout()
    plt.savefig(f"{OUT}/codon_usage_GC123_skew.png", dpi=200, bbox_inches="tight")
    print("Saved codon_usage_GC123_skew.png")

    # --- RSCU ---
    rows = []
    for aa, codons in AA_TO_CODONS.items():
        if aa == "*":
            continue
        counts = {c: codon_counts.get(c, 0) for c in codons}
        total = sum(counts.values())
        n_syn = len(codons)
        for c in codons:
            rscu = (counts[c] / total * n_syn) if total > 0 else 0
            rows.append({
                "amino_acid": aa, "codon": c, "count": counts[c],
                "rscu": round(rscu, 3), "n_synonymous_codons": n_syn,
                "pct_usage_within_aa": round(100 * counts[c] / total, 1) if total > 0 else 0,
            })

    df = pd.DataFrame(rows).sort_values(["amino_acid", "rscu"], ascending=[True, False])
    df.to_csv(f"{OUT}/codon_usage_RSCU.tsv", sep="\t", index=False)

    print("\n=== Preferred (highest-RSCU) codon per amino acid ===")
    preferred = df.loc[df.groupby("amino_acid")["rscu"].idxmax()]
    print(preferred[["amino_acid", "codon", "rscu", "pct_usage_within_aa"]].to_string(index=False))

    # --- Figure: RSCU bar chart, grouped by amino acid ---
    plt.rcParams.update({"font.size": 10})
    df_plot = df[df.n_synonymous_codons > 1].copy()  # skip Met/Trp (single codon, RSCU always 1)
    df_plot["label"] = df_plot["codon"] + " (" + df_plot["amino_acid"] + ")"

    fig, ax = plt.subplots(figsize=(16, 5))
    colors = plt.cm.tab20(np.linspace(0, 1, df_plot["amino_acid"].nunique()))
    aa_color = dict(zip(sorted(df_plot["amino_acid"].unique()), colors))
    bar_colors = [aa_color[aa] for aa in df_plot["amino_acid"]]

    ax.bar(range(len(df_plot)), df_plot["rscu"], color=bar_colors, edgecolor="black", linewidth=0.3)
    ax.axhline(1.0, color="grey", linestyle="--", linewidth=0.8)
    ax.set_xticks(range(len(df_plot)))
    ax.set_xticklabels(df_plot["label"], rotation=90, fontsize=8)
    ax.set_ylabel("RSCU")
    ax.set_title(f"Relative Synonymous Codon Usage (RSCU) in {n_used_genes} C. bombi CDS "
                 f"(longest isoform per gene)\nGC3 = {100*gc3_sum/n_codons_total:.1f}%")
    plt.tight_layout()
    plt.savefig(f"{OUT}/codon_usage_RSCU.png", dpi=200, bbox_inches="tight")
    print("\nSaved codon_usage_RSCU.png")


if __name__ == "__main__":
    main()
