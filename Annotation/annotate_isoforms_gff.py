#!/usr/bin/env python3
"""
Add gene-level product/ortholog annotation to EVERY isoform of a gene,
without collapsing/dropping isoform records.

The correct model here (trypanosomatids are effectively intron-less) is: one gene = one ORF,
and different assembled "isoforms" mostly reflect different trans-splice
acceptor / polyadenylation site choices -- i.e. different 5'/3' UTR
boundaries around the SAME coding sequence, not different proteins. So
isoform records (with their own real UTR/exon boundaries) are worth
keeping; only the CDS needs to be single-valued per gene, and the gene
product/ortholog annotation belongs on all of them, not just one pick.

Inputs:
    1. The CDS-harmonized GFF3 from harmonize_isoform_CDS.py -- all
        original isoform mRNA records, CDS already forced to match
       within CD-HIT-confirmed same-gene clusters (1,591 genes / 1,684
       isoforms). Left AS-IS (independently-called CDS per isoform) for
       the 352 multi-isoform genes CD-HIT did NOT cluster -- see below.
    2. CfCl/Ld/Lm/Lp/Ls/Tb-ortholog-transferred GFF3
       (gene_product + N_ortholog attributes on CDS lines), matched by
       gene ID (her file already picked one isoform per gene; the
       annotation is a gene-level fact regardless of which isoform
       carried it).

Output: every mRNA/isoform is kept, with every one of a gene's isoforms'
CDS lines carrying the same gene_product/ortholog attributes.

Usage:
    python3 annotate_isoforms_gff.py \\
        Cb_2.5-8_annotation.transdecoder.genome.CDSharmonized.gff3 \\
        transdecoder_collapsed.pep.clstr \\
        Cb_full_annotated_genome_transdecoder_8-20_Cftransfered.gff3 \\
        Cb_2.5-8_annotation.transdecoder.genome.isoforms_annotated.gff3
"""
import re
import sys
from collections import defaultdict

ID_RE = re.compile(r"ID=([^;]+)")
PARENT_RE = re.compile(r"Parent=([^;]+)")

ANNOT_TAGS = [
    "cf_ortholog", "ld_ortholog", "lm_ortholog", "lp_ortholog",
    "ls_ortholog", "tb_ortholog", "gene_product",
]


def parse_clstr(path):
    clusters = []
    cur = []
    with open(path) as f:
        for line in f:
            line = line.rstrip()
            if line.startswith(">Cluster"):
                if cur:
                    clusters.append(cur)
                cur = []
            else:
                m = re.search(r">(\S+)\.\.\.", line)
                if m:
                    cur.append((m.group(1), line.endswith("*")))
        if cur:
            clusters.append(cur)
    return clusters


def load_gff(path):
    lines = []
    mrna_gene = {}
    gene_mrnas = defaultdict(list)
    with open(path) as f:
        for idx, raw in enumerate(f):
            lines.append(raw)
            if raw.startswith("#") or not raw.strip():
                continue
            fields = raw.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue
            ftype, attrs = fields[2], fields[8]
            if ftype == "mRNA":
                idm = ID_RE.search(attrs)
                pm = PARENT_RE.search(attrs)
                if idm and pm:
                    mrna_gene[idm.group(1)] = pm.group(1)
                    gene_mrnas[pm.group(1)].append(idm.group(1))
    return lines, mrna_gene, gene_mrnas


def build_gene_annotations(student_path):
    mrna_gene = {}
    gene_annot = {}
    with open(student_path) as f:
        for line in f:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue
            ftype, attrs = fields[2], fields[8]
            if ftype == "mRNA":
                idm = ID_RE.search(attrs)
                pm = PARENT_RE.search(attrs)
                if idm and pm:
                    mrna_gene[idm.group(1)] = pm.group(1)
            elif ftype == "CDS":
                pm = PARENT_RE.search(attrs)
                if pm and pm.group(1) in mrna_gene:
                    gid = mrna_gene[pm.group(1)]
                    found = {}
                    for tag in ANNOT_TAGS:
                        m = re.search(rf"{tag}=([^;]+)", attrs)
                        if m:
                            found[tag] = m.group(1)
                    if found:
                        gene_annot[gid] = found
    return gene_annot


def main():
    if len(sys.argv) != 5:
        sys.exit(__doc__)
    harmonized_gff, clstr_path, student_gff, out_path = sys.argv[1:5]

    lines, mrna_gene, gene_mrnas = load_gff(harmonized_gff)
    clusters = parse_clstr(clstr_path)
    gene_annot = build_gene_annotations(student_gff)

    # genes covered by a same-gene, multi-member CD-HIT cluster (already
    # CDS-harmonized upstream)
    cdhit_harmonized_genes = set()
    for c in clusters:
        if len(c) < 2:
            continue
        c = [(mid, rep) for mid, rep in c if mid in mrna_gene]
        if len(c) < 2:
            continue
        gene_ids = {mrna_gene[mid] for mid, _ in c}
        if len(gene_ids) == 1:
            cdhit_harmonized_genes.add(next(iter(gene_ids)))

    pending_proteomics = [
        gid for gid, mids in gene_mrnas.items()
        if len(mids) > 1 and gid not in cdhit_harmonized_genes
    ]

    n_isoforms_annotated = 0
    n_genes_annotated = 0
    n_genes_no_annotation = 0

    out_lines = list(lines)  # copy; we'll mutate CDS lines in place by index
    for idx, raw in enumerate(lines):
        if raw.startswith("#") or not raw.strip():
            continue
        fields = raw.rstrip("\n").split("\t")
        if len(fields) < 9 or fields[2] != "CDS":
            continue
        pm = PARENT_RE.search(fields[8])
        if not pm or pm.group(1) not in mrna_gene:
            continue
        gid = mrna_gene[pm.group(1)]
        annot = gene_annot.get(gid)
        if not annot:
            continue
        extra = "".join(f";{tag}={annot[tag]}" for tag in ANNOT_TAGS if tag in annot)
        fields[8] = fields[8].rstrip() + extra
        out_lines[idx] = "\t".join(fields) + "\n"
        n_isoforms_annotated += 1

    for gid in gene_mrnas:
        if gene_annot.get(gid):
            n_genes_annotated += 1
        else:
            n_genes_no_annotation += 1

    with open(out_path, "w") as f:
        f.writelines(out_lines)

    with open(out_path + ".pending_proteomics_resolution.log", "w") as f:
        f.write("# Genes with >1 isoform where CD-HIT did NOT confirm the\n")
        f.write("# isoforms share one ORF (<98% identity/95% coverage between\n")
        f.write("# candidate ORFs). CDS left independently-called per isoform.\n")
        f.write("# Once proteomics peptide-to-protein-position results are back,\n")
        f.write("# check these first: single ORF with peptide support -> use it;\n")
        f.write("# distinct non-overlapping peptide support across candidates ->\n")
        f.write("# real evidence they're genuinely different proteins.\n")
        for gid in sorted(pending_proteomics):
            f.write(f"{gid}\t{','.join(gene_mrnas[gid])}\n")

    print(f"Total genes:                              {len(gene_mrnas)}")
    print(f"Genes with >1 isoform:                     {sum(1 for v in gene_mrnas.values() if len(v)>1)}")
    print(f"  CDS harmonized (CD-HIT confirmed):       {len(cdhit_harmonized_genes)}")
    print(f"  left unharmonized, pending proteomics:   {len(pending_proteomics)}")
    print(f"Isoform CDS records annotated:             {n_isoforms_annotated}")
    print(f"Genes with any product/ortholog annotation:{n_genes_annotated}")
    print(f"Genes with no annotation from student file: {n_genes_no_annotation}")
    print(f"Wrote: {out_path}")
    print(f"Wrote: {out_path}.pending_proteomics_resolution.log ({len(pending_proteomics)} genes)")


if __name__ == "__main__":
    main()
