#!/usr/bin/env python3
"""
Harmonize CDS coordinates across isoforms that CD-HIT identified as
encoding (near-)identical proteins, using the genome-mapped TransDecoder
GFF3 (post cdna_alignment_orf_to_genome_orf.pl) + the CD-HIT .clstr file.

Why this is needed: TransDecoder.Predict --single_best_only picks one ORF
PER TRANSCRIPT, so different assembled isoforms of the same gene can end
up with slightly different CDS start/stop calls even when they really
encode the same protein (assembly/prediction noise in transcript
boundaries). CD-HIT clustering on the .pep file tells you which ORFs are
essentially the same protein; this script uses that to force all same-gene isoforms in a
cluster onto the CDS coordinates of CD-HIT's chosen representative.


This dataset is effectively intron-less (exon count == CDS count == mRNA
count in the source GFF3), so each mRNA has exactly one CDS interval --
the script assumes this and will flag (not silently mishandle) any mRNA
with more than one CDS line.

Usage:
    python3 harmonize_isoform_CDS.py \\
        Cb_2.5-8_annotation.transdecoder.genome.gff3 \\
        transdecoder_collapsed.pep.clstr \\
        Cb_2.5-8_annotation.transdecoder.genome.CDSharmonized.gff3
"""
import re
import sys
from collections import defaultdict

ID_RE = re.compile(r"ID=([^;]+)")
PARENT_RE = re.compile(r"Parent=([^;]+)")


def parse_clstr(path):
    """Return list of clusters; each cluster is a list of (mrna_id, is_rep)."""
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
    """Parse the GFF3 into per-mRNA records, preserving raw lines/order."""
    lines = []
    mrna_gene = {}          # mrna_id -> gene_id (from mRNA's Parent=)
    mrna_lines = defaultdict(dict)  # mrna_id -> {feature_type: [line_idx,...]}
    with open(path) as f:
        for idx, raw in enumerate(f):
            lines.append(raw)
            if raw.startswith("#") or not raw.strip():
                continue
            fields = raw.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue
            ftype = fields[2]
            attrs = fields[8]
            if ftype == "mRNA":
                idm = ID_RE.search(attrs)
                parentm = PARENT_RE.search(attrs)
                if idm and parentm:
                    mrna_gene[idm.group(1)] = parentm.group(1)
            elif ftype in ("CDS", "exon", "five_prime_UTR", "three_prime_UTR"):
                parentm = PARENT_RE.search(attrs)
                if parentm:
                    mrna_lines[parentm.group(1)].setdefault(ftype, []).append(idx)
    return lines, mrna_gene, mrna_lines


def field(line, i):
    return line.rstrip("\n").split("\t")[i]


def set_start_end(line, start, end):
    parts = line.rstrip("\n").split("\t")
    parts[3] = str(start)
    parts[4] = str(end)
    return "\t".join(parts) + "\n"


def main():
    if len(sys.argv) != 4:
        sys.exit(__doc__)
    gff_path, clstr_path, out_path = sys.argv[1:4]

    lines, mrna_gene, mrna_lines = load_gff(gff_path)
    clusters = parse_clstr(clstr_path)

    n_multi = 0
    n_same_gene = 0
    n_mixed_gene = 0
    n_harmonized_isoforms = 0
    n_skipped_multi_cds = 0
    n_skipped_missing = 0
    mixed_clusters_log = []
    nofit_log = []

    for cluster in clusters:
        if len(cluster) < 2:
            continue
        n_multi += 1

        # drop members not present in this GFF3 (e.g. filtered elsewhere)
        cluster = [(mid, rep) for mid, rep in cluster if mid in mrna_gene]
        if len(cluster) < 2:
            continue

        gene_ids = {mrna_gene[mid] for mid, _ in cluster}
        if len(gene_ids) > 1:
            n_mixed_gene += 1
            mixed_clusters_log.append(cluster)
            continue
        n_same_gene += 1

        rep_id = next((mid for mid, rep in cluster if rep), None)
        if rep_id is None:
            rep_id = cluster[0][0]  # fallback, shouldn't normally happen

        rep_cds_idx = mrna_lines[rep_id].get("CDS", [])
        if len(rep_cds_idx) != 1:
            n_skipped_multi_cds += 1
            continue
        rep_cds_line = lines[rep_cds_idx[0]]
        rep_start, rep_end = int(field(rep_cds_line, 3)), int(field(rep_cds_line, 4))
        rep_phase = field(rep_cds_line, 7)
        rep_strand = field(rep_cds_line, 6)

        for mid, is_rep in cluster:
            if is_rep:
                continue
            cds_idx = mrna_lines[mid].get("CDS", [])
            exon_idx = mrna_lines[mid].get("exon", [])
            if len(cds_idx) != 1 or len(exon_idx) != 1:
                n_skipped_multi_cds += 1
                continue

            exon_line = lines[exon_idx[0]]
            exon_start, exon_end = int(field(exon_line, 3)), int(field(exon_line, 4))
            if not (exon_start <= rep_start and rep_end <= exon_end):
                # representative's CDS doesn't fit inside this isoform's own
                # transcript span -- don't force it, flag for manual review
                n_skipped_missing += 1
                nofit_log.append(
                    f"{mid}: exon {exon_start}-{exon_end} does not contain "
                    f"rep {rep_id} CDS {rep_start}-{rep_end}"
                )
                continue

            # overwrite this isoform's CDS coordinates + phase to match rep
            old_cds_line = lines[cds_idx[0]]
            new_cds_line = set_start_end(old_cds_line, rep_start, rep_end)
            parts = new_cds_line.rstrip("\n").split("\t")
            parts[7] = rep_phase
            if "cds_harmonized_from=" not in parts[8]:
                parts[8] = parts[8].rstrip() + f";cds_harmonized_from={rep_id}"
            lines[cds_idx[0]] = "\t".join(parts) + "\n"

            # recompute UTRs against the new CDS boundaries within this
            # isoform's own (unchanged) exon span
            strand = field(exon_line, 6)
            utr5_idx = mrna_lines[mid].get("five_prime_UTR", [])
            utr3_idx = mrna_lines[mid].get("three_prime_UTR", [])

            if strand == "+":
                new5 = (exon_start, rep_start - 1)
                new3 = (rep_end + 1, exon_end)
            else:
                new5 = (rep_end + 1, exon_end)
                new3 = (exon_start, rep_start - 1)

            for utr_idx, (s, e) in ((utr5_idx, new5), (utr3_idx, new3)):
                if utr_idx and s <= e:
                    lines[utr_idx[0]] = set_start_end(lines[utr_idx[0]], s, e)
                elif utr_idx and s > e:
                    # new CDS consumes the whole exon on this side -- no UTR
                    lines[utr_idx[0]] = None  # mark for removal

            n_harmonized_isoforms += 1

    with open(out_path, "w") as f:
        for line in lines:
            if line is not None:
                f.write(line)

    print(f"Multi-member CD-HIT clusters:              {n_multi}")
    print(f"  same-gene (harmonized):                   {n_same_gene}")
    print(f"  mixed-gene (left untouched, likely paralogs): {n_mixed_gene}")
    print(f"Isoform CDS records rewritten:              {n_harmonized_isoforms}")
    print(f"Skipped (multi-CDS/exon mRNA, needs review): {n_skipped_multi_cds}")
    print(f"Skipped (rep CDS doesn't fit isoform exon):  {n_skipped_missing}")
    print(f"Wrote: {out_path}")

    with open(out_path + ".mixed_gene_clusters.log", "w") as f:
        f.write("# CD-HIT clusters spanning >1 gene ID -- NOT harmonized, review manually\n")
        f.write("# (likely near-identical paralogs at different loci)\n")
        for cluster in mixed_clusters_log:
            f.write(", ".join(f"{mid}{'*' if rep else ''}" for mid, rep in cluster) + "\n")
    print(f"Wrote: {out_path}.mixed_gene_clusters.log ({len(mixed_clusters_log)} clusters)")

    with open(out_path + ".nofit_isoforms.log", "w") as f:
        f.write("# Isoforms where the representative's CDS didn't fit inside\n")
        f.write("# this isoform's own exon span -- NOT harmonized, review manually\n")
        for entry in nofit_log:
            f.write(entry + "\n")
    print(f"Wrote: {out_path}.nofit_isoforms.log ({len(nofit_log)} entries)")


if __name__ == "__main__":
    main()
