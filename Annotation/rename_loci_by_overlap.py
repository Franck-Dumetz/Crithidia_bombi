#!/usr/bin/env python3
"""
Rename transcripts in a GTF to a "<prefix>.<chrom_label>.<locus_number>.<isoform>"
scheme, where loci are defined by TRUE genomic overlap rather than trusting whatever
gene_id the file already has.

What it does, in one pass:
  1. Groups transcripts into loci: same chromosome + same strand, and genomically
     overlapping (touching/adjacent does NOT count as overlap -- there must be a
     real base-pair overlap). This is the important part: it does not matter what
     gene_id a transcript already has -- if two transcripts don't actually overlap,
     they will never end up in the same new locus, and if they do overlap they
     will always end up in the same locus, even if their old gene_ids differed.
  2. Within each locus, drops any transcript that is an exact coordinate duplicate
     of another one in the same locus (same start/end) -- keeping the first one
     seen. This is reported so you can see what was collapsed.
  3. Numbers loci along each chromosome (in increments of --step, starting at
     --start) in order of their leftmost coordinate, and numbers surviving
     isoforms within a locus 1, 2, 3... in order of their start coordinate.

Only single-exon-per-transcript GTFs are supported (i.e. the "exon" feature's
coordinates are always identical to its parent "transcript" feature's coordinates
-- this is what StringTie produces for these kinetoplastid genomes). The script
will refuse to run (and tell you) if it finds a transcript with more than one exon,
since the renumbering logic here doesn't try to handle introns.

Chromosome labels:
  By default, a chromosome named like "..._chr7" or "..._chr07" becomes a
  zero-padded number ("07"), and anything else (e.g. a contig or the mitochondrial
  genome) becomes whatever comes after the last underscore, e.g.
  "CbWHA1_contig_a" -> "a", "CbWHA1_maxicircle" -> "maxicircle".
  If you want different/shorter labels for non-numbered sequences, pass
  --contig-labels a CSV file of "chromosome_name,label" pairs (no header) to
  override the auto-derived label for those specific sequences.

Usage:
    python rename_loci_by_overlap.py --gtf input.gtf --out output.gtf \
        --map correspondence_key.csv

    # custom labels for non-chr sequences, custom prefix/step:
    python rename_loci_by_overlap.py --gtf input.gtf --out output.gtf \
        --contig-labels my_labels.csv --prefix CbWHA1 --start 10 --step 10
"""
import argparse
import csv
import re
from collections import defaultdict


def build_chrom_label_fn(gtf_path, chr_pattern, pad_width, overrides):
    # figure out how many digits to zero-pad numbered chromosomes to, if not given explicitly
    numbers = []
    chroms_seen = set()
    with open(gtf_path) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 1:
                continue
            chroms_seen.add(fields[0].strip())
    for c in chroms_seen:
        m = re.search(chr_pattern, c)
        if m:
            numbers.append(m.group(1))
    auto_width = max((len(n) for n in numbers), default=2)
    width = pad_width if pad_width else auto_width

    def label_fn(chrom):
        chrom = chrom.strip()
        if chrom in overrides:
            return overrides[chrom]
        m = re.search(chr_pattern, chrom)
        if m:
            return m.group(1).zfill(width)
        return chrom.rsplit("_", 1)[-1]

    return label_fn


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("--gtf", required=True, help="Input GTF file")
    ap.add_argument("--out", required=True, help="Output GTF file")
    ap.add_argument("--map", default=None, help="Optional path for a CSV correspondence key: old_transcript_id,new_transcript_id")
    ap.add_argument("--report", default=None, help="Optional path for a text summary report")
    ap.add_argument("--prefix", default="CbWHA1", help="ID prefix (default: CbWHA1)")
    ap.add_argument("--start", type=int, default=10, help="First locus number per chromosome (default: 10)")
    ap.add_argument("--step", type=int, default=10, help="Increment between locus numbers (default: 10)")
    ap.add_argument("--width", type=int, default=6, help="Zero-padded width of the locus number (default: 6)")
    ap.add_argument("--chr-pattern", default=r"chr0*(\d+)$",
                     help=r'Regex (with one capture group) used to pull the chromosome number out of a '
                          r'sequence name, e.g. default r"chr0*(\d+)$" matches "..._chr07" -> "7"')
    ap.add_argument("--chr-pad-width", type=int, default=None,
                     help="Force a specific zero-pad width for chromosome numbers (default: auto-detected)")
    ap.add_argument("--contig-labels", default=None,
                     help="Optional CSV file of chromosome_name,label (no header) to override "
                          "auto-derived labels for specific sequences (e.g. contigs, mitochondrial genome)")
    args = ap.parse_args()

    overrides = {}
    if args.contig_labels:
        with open(args.contig_labels, newline="") as f:
            for row in csv.reader(f):
                if not row or not row[0].strip():
                    continue
                overrides[row[0].strip()] = row[1].strip()

    chrom_label = build_chrom_label_fn(args.gtf, args.chr_pattern, args.chr_pad_width, overrides)

    # 1. Parse transcripts, and fail loudly if any transcript has more than one exon
    records = []  # [chrom, start, end, strand, transcript_id]
    exon_count = defaultdict(int)
    with open(args.gtf) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue
            attrs = fields[8]
            m = re.search(r'transcript_id "([^"]+)"', attrs)
            if not m:
                continue
            tid = m.group(1)
            if fields[2] == "transcript":
                records.append([fields[0].strip(), int(fields[3]), int(fields[4]), fields[6], tid])
            elif fields[2] == "exon":
                exon_count[tid] += 1

    multiexon = [t for t, n in exon_count.items() if n > 1]
    if multiexon:
        raise SystemExit(
            f"ERROR: {len(multiexon)} transcript(s) have more than one exon, e.g. {multiexon[0]}. "
            "This script only supports single-exon-per-transcript GTFs. Aborting without writing output."
        )

    # 2. Cluster by (chrom, strand) via true genomic overlap
    by_chrom_strand = defaultdict(list)
    for r in records:
        by_chrom_strand[(r[0], r[3])].append(r)

    clusters_by_chrom = defaultdict(list)
    for (chrom, strand), recs in by_chrom_strand.items():
        recs_sorted = sorted(recs, key=lambda x: x[1])
        clusters = []
        cur_cluster = [recs_sorted[0]]
        cur_max_end = recs_sorted[0][2]
        for r in recs_sorted[1:]:
            if r[1] <= cur_max_end:
                cur_cluster.append(r)
                cur_max_end = max(cur_max_end, r[2])
            else:
                clusters.append(cur_cluster)
                cur_cluster = [r]
                cur_max_end = r[2]
        clusters.append(cur_cluster)
        for cl in clusters:
            min_start = min(x[1] for x in cl)
            clusters_by_chrom[chrom].append((min_start, cl))

    # 3. Within each cluster: drop exact-coordinate duplicates, number the locus + survivors
    tid_new_tid = {}
    removed_dup_tids = {}

    for chrom, cluster_list in clusters_by_chrom.items():
        label = chrom_label(chrom)
        cluster_list.sort(key=lambda x: x[0])
        locus_num = args.start
        for min_start, cl in cluster_list:
            new_gid = f"{args.prefix}.{label}.{locus_num:0{args.width}d}"
            locus_num += args.step
            cl_sorted = sorted(cl, key=lambda x: x[1])
            seen_coords = {}
            survivors = []
            for _, start, end, _, old_tid in cl_sorted:
                key = (start, end)
                if key in seen_coords:
                    removed_dup_tids[old_tid] = seen_coords[key]
                    continue
                seen_coords[key] = old_tid
                survivors.append(old_tid)
            for new_iso, old_tid in enumerate(survivors, start=1):
                tid_new_tid[old_tid] = f"{new_gid}.{new_iso}"

    # 4. Rewrite GTF with new gene_id / transcript_id, dropping removed duplicates
    out_lines = []
    with open(args.gtf) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                out_lines.append(line)
                continue
            fields = line.rstrip("\n").split("\t")
            fields[0] = fields[0].strip()
            attrs = fields[8]
            mt = re.search(r'transcript_id "([^"]+)"', attrs)
            if not mt:
                out_lines.append(line)
                continue
            tid = mt.group(1)
            if tid in removed_dup_tids:
                continue
            new_tid = tid_new_tid[tid]
            new_gid = new_tid.rsplit(".", 1)[0]
            mg = re.search(r'gene_id "([^"]+)"', attrs)
            old_gid = mg.group(1) if mg else None
            new_attrs = attrs.replace(f'transcript_id "{tid}"', f'transcript_id "{new_tid}"')
            if old_gid:
                new_attrs = new_attrs.replace(f'gene_id "{old_gid}"', f'gene_id "{new_gid}"')
            fields[8] = new_attrs
            out_lines.append("\t".join(fields) + "\n")

    with open(args.out, "w") as f:
        f.writelines(out_lines)

    # 5. Correspondence key + report
    if args.map:
        rows_out = [[old, new] for old, new in tid_new_tid.items()]
        for removed_tid, surviving_old in removed_dup_tids.items():
            rows_out.append([removed_tid, f"REMOVED (exact duplicate of {surviving_old} -> now {tid_new_tid[surviving_old]})"])
        rows_out.sort(key=lambda r: r[0])
        with open(args.map, "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["old_transcript_id", "new_transcript_id"])
            w.writerows(rows_out)

    n_loci = len(set(v.rsplit(".", 1)[0] for v in tid_new_tid.values()))
    report_text = (
        f"Input transcripts: {len(records)}\n"
        f"Exact-coordinate duplicates removed: {len(removed_dup_tids)}\n"
        f"Final transcript count: {len(tid_new_tid)}\n"
        f"Total loci after overlap-based clustering: {n_loci}\n"
    )
    if args.report:
        with open(args.report, "w") as f:
            f.write(report_text)
    print(report_text)


if __name__ == "__main__":
    main()
