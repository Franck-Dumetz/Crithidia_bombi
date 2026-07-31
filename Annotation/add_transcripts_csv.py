#!/usr/bin/env python3
"""
Add new transcripts to a GTF, given a CSV of new entries.

CSV format: transcript_id, chrom, start, end, strand  -- no header, e.g.:
    CbWHA1.31.003480.2,CbWHA1_chr31,1321417,1322613,-
    CbWHA1.31.100000.1,CbWHA1_chr31,1350022,1351446,-

transcript_id is expected to look like "<gene_id>.<isoform_number>" (gene_id being
everything before the last dot). If a given ID already exists in the GTF, this
script will NOT overwrite it -- instead it renames the new entry to the next free
isoform number for that gene (e.g. .1 -> .2) and reports the rename, so you never
silently lose or clobber an existing transcript.

Only single-exon transcripts are produced (matching how this GTF was built): each
new entry gets one "transcript" line and one "exon" line with identical coordinates.
New lines are written with source "manual_annotation" and no cov/FPKM/TPM values.

Usage:
    python add_transcripts_csv.py --gtf input.gtf --csv new_entries.csv --out output.gtf
"""
import argparse
import csv
import re
from collections import defaultdict


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--gtf", required=True, help="Input GTF file")
    ap.add_argument("--csv", required=True, help="CSV file: transcript_id,chrom,start,end,strand (no header)")
    ap.add_argument("--out", required=True, help="Output GTF file")
    ap.add_argument("--report", default=None, help="Optional path for a text summary report")
    args = ap.parse_args()

    # 1. Load new entries
    rows = []
    with open(args.csv, newline="") as f:
        for row in csv.reader(f):
            if not row or not row[0].strip():
                continue
            sid, chrom, start, end, strand = [x.strip() for x in row[:5]]
            rows.append((sid, chrom, int(start), int(end), strand))

    # 2. Load current GTF state: existing transcript IDs + per-gene max isoform number
    existing_tx = set()
    gene_max_iso = defaultdict(int)
    gtf_lines = []
    with open(args.gtf) as fh:
        for line in fh:
            gtf_lines.append(line)
            if line.startswith("#") or not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != "transcript":
                continue
            attrs = fields[8]
            mg = re.search(r'gene_id "([^"]+)"', attrs)
            mt = re.search(r'transcript_id "([^"]+)"', attrs)
            if not (mg and mt):
                continue
            gid, tid = mg.group(1), mt.group(1)
            existing_tx.add(tid)
            m = re.match(re.escape(gid) + r"\.(\d+)$", tid)
            if m:
                gene_max_iso[gid] = max(gene_max_iso[gid], int(m.group(1)))

    # 3. Build add-list, renaming on collision with an existing transcript_id
    to_add = []
    renamed = []
    for sid, chrom, start, end, strand in rows:
        gene_id = sid.rsplit(".", 1)[0]
        final_id = sid
        if sid in existing_tx:
            gene_max_iso[gene_id] += 1
            final_id = f"{gene_id}.{gene_max_iso[gene_id]}"
            while final_id in existing_tx:
                gene_max_iso[gene_id] += 1
                final_id = f"{gene_id}.{gene_max_iso[gene_id]}"
            renamed.append((sid, final_id, chrom, start, end, strand))
        existing_tx.add(final_id)
        to_add.append((final_id, gene_id, chrom, start, end, strand))

    # 4. Write new GTF: original content + new transcript/exon lines appended
    new_lines = []
    for final_id, gene_id, chrom, start, end, strand in to_add:
        attrs_tx = f'gene_id "{gene_id}"; transcript_id "{final_id}";'
        attrs_ex = f'gene_id "{gene_id}"; transcript_id "{final_id}"; exon_number "1";'
        new_lines.append(f"{chrom}\tmanual_annotation\ttranscript\t{start}\t{end}\t.\t{strand}\t.\t{attrs_tx}\n")
        new_lines.append(f"{chrom}\tmanual_annotation\texon\t{start}\t{end}\t.\t{strand}\t.\t{attrs_ex}\n")

    with open(args.out, "w") as f:
        f.writelines(gtf_lines)
        f.writelines(new_lines)

    report_lines = [
        f"New entries in CSV: {len(rows)}",
        f"Added: {len(to_add)} new transcripts (transcript + exon line each)",
        f"Renamed due to ID collision with an existing transcript: {len(renamed)}",
    ]
    for orig, final, chrom, start, end, strand in renamed:
        report_lines.append(f"  {orig} -> {final}  ({chrom}:{start}-{end}{strand})")
    report_text = "\n".join(report_lines)

    if args.report:
        with open(args.report, "w") as f:
            f.write(report_text + "\n")

    print(report_text)


if __name__ == "__main__":
    main()
