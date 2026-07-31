#!/usr/bin/env python3
"""
Remove transcripts (and their exon lines) from a GTF, given a CSV list of transcript IDs.

CSV format: one transcript_id per row, no header, e.g.:
    CbWHA1.31.003420.1
    CbWHA1.31.003440.2
    ...

Usage:
    python remove_transcripts_csv.py --gtf input.gtf --csv to_remove.csv --out output.gtf
"""
import argparse
import csv
import re


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--gtf", required=True, help="Input GTF file")
    ap.add_argument("--csv", required=True, help="CSV file with one transcript_id per row (no header)")
    ap.add_argument("--out", required=True, help="Output GTF file")
    ap.add_argument("--report", default=None, help="Optional path for a text summary report")
    args = ap.parse_args()

    # 1. Load removal list
    targets = set()
    with open(args.csv, newline="") as f:
        for row in csv.reader(f):
            if not row or not row[0].strip():
                continue
            targets.add(row[0].strip())

    # 2. Find which target IDs actually exist as transcript_id in the GTF
    tx_ids = set()
    with open(args.gtf) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != "transcript":
                continue
            m = re.search(r'transcript_id "([^"]+)"', fields[8])
            if m:
                tx_ids.add(m.group(1))

    matched = targets & tx_ids
    unmatched = sorted(targets - tx_ids)

    # 3. Rewrite GTF, dropping transcript + exon lines for matched IDs
    out_lines = []
    removed_counts = {}
    with open(args.gtf) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                out_lines.append(line)
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] not in ("transcript", "exon"):
                out_lines.append(line)
                continue
            m = re.search(r'transcript_id "([^"]+)"', fields[8])
            tid = m.group(1) if m else None
            if tid in matched:
                removed_counts[fields[2]] = removed_counts.get(fields[2], 0) + 1
                continue
            out_lines.append(line)

    with open(args.out, "w") as f:
        f.writelines(out_lines)

    report_lines = [
        f"Removal list: {len(targets)} unique IDs",
        f"Matched and removed: {len(matched)} transcripts",
    ]
    for feat, n in sorted(removed_counts.items()):
        report_lines.append(f"  {feat} lines removed: {n}")
    report_lines.append(f"\nUnmatched (not found as a transcript_id in the GTF): {len(unmatched)}")
    for x in unmatched:
        report_lines.append(f"  {x}")
    report_text = "\n".join(report_lines)

    if args.report:
        with open(args.report, "w") as f:
            f.write(report_text + "\n")

    print(report_text)


if __name__ == "__main__":
    main()
