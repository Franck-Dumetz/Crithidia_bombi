#!/usr/bin/env python3
"""
Correct transcript start/end coordinates in a GTF, given a CSV of corrections.

CSV format: transcript_id, new_start, new_end  -- no header. Leave a field empty to
leave that side unchanged (only correcting start or only end is fine), e.g.:
    CbWHA1.01.000440.1,187914,
    CbWHA1.31.003380.1,,1285779
    CbWHA1.31.003650.1,1369505,1371439

Only single-exon transcripts are supported (exon coords are kept identical to the
transcript's, matching how this GTF was built). For each row this script:
  - skips rows with no values given
  - skips rows where the ID isn't found as a transcript_id in the GTF (reported as unmatched)
  - skips rows where the new coordinate would put start >= end relative to the
    transcript's current bounds (reported as an anomaly, NOT applied)
  - skips rows where the same ID appears more than once with different values
    (reported as a conflict, NOT applied)
  - otherwise applies the correction to both the transcript and exon line

Usage:
    python correct_transcripts_csv.py --gtf input.gtf --csv to_correct.csv --out output.gtf
"""
import argparse
import csv
import re
from collections import defaultdict


def parse_coord(val):
    val = val.strip() if val is not None else ""
    if val == "":
        return None
    return int(val)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--gtf", required=True, help="Input GTF file")
    ap.add_argument("--csv", required=True, help="CSV file: transcript_id,new_start,new_end (no header)")
    ap.add_argument("--out", required=True, help="Output GTF file")
    ap.add_argument("--report", default=None, help="Optional path for a text summary report")
    args = ap.parse_args()

    # 1. Load raw correction rows
    raw_rows = []
    with open(args.csv, newline="") as f:
        for row in csv.reader(f):
            if not row or not row[0].strip():
                continue
            sid = row[0].strip()
            b = parse_coord(row[1]) if len(row) > 1 else None
            c = parse_coord(row[2]) if len(row) > 2 else None
            raw_rows.append((sid, b, c))

    # 2. Drop "no values given" rows, then dedup / flag conflicting duplicates
    rows = [r for r in raw_rows if not (r[1] is None and r[2] is None)]
    by_id = defaultdict(set)
    for sid, b, c in rows:
        by_id[sid].add((b, c))
    conflicts = {sid: vals for sid, vals in by_id.items() if len(vals) > 1}
    rows_dedup = [(sid, list(vals)[0][0], list(vals)[0][1])
                  for sid, vals in by_id.items() if sid not in conflicts]

    # 3. Load current GTF transcript coords
    tx_current = {}
    with open(args.gtf) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != "transcript":
                continue
            m = re.search(r'transcript_id "([^"]+)"', fields[8])
            if m:
                tx_current[m.group(1)] = (int(fields[3]), int(fields[4]))

    # 4. Validate each correction against current bounds
    corrections = {}
    unmatched = []
    anomalies = []
    for sid, b, c in rows_dedup:
        if sid not in tx_current:
            unmatched.append((sid, b, c))
            continue
        cur_s, cur_e = tx_current[sid]
        problem = None
        if b is not None and c is not None and b >= c:
            problem = "new_start>=new_end"
        elif b is not None and c is None and b >= cur_e:
            problem = "new_start>=current_end"
        elif c is not None and b is None and c <= cur_s:
            problem = "new_end<=current_start"
        if problem:
            anomalies.append((sid, b, c, cur_s, cur_e, problem))
            continue
        corrections[sid] = (b, c)

    # 5. Rewrite GTF: apply corrections to transcript + exon lines
    out_lines = []
    applied = 0
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
            if tid in corrections:
                new_s, new_e = corrections[tid]
                if new_s is not None:
                    fields[3] = str(new_s)
                if new_e is not None:
                    fields[4] = str(new_e)
                out_lines.append("\t".join(fields) + "\n")
                if fields[2] == "transcript":
                    applied += 1
            else:
                out_lines.append(line)

    with open(args.out, "w") as f:
        f.writelines(out_lines)

    report_lines = [
        f"Correction list: {len(raw_rows)} rows",
        f"Applied corrections: {applied} transcripts (transcript + exon line updated for each)",
        f"Rows with no start/end value given: {len(raw_rows) - len(rows)}",
        f"\nUnmatched (not found as transcript_id): {len(unmatched)}",
    ]
    for x in unmatched:
        report_lines.append(f"  {x}")
    report_lines.append(f"\nAnomalies (new coord inconsistent with current bounds, NOT applied): {len(anomalies)}")
    for x in anomalies:
        report_lines.append(f"  {x}")
    if conflicts:
        report_lines.append(f"\nConflicting duplicate rows (same ID, different values, NOT applied): {len(conflicts)}")
        for sid, vals in conflicts.items():
            report_lines.append(f"  {sid}: {vals}")
    report_text = "\n".join(report_lines)

    if args.report:
        with open(args.report, "w") as f:
            f.write(report_text + "\n")

    print(report_text)


if __name__ == "__main__":
    main()
