#!/usr/bin/env python3
import argparse, re, sys
from typing import Dict, Set

ATTR_RE = re.compile(r'([A-Za-z0-9_]+)\s+"?([^";]+)"?')

def parse_attrs(attr: str) -> Dict[str, str]:
    out = {}
    for part in attr.strip().split(";"):
        part = part.strip()
        if not part:
            continue
        m = ATTR_RE.match(part)
        if m:
            out[m.group(1)] = m.group(2)
    return out

def main():
    ap = argparse.ArgumentParser(
        description="Filter StringTie GTF by TPM and cov thresholds; outputs only kept transcripts and their features."
    )
    ap.add_argument("-i", "--input", required=True, help="Input GTF (StringTie)")
    ap.add_argument("-o", "--output", required=True, help="Output GTF")
    ap.add_argument("--tpm", type=float, required=True, help="Minimum TPM (keep if TPM >= this)")
    ap.add_argument("--cov", type=float, required=True, help="Minimum coverage (keep if cov >= this)")
    ap.add_argument("--keep-missing", action="store_true",
                    help="If set, transcripts missing TPM or cov are KEPT; default is to DROP if missing.")
    ap.add_argument("--feature-keep", default="transcript,exon",
                    help="Comma-sep list of feature types to emit for kept transcripts (default: transcript,exon).")
    args = ap.parse_args()

    features_to_emit = set(x.strip() for x in args.feature_keep.split(",") if x.strip())
    kept_tids: Set[str] = set()
    kept_gids: Set[str] = set()

    # ---- Pass 1: decide which transcripts to keep ----
    n_tx_total = n_tx_kept = 0
    with open(args.input) as fin:
        for line in fin:
            if not line or line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9:
                continue
            if cols[2] != "transcript":
                continue
            n_tx_total += 1
            attrs = parse_attrs(cols[8])
            tid = attrs.get("transcript_id")
            gid = attrs.get("gene_id")
            # Parse metrics (if present)
            tpm = attrs.get("TPM")
            cov = attrs.get("cov")
            tpm_val = float(tpm) if tpm is not None and tpm.replace(".","",1).isdigit() else None
            try:
                cov_val = float(cov) if cov is not None else None
            except ValueError:
                cov_val = None

            if tpm_val is None or cov_val is None:
                keep = args.keep-missing
            else:
                keep = (tpm_val >= args.tpm) and (cov_val >= args.cov)

            if keep and tid:
                kept_tids.add(tid)
                if gid: kept_gids.add(gid)
                n_tx_kept += 1

    # ---- Pass 2: write out kept items ----
    w_tx = w_feat = 0
    with open(args.input) as fin, open(args.output, "w") as fout:
        for line in fin:
            if line.startswith("#") or line.strip() == "":
                fout.write(line)
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9:
                fout.write(line)  # pass through malformed lines
                continue
            feature = cols[2]
            attrs = parse_attrs(cols[8])
            tid = attrs.get("transcript_id")

            if feature == "transcript":
                if tid in kept_tids:
                    fout.write(line)
                    w_tx += 1
            elif feature in features_to_emit:
                if tid in kept_tids:
                    fout.write(line)
                    w_feat += 1
            else:
                # For other feature types without transcript_id, you could add logic keyed by gene_id if needed.
                pass

    sys.stderr.write(
        f"[filter_stringtie_gtf] transcripts: kept {n_tx_kept}/{n_tx_total}; "
        f"wrote transcripts={w_tx}, features={w_feat}; out={args.output}\n"
    )

if __name__ == "__main__":
    main()
