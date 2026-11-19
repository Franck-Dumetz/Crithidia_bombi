"""Step 1 - convert the FragPipe/MaxLFQ report to a plain table and record what was dropped.

Input : the proteomics workbook (one sheet, one row per searched protein/ORF)
Output: data/prot_raw.tsv, used by every later step

Two things are removed here and nowhere else, so the decision is visible:

  * contaminant entries (Description contains "__CONT")
  * protein groups whose accession does not begin "CbWHA1". In this dataset those are
    566 entries carrying "Cb.1.*" accessions that match no gene in the CbWHA1
    annotation and appear to come from a second strain's gene set in the search
    database. They are far more often absent from WHA1 samples (25.8% have zero
    MaxLFQ in all six, against 1.5% of CbWHA1-named entries).

Note that the gene accession is NOT unique: a handful of genes were searched as two
ORFs (e.g. CbWHA1.01.000700.1:213-1466(+) and .3:292-609(+)). The full "Protein"
string is the unique key and is carried through as `orf`.

Usage:  python3 01_prepare_inputs.py <workbook.xlsx> [sheet]
"""
import os
import sys

import numpy as np
import pandas as pd

DATA = os.environ.get("CB_DATA", os.path.join(os.path.dirname(__file__), "..", "data")) + os.sep

xlsx = sys.argv[1]
sheet = sys.argv[2] if len(sys.argv) > 2 else 0
d = pd.read_excel(xlsx, sheet_name=sheet)
print("workbook: %d rows, %d columns" % d.shape)

cont = d.Description.astype(str).str.contains("__CONT", na=False)
d = d[~cont]
print("  removed %d contaminant entries" % cont.sum())

d["Protein ID"] = d["Protein ID"].astype(str).str.strip()
d["Protein"] = d["Protein"].astype(str).str.strip()
other = ~d["Protein ID"].str.startswith("CbWHA1")
if other.any():
    W = [c for c in d.columns if c.startswith("WHA1") and "MaxLFQ" in c]
    w = d[W].astype(float)
    print("  removed %d non-CbWHA1 entries (%.1f%% of them have zero MaxLFQ in all "
          "WHA1 samples, vs %.1f%% of the CbWHA1 entries)"
          % (other.sum(), 100 * (w[other].sum(1) == 0).mean(),
             100 * (w[~other].sum(1) == 0).mean()))
d = d[~other]

assert d["Protein"].is_unique, "the ORF-level 'Protein' column should be unique"
dup = d["Protein ID"].duplicated().sum()
print("  kept %d protein groups (%d gene accessions appear more than once, as separate ORFs)"
      % (len(d), dup))

os.makedirs(DATA, exist_ok=True)
d.to_csv(DATA + "prot_raw.tsv", sep="\t", index=False)
print("wrote %sprot_raw.tsv" % DATA)
