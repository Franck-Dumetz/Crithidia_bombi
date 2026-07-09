# This script takes one of the output files of cnvkit and creates a new file with only
# significant log2 values (log2 > 1.0 and log2 < -1.0).


with open("Cbombi_swiss.cns", "r") as f, open("Cbombi_swiss_filt.cns", "w") as out:
    i = 0
    for line in f:
        if i == 0:
            i = 1
            continue
        fields = line.split('\t')
        log2 = float(fields[4])
        if log2 <= -1.0 or log2 >= 1.0:
            out.write(line)
