# This script first creates a list of all the genes present in the genome annotation file for
# C.bombi. After that, it looks through the output of the bedtools intersect file and calculates
# the log2 for each gene. It outputs a file with the logs of all the genes and another file with
# the CNVs, indicating whether they are losses or gains.

import subprocess

all_genes = []
missing_genes = []

with open("../../../fastas-ortho/Cb_Cfannotated_genome_transdecoder_6-22.gff3", "r") as f:
    for line in f:
        fields = line.split("\t")
        if len(fields) != 9:
            continue
        if fields[2] != "CDS":
            continue
        all_genes.append(fields[8].split("Parent=")[1].split(";")[0].split("\n")[0])

with open("all_logs.txt", "w") as all_logs, open("cnvs_only.txt", "w") as sig:
    for gene in all_genes:
        cmd = ["grep", f"{gene}", "gene_overlap_CDS.txt"]
        try:
            result = subprocess.run(cmd, stdout=subprocess.PIPE, universal_newlines=True, check=True)
        except subprocess.CalledProcessError as e:
            missing_genes.append(gene)
            continue
        lines = result.stdout.strip().split("\n")
        num = 0
        denom = 0
        for line in lines:
            fields = line.split("\t")
            log = float(fields[3])
            bin_start = int(fields[1])
            bin_end = int(fields[2])
            w = float(fields[4])
            gene_start = int(fields[8])
            gene_end = int(fields[9])
            # calculate percentage of gene that is covered by bin
            length = gene_end - gene_start
            overlap_end = gene_end
            overlap_start = gene_start
            if gene_end > bin_end: # if the gene ends after the bin ends
                overlap_end = bin_end
            if gene_start < bin_start: # if the gene starts before the bin starts
                overlap_start = bin_start
            percent = (overlap_end - overlap_start)/length
            num += log * percent * w
            denom += percent * w
        final_log = num/denom
        print(final_log)
        all_logs.write(f"{gene}\t{final_log}\n")
        if final_log >= 1.0:
            sig.write(f"{gene}\t{final_log}\tgain\n")
        if final_log <= -1.0:
            sig.write(f"{gene}\t{final_log}\tloss\n")

print(missing_genes)
