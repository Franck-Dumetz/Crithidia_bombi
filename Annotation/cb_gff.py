# This script maps Cbombi TransDecoder genes to Cfasciculata orthologs using a TSV lookup, then parses a Cfasciculata FASTA header file
# to build gene product annotations in a dictionary. It finally updates a genome-aligned GFF3 by adding CF ortholog and gene product
# metadata to CDS features when a matching gene exists in the orthology mapping.


from collections import defaultdict
import subprocess

cb_cf = {}
cf_desc = {}

with open("CbombiCfasic_p1.tsv", "r") as table:
    for line in table:
        fields = line.strip().split('\t')
        cb_cf[fields[0]] = fields[1]

#with open("TriTrypDB-68_CfasciculataCfCl_AnnotatedProteins.fasta", "r") as fasta
cmd = ["grep", ">", "ortho-input/TriTrypDB-68_CfasciculataCfCl_AnnotatedProteins.fasta"]
lines = subprocess.run(cmd, stdout=subprocess.PIPE, universal_newlines=True, check=True)
lines = lines.stdout.strip().splitlines()
for line in lines:
    gene = line.split(">")[1].split(" ")[0] # before first space
    transcript = line.split("transcript=")[1].split(" ")[0]
    gene_product = line.split("gene_product=")[1].split(" |")[0]
    cf_desc[gene] = f"cf_ortholog={transcript};gene_product={gene_product}"


with open("Cb_genome_transdecoder_p1_6-12.gff3", "r") as inp, open("Cb_Cfannotated_genome_transdecoder_6-12.gff3", "w") as out:
    for line in inp:
        fields = line.strip().split('\t')
        #print(cb_cf[fields[0]])
        if fields[2] == "CDS":
            gene = fields[8].split("Parent=")[1]
            if gene in cb_cf.keys():
                fields[8] = f"{fields[8]};{cf_desc[cb_cf[gene]]}"
            else:
                print(gene)
        out.write("\t".join(fields) + "\n")
