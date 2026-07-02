# DGE Analysis: Swimming vs Adherent Samples
- bwa 0.7.17-r1188
- samtools 1.20
- DESeq2 1.52.0
- edgeR 4.10.1

Aligned reads to genome using:
```
bwa mem -t 24 genome.fasta sample.fastq > sample.sam
```
Converted to bam using samtools, then ran count-reads.py to get the raw counts.
Raw counts were used as input to DE-Cbombi.R

GO analysis was performed using edgeR results obtained from DE-Cbombi.R. Running GO-Cbombi.R produced a list of the top GO terms from the overexpressed genes in swimming or adherent samples, and which genes those GO terms are associated with. It also generated a bar plot showing the number of overexpressed C. bombi genes associated with each of the top GO terms, as well as a corresponding plot with fold enrichment on the x-axis.
