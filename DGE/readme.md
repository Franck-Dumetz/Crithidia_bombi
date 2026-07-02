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
