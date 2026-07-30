# DGE Analysis: Swimming vs Adherent Samples
- hisat2 2.2.1
- edgeR 4.10.1
- samtools 1.20
- topGO 2.64.0

Aligned reads to genome using the following commands:
```
# Create index
hisat2-build -p 8 CbWHA1.fasta Cbombi_index

# Align reads
hisat2 -x Cbombi_index -1 reads/read1.fastq.gz -2 reads/read2.fastq.gz -S alignment.sam -p 12 --summary-file "hisat_summary.txt"

# Convert SAM to BAM
samtools view -S -b alignment.sam > alignment_unsorted.bam

# Sort BAM file
samtools sort -o alignment_sorted.bam alignment_unsorted.bam

# Remove duplicates
java -Xmx44g -jar picard-2.9.4/picard.jar MarkDuplicates INPUT=alignment_sorted.bam OUTPUT=alignment_dup.bam METRICS_FILE=dedup_metrics.txt REMOVE_DUPLICATES=true ASSUME_SORTED=true VALIDATION_STRINGENCY=SILENT

# Keep only primary alignments
samtools view -F 2308 -b alignment_dup.bam > alignment_primary.bam
```
Then ran count-reads.py to get the raw counts.
Raw counts were used as input to cpm.py, which creates a normalized counts matrix.
This matrix is used as the input to EdgeR in DE-Cbombi.R

GO analysis was performed using edgeR results obtained from DE-Cbombi.R. Running GO-Cbombi.R produced a list of the top GO terms from the overexpressed genes in swimming or adherent samples, and which genes those GO terms are associated with. It also generated a bar plot showing the number of overexpressed C. bombi genes associated with each of the top GO terms, as well as a corresponding plot with fold enrichment on the x-axis.
