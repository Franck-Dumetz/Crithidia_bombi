# Strain Comparison and Variant Analysis Between CbWHA1 and 08.076

## Software
- bwa 0.7.17-r1188
- bcftools 1.11
- bedtools v2.27.1
- snpEff 5.4c
- cnvkit 0.9.3

## Alignment
Reads from Swiss strain 08.076 were obtained from SRA (Accession numbers ERR2214503, ERR2214504, ERR2214505). All reads were combined into one fastq file and aligned to our WHA1 assembled genome using the Burrows-Wheeler Alignment.
```
bwa mem -t 8 CbWHA1_genome.fasta fastqs/08706.fastq > alignment.sam
```
The SAM file was converted to a BAM file for the next steps.

## Variant Calling
The following command was used to get variants from the alignment:
```
bcftools mpileup -f CbWHA1_genome.fasta alignment.bam | bcftools call -mv -Oz -o variants.vcf.gz
```
Variants with a read depth < 25 and a quality score < 20 were filtered out using the following command:
```
bcftools filter -e 'QUAL<20 || DP<25' variants.vcf.gz -Oz -o variants_filtered.vcf.gz
```
Get the number of variants:
```
bcftools view -H variants_filt.vcf.gz | wc -l
# By the genotype
bcftools view -H -g hom variants_filt.vcf.gz | wc -l
bcftools view -H -g het variants_filt.vcf.gz | wc -l
```
Get the indels and snps:
```
bcftools view -v indels variants_filt.vcf.gz -Oz -o indels.vcf.gz
bcftools view -v snps variants_filt.vcf.gz -Oz -o snps.vcf.gz
```

## Variant Annotation
snpEff was downloaded using conda. A directory Cbombi was created in ~/.conda/envs/snpEff-env/share/snpeff-5.4.0c-0/data/, and the gff and genome fasta files were uploaded to that folder and renamed "genes.gff" and "sequences.fa". In the file snpEff.config, this line was added at the bottom:
```
Cbombi.genome : Crithidia_bombi
```
The following command was used to annotate the snps:
```
snpEff ann -v Cbombi snps.vcf.gz > snps_ann.vcf
```
To isolate the number of missense mutations, this command was used:
```
bcftools view -i 'INFO/ANN ~ "missense_variant"' snps_ann.vcf -Oz -o missense.vcf.gz
```
And the same command was used to isolate the other types of SNPs

The same steps were used to annotate the indels.

## CNVs
Before running cnvkit, WHA1 ONT reads were mapped to the WHA1 genome, and the bam file is used as a normal reference for cnvkit.

The following command was run using cnvkit:
```
cnvkit.py batch Cbombi_swiss.bam -n Cbombi_WHA1.bam -f CbWHA1_genome.fasta -m wgs -d cnvs
```
### Primary Analysis
In the output folder, there is a file Cbombi_swiss.cns containing genomic regions where adjacent bins were merged into high-confidence copy-number segments. To keep only the segments that represented true CNVs, we ran the script filter_cns.py. We used the log2 value for each segment to determine whether they were true CNVs, using log2 < -1.0 and log2 > 1.0 as cutoffs. 

To get further information on whether these genomic regions overlapped with any coding regions, the following command was run:
```
bedtools intersect -a Cbombi_swiss_filt.cns -b annotation.gff3 -wa -wb > cns_gene_overlap.txt
```

### Secondary Analysis
In the output folder, there is a file called Cbombi_swiss.cnr containing bin-level statistics, including log2, for the entire genome (bin size was set to 94bp). The same bedtools intersect command from above was run using the .cnr file, and the script cnr_cnvs.py was used to calculate the log2 for each gene and determine which ones are CNVs. The cnr file provides a log2 and a confidence value for each bin, so the log2 for each gene was calculated using the equation:

```math
gene\_log2 = \frac{\sum^{gene\_bins}log2 \times confidence \times percent}{\sum^{gene\_bins} confidence \times percent}
```

where gene_bins refers to the bins that overlap a gene and percent refers to how much of the gene each bin covers.

Before comparison with the primary analysis, genes that were identified as CNVs that has less than 50% coverage by the bins were removed from the secondary analysis results. Similarly, genes with less than 50% coverage by the CNV genomic regions were removed from the primary analysis results.

The results from both analyses were compared to find overlapping genes and find discrepancies, including whether genes were near the CNV calling thresholds, had low coverage by CNV bins, or differed due to other factors.
