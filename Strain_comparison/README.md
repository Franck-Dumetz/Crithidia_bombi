# Strain Comparison and Variant Analysis Between CbWHA1 and 08.076

## Software
- bwa 0.7.17-r1188
- bcftools 1.11
- snpEff 5.4c

## Alignment
Reads from Swiss strain 08.076 were obtained from SRA (Accession numbers ERR2214503, ERR2214504, ERR2214505). All reads were combined into one fastq file and aligned to our WHA1 assembled genome using the Burrows-Wheeler Alignment.
```
bwa mem -t 8 ../cbombi/CbWHA1_assembly.noBDEF.fasta fastqs/08706.fastq > alignment.sam
```
The SAM file was converted to a BAM file for the next steps.

## Variant Calling
The following command was used to get variants from the alignment:
```
bcftools mpileup -f ../cbombi/CbWHA1_assembly.noBDEF.fasta alignment.bam | bcftools call -mv -Oz -o variants.vcf.gz
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
*Add where the gff and fasta files were put
The following command was used to annotate the snps:
```
snpEff ann -v Cbombi swiss_indels.vcf.gz > swiss_indels_ann.vcf
```
