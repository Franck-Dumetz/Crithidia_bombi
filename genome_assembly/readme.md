# Crithidia bombi genome assembly

BioProject PRJNA1231061


Software requirements: <br />
- bedtools-2.27.1 <br />
- busco-5.4.3 <br />
- filtlong-0.2.0 <br />
- hifiasm-0.24.0 <br />
- minimap2-2.2.28 <br />
- mummer-3.23 <br />
- quast-5.2.0 <br />
- samtools-1.20 <br />
- seqkit-0.7.2 <br />
- seqtk-1.0-r63 <br />

## Preparing the data by filtering for reads longer than 1000 bp
```
/usr/local/packages/filtlong-0.2.0/bin/filtlong --min_length 1000 --keep_percent 90 EDS10_20240819_R84050_PL14293-001_1-1-A01_bc2093-bc2093.hifi_reads.fastq.gz | gzip > Cbombi_hifi_filtered.fastq.gz

zgrep '^@' EDS10_20240819_R84050_PL14293-001_1-1-A01_bc2093-bc2093.hifi_reads.fastq.gz > original_read_names.txt

zgrep '^@' Cbombi_hifi_filtered.fastq.gz > filtered_read_names.txt
```
Find the Read Names that Were Filtered Out
```
grep -v -f filtered_read_names.txt original_read_names.txt > filtered_out_read_names.txt
```
## Extract reads and create a new fastq
```
seqtk subseq EDS10_20240819_R84050_PL14293-001_1-1-A01_bc2093-bc2093.hifi_reads.fastq.gz filtered_out_read_names.txt > Cbombi_reads_Minus1000bp.fastq
```
## HiFiasm to assemble
Use the following slurm script [hifi_slurm.sh](https://github.com/Franck-Dumetz/Crithidia_bombi/blob/main/genome_assembly/hifi_slurm.sh) <br />
From the haplotype condensed fasta file, rename the T2T contigs with numbers in increasing size and the one of 1T with letters in increasing size and sort the fasta file.
```
seqkit sort Cbombi_genome_refined.fasta > CbWHA1_assembly.final.fasta
```
## BUSCO analysis
```
export PATH=/usr/local/packages/metaeuk-6-a5d39d9/bin:$PATH
export PATH=/usr/local/packages/bbtools-39.32:$PATH

/usr/local/packages/busco-5.4.3/bin/busco -m genome -i CbWHA1_assembly.final.fasta --auto-lineage-euk --out /local/projects-t3/SerreDLab-3/fdumetz/Crithidia/hifi/busco -f
/usr/local/packages/busco-5.4.3/scripts/generate_plot.py -wd /local/projects-t3/SerreDLab-3/fdumetz/Crithidia/hifi/busco/busco
```

## Determining coverage
```
minimap2 -ax map-hifi CbWHA1_assembly.final.fasta EDS10_20240819_R84050_PL14293-001_1-1-A01_bc2093-bc2093.hifi_reads.fastq.gz | samtools sort -o Cbombi_reads2hifi.bam

samtools index Cbombi_reads2hifi.bam

samtools depth -aa Cbombi_reads2hifi.bam > Cbombi_depth.txt
```
## GC content calculation per contigs
```
awk '
BEGIN{FS=""; OFS="\t"}
/^>/{
  if(id!=""){
    if(len>0) printf("%s\t%d\t%.6f\n", id, len, 100*gc/len)
  }
  id=substr($0,2)
  sub(/[ \t].*/, "", id)
  len=0
  gc=0
  next
}
{
  line=toupper($0)
  gsub(/[^ACGT]/, "", line)
  len+=length(line)
  tmp=line
  gc+=gsub(/[GC]/, "", tmp)
}
END{
  if(id!="" && len>0) printf("%s\t%d\t%.6f\n", id, len, 100*gc/len)
}
' CbWHA1_assembly.final.fasta > CbWHA1_assembly.final.fasta.GCcontent.tsv
```
## Finding telomeres
```
seqkit locate -i -p TTAGGGTTAGGG CbWHA1_assembly.final.fasta > HiFiCb_telomere.txt
```
## Finding rDNA loci
We use Leishmania donovani ribosomal RNA sequences to locate them in _C .bombi_
```
/usr/local/packages/ncbi-blast+-2.14.0/bin/blastn -query CbWHA1_assembly.final.fasta -subject Ld5.8S.fasta -outfmt 7 -out ./HiFiCbBlastLd5S.txt
```
## Identifying duplicated contigs and misassembly
```
minimap2 -x asm5 -DP CbWHA1_assembly.final.fasta CbWHA1_assembly.final.fasta > self.paf
```
## Finding local CNVs of maximum 1000bp
```
bedtools makewindows -g CbWHA1_assembly.final.fasta.fai -w 1000 > Cb_1000bin.bed
bedtools coverage -a Cb_1000bin.bed -b Cbombi_reads2hifi.bam > Cb_1000bin_cov.txt
```

## Finding the spliced leader sequence
First step: identification of the spliced leader sequence using ONT DRS reads
Extract the first 50 nucleotides of every ONT DRS
```
seqkit subseq -r 1:50 Cbombi_ONT.fasta > Cbombi_first50.fasta
```
Extract 100,000 sequences
```
seqkit split -s 100000 Cbombi_first50.fasta -o Cbombi_split
cat *.fasta > Cbombi_first50_100k.fasta
```
Using MEME to find enriched motifs
```
/usr/local/packages/meme-5.5.5/bin/meme Cbombi_first50_100k.filt.fasta -dna -oc . -mod zoops -nmotifs 10 -minw 20 -maxw 50 -maxsize 100000000000
```
## Ploidy determination using allele frequency
We used a program called [ploidyNGS.py](https://github.com/diriano/ploidyNGS?tab=readme-ov-file) <br />
Whole genome, standard run using only 10% of the reads
```
./ploidyNGS.py --out ploidyNGS/Cbombi_ploidyNGS_guessPloidy --bam hifi/Cbombi_reads2hifi10%.bam

```

Whole genome, all the reads, using ploidy guess
```
./ploidyNGS.py --out ploidyNGS/Cbombi_ploidyNGS_guessPloidy --bam hifi/Cbombi_reads2hifi.bam -u 10 -d 200 -m 0.95 -g

```
Per chromosome (using slurm)
```
N=$(find /local/projects-t3/SerreDLab-3/fdumetz/Crithidia/hifi/ploidy/ploidyNGS_perChr -maxdepth 1 -name "*.bam" | wc -l)

sbatch --array=1-"$N" /home/fdumetz/ploidyNGS/ploidyNGS_perChr.sh
```
