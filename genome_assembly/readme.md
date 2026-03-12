# Crithidia bombi genome assembly

BioProject PRJNA1231061

Table of content: <br />


Software requirements: <br />
- bedtools-2.27.1 <br />
- busco-5.4.3 <br />
- filtlong-0.2.0 <br />
- hifiasm-0.24.0 <br />
- minimap2.1 <br />
- mummer-3.23 <br />
- quast-5.2.0 <br />
- samtools-1.20 <br />
- seqkit-0.7.2 <br />
- seqtk-1.0-r63 <br />
- trnascan-se-2.0.3 <br />

## Preparing the data by filtering for reads longer than 1000 bp
```
/usr/local/packages/filtlong-0.2.0/bin/filtlong --min_length 1000 --keep_percent 90 /local/projects-t3/EDS10/Cbombi_WHA1_20240805/PACBIO_DATA/EDS10_20240819_R84050_PL14293-001_1-1-A01_bc2093-bc2093.hifi_reads.fastq.gz | gzip > Cbombi_hifi_filtered.fastq.gz

zgrep '^@' /local/projects-t4/aberdeen2ro/SerreDLab-4/raw_reads/2024-08-22_Pacbio/Cbombi_WHA1_20240805/PACBIO_DATA/EDS10_20240819_R84050_PL14293-001_1-1-A01_bc2093-bc2093.hifi_reads.fastq.gz > original_read_names.txt

zgrep '^@' Cbombi_hifi_filtered.fastq.gz > filtered_read_names.txt
```
Find the Read Names that Were Filtered Out
```
grep -v -f filtered_read_names.txt original_read_names.txt > filtered_out_read_names.txt
```
## Extract reads and create a new fastq
```
seqtk subseq /local/projects-t4/aberdeen2ro/SerreDLab-4/raw_reads/2024-08-22_Pacbio/Cbombi_WHA1_20240805/PACBIO_DATA/EDS10_20240819_R84050_PL14293-001_1-1-A01_bc2093-bc2093.hifi_reads.fastq.gz filtered_out_read_names.txt > Cbombi_reads_Minus1000bp.fastq
```
## HiFiasm to assemble
Use the following slurm script [hifi_slurm.sh](https://github.com/Franck-Dumetz/Crithidia_bombi/blob/main/genome_assembly/hifi_slurm.sh) <br />
From the haplotype condensed fasta file, rename the T2T contigs with numbers in increasing size and the one of 1T with letters in increasing size and sort the fasta file.
```
seqkit sort Cbombi_genome_refined.fasta > Cbombi_genome_refined.sorted.fasta
```
## BUSCO analysis
```
/usr/local/packages/busco-5.4.3/bin/busco -m genome -i /local/projects-t3/SerreDLab-3/fdumetz/Crithidia/hifi/Cb_hifi.asm.bp.p_ctg.fasta --auto-lineage-euk --out /local/projects-t3/SerreDLab-3/fdumetz/Crithidia/hifi/busco --long -f
/usr/local/packages/busco-5.4.3/scripts/generate_plot.py -wd /local/projects-t3/SerreDLab-3/fdumetz/Crithidia/hifi/busco/busco
```

## Determining coverage
```
minimap2 -ax map-hifi /local/projects-t3/SerreDLab-3/fdumetz/Crithidia/hifi/Cb_hifi.asm.bp.p_ctg.fasta /local/projects-t4/aberdeen2ro/SerreDLab-4/raw_reads/2024-08-22_Pacbio/Cbombi_WHA1_20240805/PACBIO_DATA/EDS10_20240819_R84050_PL14293-001_1-1-A01_bc2093-bc2093.hifi_reads.fastq.gz | samtools sort -o Cbombi_reads2hifi.bam

samtools index Cbombi_reads2hifi.bam

samtools depth Cbombi_reads2hifi.bam > Cbombi_depth.txt
awk '{print $3}' Cbombi_depth.txt
```
## Indentifying duplicated contings and missassembly
```
minimap2 -x asm5 -DP Cb_hifi.asm.bp.p_ctg.fasta Cb_hifi.asm.bp.p_ctg.fasta > self.paf
```
## GC content calculation per contigs
```
awk '
>   BEGIN{FS=""; OFS="\t"}
>   /^>/{if(id!=""){
>          if(len>0) printf("%s\t%d\t%.6f\n", id, len, 100*gc/len);
>        }
>        id=substr($0,2); sub(/[ \t].*/,"",id);
>        len=0; gc=0; next}
>   {line=toupper($0);
>    gsub(/[^ACGT]/,"",line);
>    len+=length(line);
>    tmp=line; gc+=gsub(/[GC]/,"",tmp)}
>   END{if(id!="" && len>0) printf("%s\t%d\t%.6f\n", id, len, 100*gc/len)}
> ' Cbombi_genome_refined.fasta > gc_per_contig.tsv
```
## Finding telomeres
```
seqkit locate -i -p TTAGGGTTAGGG Cb_hifi.asm.bp.p_ctg.fasta > HiFiCb_telomere.txt
```

## Finding rDNA loci
We use Leishmania donovani ribosomal RNA sequences to locate them in _C .bombi_
```
/usr/local/packages/ncbi-blast+-2.14.0/bin/blastn -query Cbombi_genome_refined_renamed.fasta -subject Ld5.8S.fasta -outfmt 7 -out ./HiFiCbBlastLd5S.txt
```
## Finding local CNVs of maximum 1000bp
```
bedtools makewindows -g /local/projects-t3/SerreDLab-3/fdumetz/Crithidia/hifi/Cb_hifi.asm.bp.hap1.p_ctg.fasta.fai -w 1000 > Cb_1000bin.bed
bedtools coverage -a Cb_1000bin.bed -b /local/projects-t3/SerreDLab-3/fdumetz/Crithidia/hifi/Cbombi_reads2hifi.bam > Cb_1000bin_cov.txt
```

## Finding the spliced leader sequence
First step: identification of the spliced leader sequence using ONT DRS reads
Extract the 50 first nucleotides of every ONT DRS
```
seqkit subseq -r 1:50 Cbombi_ONT.fasta > Cbombi_first50.fasta
```
Extract 100,000 sequences
```
seqkit split -s 100000 Cbombi_first50.fasta -o Cbombi_split
cat *.fasta > Cbombi_first50_100k.fasta
```
Using meme to find enriched motifs
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
