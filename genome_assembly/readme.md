# Crithidia bombi genome assembly

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
use the following slurm script [hifi_slurm.sh](https://github.com/Franck-Dumetz/Crithidia_bombi/blob/main/genome_assembly/hifi_slurm.sh) <br />

## Determining coverage
```
minimap2 -ax map-hifi /local/projects-t3/SerreDLab-3/fdumetz/Crithidia/hifi/Cb_hifi.asm.bp.p_ctg.fasta /local/projects-t4/aberdeen2ro/SerreDLab-4/raw_reads/2024-08-22_Pacbio/Cbombi_WHA1_20240805/PACBIO_DATA/EDS10_20240819_R84050_PL14293-001_1-1-A01_bc2093-bc2093.hifi_reads.fastq.gz | samtools sort -o Cbombi_reads2hifi.bam

samtools index Cbombi_reads2hifi.bam

samtools depth Cbombi_reads2hifi.bam > Cbombi_depth.txt
awk '{print $3}' Cbombi_depth.txt
```

## Finding telomeres
```
seqkit locate -i -p TTAGGGTTAGGG Cb_hifi.asm.bp.p_ctg.fasta > HiFiCb_telomere.txt
```

## Finding rDNA loci
We use Leishmania donovani ribosomal RNA sequences to locate them in C .bombi
```
/usr/local/packages/ncbi-blast+-2.14.0/bin/blastn -query /local/projects-t3/SerreDLab-3/fdumetz/Crithidia/hifi/Cb_hifi.asm.bp.p_ctg.fasta -subject /local/projects-t3/SerreDLab-3/fdumetz/Crithidia/canu/rRNA/Ld5SrRNA.fasta -outfmt 7 -out ./HiFiCbBlastLd5S.txt
```
