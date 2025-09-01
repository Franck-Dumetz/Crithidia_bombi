# Using Direct RNA sequencing from Oxford Nanopore Technology to annotate the _C. bombi_ genome

BioSample

Table of content: <br />


Software requirements: <br />
- dorado-0.8.1
- minimap-2.1
- samtools-1.20
- stringtie-2.2.1
- transdecoder-5.7.1

## Basecalling using Dorado
```
dorado-0.8.1/bin/dorado basecaller --device cuda:1 --emit-fastq --min-qscore 7 /usr/local/packages/dorado-0.8.1/models/rna004_130bps_sup@v5.1.0 /local/projects-t2/RDMIN/SEQUENCE/20250513-MN23690_Cbombi_swimming_FranckDumetz-Blyss/20250513-MN23690/Cbombi_swimming/20250513_1621_MN23690_FBC11489_6cc52add/pod5_skip > Cbombi_swimming_ONT.fastq
```

## Aligning reads to the reference genome
Minimap was used as described in [Minimap_ONT.sh](https://github.com/Franck-Dumetz/Crithidia_bombi/blob/main/Annotation/Minimap_ONT.sh) <br />

Parsing the reads
```
samtools view -bhF 2308 Cbombi_ONT.sam | samtools sort -o Cbombi_ONT.bam
samtools index Cbombi_ONT.bam
```
## Filtering for reads with a Spliced-leader and a polyA tail of at least 30 As 
Calling for polyA tail
```
dorado-0.8.1/bin/dorado basecaller --device cuda:all --emit-sam --estimate-poly-a --min-qscore 7 /usr/local/packages/dorado-0.8.1/models/rna004_130bps_sup@v5.1.0 /local/projects-t2/RDMIN/SEQUENCE/20250513-MN23690_Cbombi_swimming_FranckDumetz-Blyss/20250513-MN23690/Cbombi_swimming/20250513_1621_MN23690_FBC11489_6cc52add/pod5_skip > Cbombi_swimming_polyA_ONT.sam
```

## Semi-automated transcript evidence generation 
```
stringtie-2.2.3/stringtie /local/projects-t3/SerreDLab-3/fdumetz/Crithidia/ONT/Cbombi_ONT.bam --fr -f 0.5 -c 3 -l CbWHA -p 8 -L -R -o Cbombi_transcript.gtf
```
Stringtie is "isoform aware" but doesn't really handle polycistron well. <br />
First filter the transcript by TPM and coverage. To ensure a medium filtering that will remove polycistron and low abundance transcripts, set TPM to 2.5 and cov to 8 using [filter_stringtie_gtf.py](https://github.com/Franck-Dumetz/Crithidia_bombi/blob/main/Annotation/filter_stringtie_gtf.py). <br />
A manual step to remove the wrong transcript evidence annotation was added using []() <br />
To perform this step, open IGV with the bam file used to create the stringtie gft and the transdecoder gff3 output. Then identify every transcript evidence that is wrong based on ORF presence and soft clipping on reads on the right AND the left of the transcript indicating of the presence of the spliced-leader and of the polyA tail. Remove every transcript evidence that is overlapping, this is indicative of a polycistron. Record the transcript id in a txt file and use []() to remove them from the gft file. <br />
Another step is to add missing annatation was aded using []() and to correct start/end using []()<br /> 
## Identification of the longest ORF
```
transdecoder-5.7.1/util/gtf_genome_to_cdna_fasta.pl /local/projects-t3/SerreDLab-3/fdumetz/Crithidia/stringtie/Cbombi_transcript.gtf /local/projects-t3/SerreDLab-3/fdumetz/Crithidia/hifi/Cbombi_genome_refined.fasta > Transdecoder_transcripts_2.5-8_Cb.fasta
```
```
transdecoder-5.7.1/TransDecoder.LongOrfs -t Transdecoder_transcripts_2.5-8_Cb.fasta
```
```
transdecoder-5.7.1/TransDecoder.Predict -t Transdecoder_transcripts_2.5-8_Cb.fasta --single_best_only
```
```
/usr/local/packages/transdecoder-5.7.1/util/gtf_to_alignment_gff3.pl /local/projects-t3/SerreDLab-3/fdumetz/Crithidia/stringtie/Cbombi_transcript.gtf > Cb_2.8-5.gff3
```
```
transdecoder-5.7.1/util/cdna_alignment_orf_to_genome_orf.pl Transdecoder_transcripts_2.5-8_Cb.fasta.transdecoder.gff3 Cb_2.8-5.gff3 Transdecoder_transcripts_2.5-8_Cb.fasta > Cb_2.5-8_annotation.transdecoder.genome.gff3
```
## 
