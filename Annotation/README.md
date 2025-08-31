# Using Direct RNA sequencing from Oxford Nanopore Technology to annotate the _C. bombi_ genome

BioSample

Table of content: <br />


Software requirements: <br />
- minimap-2.1
- samtools-1.20
- stringtie-2.2.1
- transdecoder-5.7.1

## Aligning reads to the reference genome
Minimap was used as described in [Minimap_ONT.sh](https://github.com/Franck-Dumetz/Crithidia_bombi/blob/main/Annotation/Minimap_ONT.sh) <br />

Parsing the reads
```
samtools view -bhF 2308 Cbombi_ONT.sam | samtools sort -o Cbombi_ONT.bam
samtools index Cbombi_ONT.bam
```
## Spliced-leader identification 

## Transcript evidence using DRS 
```
stringtie-2.2.3/stringtie /local/projects-t3/SerreDLab-3/fdumetz/Crithidia/ONT/Cbombi_ONT.bam --fr -f 0.5 -c 3 -l Ld1S -p 8 -L -R -o Cbombi_transcript.gtf
```
Stringtie is "isoform aware" but doesn't really handle polycistron well. <br />
A manual step to remove the wrong annotation was added using []() and another step to add missing annatation was aded using []() <br /> 

## Identification of the longest ORF
```
transdecoder-5.7.1/util/gtf_genome_to_cdna_fasta.pl /local/projects-t3/SerreDLab-3/fdumetz/Crithidia/stringtie/Cbombi_transcript.gtf /local/projects-t3/SerreDLab-3/fdumetz/Crithidia/hifi/Cbombi_genome_refined.fasta > Transdecoder_transcripts_Cb.fasta
```
```
transdecoder-5.7.1/TransDecoder.LongOrfs -t Transdecoder_transcripts_Cb.fasta
```
```
transdecoder-5.7.1/TransDecoder.Predict -t Transdecoder_transcripts_Cb.fasta --single_best_only
```
```
transdecoder-5.7.1/util/cdna_alignment_orf_to_genome_orf.pl Transdecoder_transcripts_Cb.fasta.transdecoder.gff3 Cb_trasndecoder1st.gff3 Transdecoder_transcripts_Cb.fasta > Cb_full_annotation.transdecoder.genome.gff3
```
## 
