# Using Direct RNA sequencing from Oxford Nanopore Technology to annotate the _C. bombi_ genome

BioSample

Table of content: <br />


Software requirements: <br />
- minimap-2.1
- samtools-1.20
- stringtie-2.2.1
- 

## Aligning reads to the reference genome
Minimap was used as described in [Minimap_ONT.sh](https://github.com/Franck-Dumetz/Crithidia_bombi/blob/main/Annotation/Minimap_ONT.sh) <br />

Parsing the reads
```
samtools view -bhF 2308 Cbombi_ONT.sam | samtools sort -o Cbombi_ONT.bam
samtools index Cbombi_ONT.bam
```
## Spliced-leader identification 

## Transcript evidence using DRS 
