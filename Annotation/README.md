# Using Direct RNA sequencing from Oxford Nanopore Technology to annotate the _C. bombi_ genome

BioSample

Table of content: <br />


Software requirements: <br />
- dorado-0.8.1
- EMBOSS
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
## Filtering for reads with a Spliced-leader and a polyA tail of at least 20 As 
Filtering for the reads with a spliced-leader sequence in 5' <br />
Use the [SL_finding.sh](https://github.com/Franck-Dumetz/Crithidia_bombi/blob/main/Annotation/SL_finding.sh) script to identify the reads with a SL. Basically, that script executes the EMBOSS subpackage fuzznuc to find motifs with a certain number of mismatches (to accommodate ONT miscalling). Then it filters for only the motifs that are in the 30 first nucleotide of the read. <br />
Use a U-to-T converted fasta file to work with fuzznuc. <br />
```
/SL_finding.sh -s ATAAGTATCAGTTTCTGTACTTTATTG -F /path/to/directory -m 4 -l 10 -o SL_calls_4 ONT_reads.fasta
```
Extracting read names from 
```
awk '{print $3}' Cbombi_swimming_ONT_SLlen19.tsv > SL19_read_names.txt
```
Extracting reads with at least 19 nucleotides from the SL sequence
```
SAMTOOLS=/usr/local/packages/samtools-1.20/bin/samtools

BAM="Cbombi_ONT.bam"

NAMES="SL10_read_names.txt"

# output names

OUT_PREFIX="Cbombi_ONT.SL10"

OUT_BAM="${OUT_PREFIX}.sorted.bam"

# filter -> sort -> index

$SAMTOOLS view -b -N "$NAMES" "$BAM" | $SAMTOOLS sort -o "$OUT_BAM"

$SAMTOOLS index "$OUT_BAM"
```

## Semi-automated transcript evidence generation 
```
/usr/local/packages/stringtie-2.2.3/stringtie Cbombi_ONT.SL10.sorted.bam -L -p 8 -f 0.05 -l CbWHA1 -c 1 -g 10 -m 50 -o Cbombi_SL10_fullOption.gtf
```
Stringtie is "isoform aware" but doesn't really handle polycistron well. <br />
First, filter the transcript by TPM and coverage. To ensure a medium filtering that will remove polycistron and low-abundance transcripts, set TPM to 2.5 and cov to 8 using [filter_stringtie_gtf.py](https://github.com/Franck-Dumetz/Crithidia_bombi/blob/main/Annotation/filter_stringtie_gtf.py). <br />
A manual step to remove the wrong transcript evidence annotation was added using []() <br />
To perform this step, open IGV with the BAM file used to create the StringTie GFF and the TransDecoder GFF3 output. Then identify every transcript evidence that is wrong based on ORF presence and soft clipping on reads on the right AND the left of the transcript, indicating the presence of the spliced-leader and of the polyA tail. Remove any overlapping transcript evidence; this indicates a polycistron. Record the transcript ID in a to_remove.txt file and use the following line to remove them from the gtf. <br />
```
awk '
  NR==FNR { gsub(/\r/,"",$1); if ($1!="") bad[$1]=1; next }

  # parse transcript_id and gene_id from 9th column
  {
    tid=""; gid=""
    if (match($0, /transcript_id "[^"]+"/)) { tid=substr($0, RSTART+14, RLENGTH-15) }
    if (match($0, /gene_id "[^"]+"/))       { gid=substr($0, RSTART+9,  RLENGTH-10) }

    if ((tid != "" && bad[tid]) || (gid != "" && bad[gid])) next
    print
  }
' to-remove_stg_020226.txt Cbombi_transcript_TPM2.5_cov8.CbWHA1.renamed.gtf > filtered.gtf
```
Now comes the manual curation. Go transcript by transcript over the entire genome, creating 3 different files, and use the following script in the following order:
  - remove inaccurate annotations using [remove_transcripts_csv.py](https://github.com/Franck-Dumetz/Crithidia_bombi/blob/main/Annotation/remove_transcripts_csv.py)
```
python remove_transcripts_csv.py --gtf your_input.gtf --csv to_remove.csv --out output.gtf
```
  - add missing annotations, or isoforms, using [add_transcripts_csv.py](https://github.com/Franck-Dumetz/Crithidia_bombi/blob/main/Annotation/add_transcripts_csv.py)
```
python add_transcripts_csv.py --gtf your_input.gtf --csv new_entries.csv --out output.gtf
```
  - Finally, correct annotations that are not correct using [correct_transcripts_csv.py](https://github.com/Franck-Dumetz/Crithidia_bombi/blob/main/Annotation/correct_transcripts_csv.py)
```
python correct_transcripts_csv.py --gtf your_input.gtf --csv to_correct.csv --out output.gtf
```
Now, we can rename the transcripts in the following format CbWHA1.chr#.000000.isoform# using [rename_loci_by_overlap.py](https://github.com/Franck-Dumetz/Crithidia_bombi/blob/main/Annotation/rename_loci_by_overlap.py). It also outputs a mapping key between the old names and the new names
```
python rename_loci_by_overlap.py --gtf input.gtf --out Cbombi_stringtie_final.gtf --map correspondence_key.csv
```

## Identification of the longest ORF
```
transdecoder-5.7.1/util/gtf_genome_to_cdna_fasta.pl Cbombi_stringtie_final.gtf Cbombi_assembly_final.fasta > Transdecoder_transcripts_2.5-8_Cb.fasta
```
```
transdecoder-5.7.1/TransDecoder.LongOrfs -t Transdecoder_transcripts_2.5-8_Cb.fasta
```
```
transdecoder-5.7.1/TransDecoder.Predict -t Transdecoder_transcripts_2.5-8_Cb.fasta --single_best_only
```
```
transdecoder-5.7.1/util/gtf_to_alignment_gff3.pl Cbombi_stringtie_final.gtf > Cb_2.8-5.gff3
```
```
transdecoder-5.7.1/util/cdna_alignment_orf_to_genome_orf.pl Transdecoder_transcripts_2.5-8_Cb.fasta.transdecoder.gff3 Cb_2.8-5.gff3 Transdecoder_transcripts_2.5-8_Cb.fasta > Cb_2.5-8_annotation.transdecoder.genome.gff3
```
## 
