# RNA modifications detection 

---

## Inosine and m6A
```
pod_directory=/path/to/Pod5
fasta=/path/to/fasta

/usr/local/packages/dorado-0.8.1/bin/dorado basecaller /usr/local/packages/dorado-0.8.1/models/rna004_130bps_sup@v5.1.0 $pod_directory --reference $fasta --modified-bases-models /usr/local/packages/dorado-0.8.1/models/rna004_130bps_sup@v5.1.0_inosine_m6A@v1 -x cuda:1,2 --emit-sam > Cbombi_inosine_m6A.sam
```
## m6A DRACH
```
pod_directory=/path/to/Pod5
fasta=/path/to/fasta

/usr/local/packages/dorado-0.8.1/bin/dorado basecaller /usr/local/packages/dorado-0.8.1/models/rna004_130bps_sup@v5.1.0 /usr/local/packages/dorado-0.8.1/models/rna004_130bps_sup@v5.1.0 $pod_directory --reference $fasta --modified-bases-models /usr/local/packages/dorado-0.8.1/models/rna004_130bps_sup@v5.1.0_m6A_DRACH@v1 -x cuda:7 --emit-sam > Cbombi_DRACH_m6A.sam
```
find DRACH motif in the genome
```
grep -o -i '[AGT][AG]AC[ACT]' $fasta | wc -l
```
## m5C
```
pod_directory=/path/to/Pod5
fasta=/path/to/fasta

/usr/local/packages/dorado-0.8.1/bin/dorado basecaller /usr/local/packages/dorado-0.8.1/models/rna004_130bps_sup@v5.1.0 $pod_directory --reference $fasta --modified-bases-models /usr/local/packages/dorado-0.8.1/models/rna004_130bps_sup@v5.1.0_m5C@v1 -x cuda:3,4 --emit-sam > Cbombi_mC5.sam
```
## pseudoU
```
pod_directory=/path/to/Pod5
fasta=/path/to/fasta

/usr/local/packages/dorado-0.8.1/bin/dorado basecaller /usr/local/packages/dorado-0.8.1/models/rna004_130bps_sup@v5.1.0 $pod_directory --reference $fasta --modified-bases-models /usr/local/packages/dorado-0.8.1/models/rna004_130bps_sup@v5.1.0_pseU@v1 -x cuda:5,6 --emit-sam > Cbombi_pseudoU.sam
```

## Parsing of the data
```
samtools view -bhF 2308 aligned.sam | samtools sort -o bam_sorted.bam
samtools index bam_sorted.bam
```

## Per-read comparison
Calling mods on modkit
| Modification Name     | modkit Code | Example Threshold Flag     | Notes                                 |
|------------------------|-------------|-----------------------------|----------------------------------------|
| N6-methyladenosine     | `m`         | `--mod-threshold m:0.9`     | Commonly used m6A motif (DRACH)        |
| 5-methylcytosine       | `c`         | `--mod-threshold c:0.9`     | Cytosine methylation (m5C)             |
| Pseudouridine (Ψ)      | `h`         | `--mod-threshold h:0.9`     | Calling ofr PseudoUridine              |
| Inosine                | `i`         | `--mod-threshold i:0.9`     | A-to-I RNA editing                     |

```
/usr/local/packages/modkit-0.4.4/modkit pileup --ref /local/projects-t3/SerreDLab-3/fdumetz/Crithidia/hifi/Cbombi_genome_refined.fasta --mod-threshold h:0.9 Cbombi_pseudoU_hac.sorted.bam Cbombi_pseudoU_hac.bed

/usr/local/packages/modkit-0.4.4/modkit summary Cbombi_pseudoU_hac.sorted.bam
```
Use pass threshold from the summary output has start to establish a threshold for later <br />

## Preparing the bed file
at least 10 reads at the position and pecentage modified different than zero
```
awk '$10 >= 10 && $11 != 0' Cbombi_inosine_m6A.bed | head
```
