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
samtools view -bhF 2308 aligned.sam > aligned.bam
samtools sort aligned.bam -o aligned.sorted.bam
samtools index aligned.sorted.bam
```
