# Genome comparison




## Spliced-Leader sequence in other Trypanosomatids (_Crithidia fasciculata_ and _Lotmaria passim_)
Using blast to find/identify the SL
```
makeblastdb -in TriTrypDB-68_CfasciculataCfCl_Genome.fasta -dbtype nucl -out Cf_db
makeblastdb -in Lpassim_GCA_037349495.1_ASM3734949v1_genomic.fna -dbtype nucl -out Lpassim_db

/usr/local/packages/ncbi-blast+-2.14.0/bin/blastn -query query.fa -db Cf_db -out SL2Cf.blastout -outfmt 6
/usr/local/packages/ncbi-blast+-2.14.0/bin/blastn -query query.fa -db Lpassim_db -out SL2Lpassim.blastout -outfmt 6
```

## Somy estimation using ploidyNGS and CNVkit
We split the reads from the BAM file into single-chromosome BAM files < br>
Then we used [ploidyNGS_perChr.slurm.sh]() to determine the allele frequency per chromosome <br >

CNVkit usage:
```
cnvkit.py batch WHA1.bam \
    -n \
    -f CbWHA1_assembly.fasta \
    -m wgs \
    -p 8 \
    -d cnvkit_out
```
