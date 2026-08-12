# Genome comparison




## Spliced-Leader sequence in other Trypanosomatids (_Crithidia fasciculata_ and _Lotmaria passim_)
Using blast to find/identify the SL
```
makeblastdb -in TriTrypDB-68_CfasciculataCfCl_Genome.fasta -dbtype nucl -out Cf_db
makeblastdb -in Lpassim_GCA_037349495.1_ASM3734949v1_genomic.fna -dbtype nucl -out Lpassim_db

/usr/local/packages/ncbi-blast+-2.14.0/bin/blastn -query CbSL.fasta -db Cf_db -out SL2Cf.blastout -outfmt 6
/usr/local/packages/ncbi-blast+-2.14.0/bin/blastn -query CbSL.fasta -db Lpassim_db -out SL2Lpassim.blastout -outfmt 6
```

