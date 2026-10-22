#!/bin/bash

#SBATCH --job-name=Cb_map                      # Job name
#SBATCH --output=path/ONT_DRS.out   # Standard output log
#SBATCH --error=path/ONT_DRS.err    # Standard error log
#SBATCH --mail-type=BEGIN,END --mail-user=fdumetz@som.umaryland.edu
#SBATCH --cpus-per-task=32				# Number of CPUs per task
#SBATCH --mem=30G                              # Memory per node
#SBATCH --account=serre-lab

minimap2 -ax map-ont -k14 -t 32 path/Cbombi_genome_refined.fasta path/Cbombi_swimming_ONT.fastq.gz | samtools sort -@ 8 -bhF 2308 -o path/Cbombi_ONT.sorted.bam 

samtools index path/Cbombi_ONT.sorted.bam
