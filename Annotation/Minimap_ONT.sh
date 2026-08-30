#!/bin/bash

#SBATCH --job-name=Cb_map                      # Job name
#SBATCH --output=/local/projects-t3/SerreDLab-3/fdumetz/Crithidia/ONT/ONT_DRS.out   # Standard output log
#SBATCH --error=/local/projects-t3/SerreDLab-3/fdumetz/Crithidia/ONT/ONT_DRS.err    # Standard error log
#SBATCH --mail-type=BEGIN,END --mail-user=fdumetz@som.umaryland.edu
#SBATCH --cpus-per-task=32				# Number of CPUs per task
#SBATCH --mem=30G                              # Memory per node
#SBATCH --account=serre-lab

minimap2 -ax map-ont -k14 -t 32 /local/projects-t3/SerreDLab-3/fdumetz/Crithidia/hifi/Cbombi_genome_refined.fasta /local/projects-t3/SerreDLab-3/fdumetz/Crithidia/ONT/Cbombi_swimming_ONT.fastq.gz | samtools sort -@ 8 -o /local/projects-t3/SerreDLab-3/fdumetz/Crithidia/ONT/Cbombi_ONT.sorted.bam 

samtools index /local/projects-t3/SerreDLab-3/fdumetz/Crithidia/ONT/Cbombi_ONT.sorted.bam
