#!/bin/bash
#SBATCH --job-name=Cbombi_asm                     # Job name
#SBATCH --output=/local/projects-t3/SerreDLab-3/fdumetz/Crithidia/hifi/Cbhifiasm.out   # Standard output log
#SBATCH --error=/local/projects-t3/SerreDLab-3/fdumetz/Crithidia/hifi/Cbhifiasm.err    # Standard error log
#SBATCH --mail-type=BEGIN,END --mail-user=fdumetz@som.umaryland.edu
#SBATCH --cpus-per-task=32				# Number of CPUs per task
#SBATCH --mem=700G                                        # Memory per node


/usr/local/packages/hifiasm-0.24.0/bin/hifiasm -o Cb_hifi.asm -t32 --telo-m TTAGGG --dual-scaf /local/projects-t3/SerreDLab-3/fdumetz/Crithidia/Cbombi_hifi_filtered.fastq.gz