#!/bin/bash
#SBATCH --job-name=Cbombi_asm                     # Job name
#SBATCH --output=Cbhifiasm.out   # Standard output log
#SBATCH --error=Cbhifiasm.err    # Standard error log
#SBATCH --mail-type=BEGIN,END --mail-user=email
#SBATCH --cpus-per-task=32				# Number of CPUs per task
#SBATCH --mem=700G                                        # Memory per node


hifiasm -o Cb_hifi.asm -t32 --telo-m TTAGGG --dual-scaf Cbombi_hifi_filtered.fastq.gz
