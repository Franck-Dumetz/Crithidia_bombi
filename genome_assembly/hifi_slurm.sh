#!/bin/bash
#SBATCH --job-name=Cbombi_asm                     # Job name
#SBATCH --output=/path/Cbhifiasm.out   # Standard output log
#SBATCH --error=/path/Cbhifiasm.err    # Standard error log
#SBATCH --mail-type=BEGIN,END --mail-user=email
#SBATCH --cpus-per-task=32				# Number of CPUs per task
#SBATCH --mem=700G                                        # Memory per node


/usr/local/packages/hifiasm-0.24.0/bin/hifiasm -o Cb_hifi.asm -t32 --telo-m TTAGGG --dual-scaf /path/Cbombi_hifi_filtered.fastq.gz
