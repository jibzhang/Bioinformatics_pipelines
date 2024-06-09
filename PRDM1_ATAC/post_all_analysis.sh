#!/bin/bash
#SBATCH --job-name=TOBIAS         # Job name
#SBATCH -n 4                          # Number of cores
#SBATCH -N 1-1                        # Min - Max Nodes
#SBATCH -p defq                        # gpu queue
#SBATCH --mem=60G            # Amount of memory in GB
#SBATCH --time=96:00:00               # Time limit hrs:min:sec
#SBATCH --output=/coh_labs/jochan/PRDM1_ATAC/smk_log/TOBIAS_%j.log    # Standard output and error log
snakemake -s PRDM1_ATAC.smk -p /coh_labs/jochan/PRDM1_ATAC/log/TOBIAS_bindetect.bmk -j4