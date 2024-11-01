#!/bin/bash
#SBATCH --job-name=allele_call
#SBATCH --mail-user=rossellwood@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2            # Number of CPU cores per task
#SBATCH --mem=3G
#SBATCH --time=12:00:00
#SBATCH --account=lien.nguyen
#SBATCH --qos=lien.nguyen-b

# This script runs my allele_call.py script. 

/usr/bin/time -v ./allele_call.py -d /blue/lien.nguyen/rossellwood/repeat_expansion_tools/sort_expansion_type/filtered_full_FASTQ -o alleleCall_filtered18.txt