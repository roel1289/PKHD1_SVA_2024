#!/bin/bash
#SBATCH --job-name=mTR
#SBATCH --mail-user=rossellwood@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1            # Number of CPU cores per task
#SBATCH --mem=2G
#SBATCH --time=12:00:00
#SBATCH --account=lien.nguyen
#SBATCH --qos=lien.nguyen-b


# By changing the file paths, the program, mTR, can be ran on different FASTQ files. 
# The script transforms them into FASTA format, then applies the mTR program to the files. 

sed -n '1~4s/^@/>/p;2~4p' /blue/lien.nguyen/rossellwood/amplicon_pool_data_Huong/demux_out/demux_fastq/demux_fastq_perSample/pool-1_A20-002_PK1sva_bc1013_Male_AD_47R1_48R2_p4.fastq > pool-1_A20-002_PK1sva_bc1013_Male_AD_47R1_48R2_p4.fasta
/usr/bin/time -v mTR pool-1_A20-002_PK1sva_bc1013_Male_AD_47R1_48R2_p4.fasta -m 0.8 > pool-1_A20-002_PK1sva_bc1013_Male_AD_47R1_48R2_p4_mTR_repeat_finder.txt
