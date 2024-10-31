#!/bin/bash
#SBATCH --job-name=mTR
#SBATCH --mail-user=rossellwood@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=3            # Number of CPU cores per task
#SBATCH --mem=10G
#SBATCH --time=12:00:00
#SBATCH --account=lien.nguyen
#SBATCH --qos=lien.nguyen



# Running on the big fastq file
# /usr/bin/time -v mTR pool-2_2060_PK1sva_bc1009_Male_AD_20R1_56R2_p1.fasta -m 0.8 > pool-2_2060_PK1sva_bc1009_Male_AD_20R1_56R2_p1_mTR_repeat_finder.txt

# Running on the short and long hexamer repeat fastq files: 

#first have to conver them into fasta files:
#short:
# sed -n '1~4s/^@/>/p;2~4p' /blue/lien.nguyen/rossellwood/repeat_expansion_tools/sort_expansion_type/pool-2_2060_PK1sva_bc1009_Male_AD_20R1_56R2_p1_short.fastq > pool-2_2060_PK1sva_bc1009_Male_AD_20R1_56R2_p1_short.fasta
# #long:
# sed -n '1~4s/^@/>/p;2~4p' /blue/lien.nguyen/rossellwood/repeat_expansion_tools/sort_expansion_type/pool-2_2060_PK1sva_bc1009_Male_AD_20R1_56R2_p1_long.fastq > pool-2_2060_PK1sva_bc1009_Male_AD_20R1_56R2_p1_long.fasta


# #Now I will run mTR:

# /usr/bin/time -v mTR pool-2_2060_PK1sva_bc1009_Male_AD_20R1_56R2_p1_short.fasta -m 0.8 > pool-2_2060_PK1sva_bc1009_Male_AD_20R1_56R2_p1_short_mTR_repeat_finder.txt
# /usr/bin/time -v mTR pool-2_2060_PK1sva_bc1009_Male_AD_20R1_56R2_p1_long.fasta -m 0.8 > pool-2_2060_PK1sva_bc1009_Male_AD_20R1_56R2_p1_long_mTR_repeat_finder.txt



#next sample
#short
# sed -n '1~4s/^@/>/p;2~4p' /blue/lien.nguyen/rossellwood/repeat_expansion_tools/sort_expansion_type/pool-1_1054_PK1sva_bc1010_Female_Control_19R1_53R2_p1.fastq > pool-1_1054_PK1sva_bc1010_Female_Control_19R1_53R2_p1.fasta
# /usr/bin/time -v mTR pool-1_1054_PK1sva_bc1010_Female_Control_19R1_53R2_p1.fasta -m 0.8 > pool-1_1054_PK1sva_bc1010_Female_Control_19R1_53R2_p1_mTR_repeat_finder.txt


# #long
# sed -n '1~4s/^@/>/p;2~4p' /blue/lien.nguyen/rossellwood/repeat_expansion_tools/sort_expansion_type/pool-1_1054_PK1sva_bc1010_Female_Control_19R1_53R2_p1_long.fastq > pool-1_1054_PK1sva_bc1010_Female_Control_19R1_53R2_p1_long.fasta
# /usr/bin/time -v mTR pool-1_1054_PK1sva_bc1010_Female_Control_19R1_53R2_p1_long.fasta -m 0.8 > pool-1_1054_PK1sva_bc1010_Female_Control_19R1_53R2_p1_long_mTR_repeat_finder.txt


#pattern 3

# sed -n '1~4s/^@/>/p;2~4p' /blue/lien.nguyen/rossellwood/amplicon_pool_data_Huong/pattern_3/pool-1_1361_PK1sva_bc1011_Female_Control_17R1_43R2_p3.fastq > pool-1_1361_PK1sva_bc1011_Female_Control_17R1_43R2_p3.fasta
# /usr/bin/time -v mTR pool-1_1361_PK1sva_bc1011_Female_Control_17R1_43R2_p3.fasta -m 0.8 > pool-1_1361_PK1sva_bc1011_Female_Control_17R1_43R2_p3_mTR_repeat_finder.txt

sed -n '1~4s/^@/>/p;2~4p' /blue/lien.nguyen/rossellwood/amplicon_pool_data_Huong/pattern_3/pool-1_A15-021_PK1sva_bc1014_Male_AD_23R1_42R2_p3.fastq > pool-1_A15-021_PK1sva_bc1014_Male_AD_23R1_42R2_p3.fasta
/usr/bin/time -v mTR pool-1_A15-021_PK1sva_bc1014_Male_AD_23R1_42R2_p3.fasta -m 0.8 > pool-1_A15-021_PK1sva_bc1014_Male_AD_23R1_42R2_p3_mTR_repeat_finder.txt
