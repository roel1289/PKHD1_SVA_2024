#!/bin/bash
#SBATCH --job-name=G4/G3_ratio
#SBATCH --mail-user=rossellwood@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1            # Number of CPU cores per task
#SBATCH --mem=2G
#SBATCH --time=12:00:00
#SBATCH --account=lien.nguyen
#SBATCH --qos=lien.nguyen-b


###########################################################
### Looking into LHex_Hexv5 alleles:
###########################################################
# making an array of the filenames that I am interested in:
file_list=("pool-5_ADRC1862_PK1sva_bc1009_Male_AD_29R1_42R2_p5.fastq_filtered.fastq"
 "pool-5_2153_PK1sva_bc1014_Male_AD_31R1_44R2_p4.fastq_filtered.fastq" 
 "pool-5_ADRC2014_PK1sva_bc1020_Male_AD_32R1_48R2_p5.fastq_filtered.fastq"
 "pool-4_A15-019_PK1sva_bc1020_Female_AD_27R1_29R2_p5.fastq_filtered.fastq"
 "pool-5_ADRC2014_PK1sva_bc1020_Male_AD_32R1_48R2_p5.fastq_filtered.fastq"
 "pool-5_2153_PK1sva_bc1014_Male_AD_31R1_44R2_p4.fastq_filtered.fastq"
 "pool-5_ADRC1862_PK1sva_bc1009_Male_AD_29R1_42R2_p5.fastq_filtered.fastq")
###########################################################

###########################################################
### Looking into LHex_hexv0 alleles:
###########################################################
# file_list=("")
###########################################################

# Loop through each file in the array
for file in "${file_list[@]}"
do
    echo "$file"
    /blue/lien.nguyen/rossellwood/amplicon_pool_data_Huong/scripts/SVA_hexv5_search.py \
    -f ../../repeat_expansion_tools/sort_expansion_type/filtered_full_FASTQ/$file \
    -o test | sort | uniq -c
done

