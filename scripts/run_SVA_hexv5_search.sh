#!/bin/bash
#SBATCH --job-name=hex_primer-dist
#SBATCH --mail-user=rossellwood@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1            # Number of CPU cores per task
#SBATCH --mem=2G
#SBATCH --time=12:00:00
#SBATCH --account=lien.nguyen
#SBATCH --qos=lien.nguyen-b


# By uncommenting certain parts of this script, different tasks can be done. 
# This script can analyze LHex_Hexv5 and Lhex_hexv0 files, 
# as well as run python script that analyzes distance between hexamer and primer


###########################################################
### Looking into LHex_Hexv5 alleles:
###########################################################
# making an array of the filenames that I am interested in:
# file_list=("pool-5_ADRC1862_PK1sva_bc1009_Male_AD_29R1_42R2_p5.fastq_filtered.fastq"
#  "pool-5_2153_PK1sva_bc1014_Male_AD_31R1_44R2_p4.fastq_filtered.fastq" 
#  "pool-5_ADRC2014_PK1sva_bc1020_Male_AD_32R1_48R2_p5.fastq_filtered.fastq"
#  "pool-4_A15-019_PK1sva_bc1020_Female_AD_27R1_29R2_p5.fastq_filtered.fastq"
#  "pool-5_ADRC2014_PK1sva_bc1020_Male_AD_32R1_48R2_p5.fastq_filtered.fastq"
#  "pool-5_2153_PK1sva_bc1014_Male_AD_31R1_44R2_p4.fastq_filtered.fastq"
#  "pool-5_ADRC1862_PK1sva_bc1009_Male_AD_29R1_42R2_p5.fastq_filtered.fastq")
# ###########################################################

# ###########################################################
# ### Looking into LHex_hexv0 alleles:
# ###########################################################
# # file_list=("")
# ###########################################################

# # Loop through each file in the array
# for file in "${file_list[@]}"
# do
#     echo "$file"
#     /blue/lien.nguyen/rossellwood/amplicon_pool_data_Huong/scripts/SVA_hexv5_search.py \
#     -f ../../repeat_expansion_tools/sort_expansion_type/filtered_full_FASTQ/$file \
#     -o test | sort | uniq -c
# done


### RUN THIS PART TO GENERATE FILE WITHE EVERY LENGTH BETWEEN HEXAMER REPEAT AND PRIMER (TO COMPARE HEXV0 AND HEXV1)  
fastq_dir="/blue/lien.nguyen/rossellwood/repeat_expansion_tools/sort_expansion_type/filtered_full_FASTQ"
cd $fastq_dir
files=$(ls -1 *.fastq)

for file in $files
do 
    echo $file
    /blue/lien.nguyen/rossellwood/amplicon_pool_data_Huong/scripts/SVA_hexv5_search.py -f $file >> hex_primer_dist.txt
done
