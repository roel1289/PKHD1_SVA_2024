#!/bin/bash
#SBATCH --job-name=rename_fastq
#SBATCH --mail-user=rossellwood@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1            # Number of CPU cores per task
#SBATCH --mem=8G
#SBATCH --time=12:00:00
#SBATCH --account=lien.nguyen

### Creating dictionary of necessary name items:

file="/blue/lien.nguyen/rossellwood/amplicon_pool_data_Huong/pool_samples_info/name_PKHD1SVA-pool_mixing-samples_detail.csv"

if [ ! -f "$file" ]; then
    echo "File '$file' does not exist or is not readable."
    exit 1
fi

# saving each part of the file  

cat "$file" | while IFS= read -r line; do
    #echo "$line"
    pool=$(echo $line | cut -d',' -f1)
    barcode=$(echo $line | cut -d',' -f2)
    sampleID=$(echo $line | cut -d',' -f3)
    cohort=$(echo $line | cut -d',' -f4)
    RU1=$(echo $line | cut -d',' -f5)
    RU2=$(echo $line | cut -d',' -f6)
    pattern=$(echo $line | cut -d',' -f7)
    #array+=("${tuple}")
done < "$file"

# adding things into the dictionary
declare -A dict

while IFS=, read -r pool barcode sampleID cohort RU1 RU2 pattern || [ -n "$pool" ]; do
    dict["$pool","$barcode"]="$pool, $barcode, $sampleID, $cohort, $RU1, $RU2, $pattern"
done < "$file"

# Now you can use dict outside the loop
# Example: Print all keys and values in dict
for key in "${!dict[@]}"; do
    echo "Key: $key,--> Value: ${dict[$key]}"
done

### Now changing file names:

# cd to demux FASTQ folder
#cd /blue/lien.nguyen/rossellwood/amplicon_pool_data_Huong/demux_out/demux_fastq/demux_fastq_perSample  


###############################
###############################

# Change directory to demux FASTQ folder
# cd /blue/lien.nguyen/rossellwood/amplicon_pool_data_Huong/demux_out/demux_fastq/TEST_FASTQ_dir

# # Loop through each file in the directory
# ls -1 | while read -r file; do
#     if [ -f "$file" ]; then
#         # Loop through each key-value pair in the dictionary
#         for sampleID in "${!dict[@]}"; do
#             info="${dict[$sampleID]}"  # Get the stored info for this sampleID
#             pool=$(echo $line | cut -d',' -f1)
#             barcode=$(echo $line | cut -d',' -f2)
#             sampleID=$(echo $line | cut -d',' -f3)
#             cohort=$(echo $line | cut -d',' -f4)
#             RU1=$(echo $line | cut -d',' -f5)
#             RU2=$(echo $line | cut -d',' -f6)
#             pattern=$(echo $line | cut -d',' -f7)
            
#             # Check if the filename contains both $barcode and $sampleID
#             if [[ "$file" == *"$barcode"* && "$file" == *"$sampleID"* ]]; then
#                 # Rename the file using the information from the dictionary
#                 new_filename="${pool}_${barcode}_${sampleID}_${cohort}_${RU1}_${RU2}_${pattern}.fastq"
#                 mv "$file" "$new_filename"
#                 echo "Renamed '$file' to '$new_filename'"
#             fi
#         done
#     fi
# done


#################################
#################################

# go through directory one file at a time
# ls -1 | while read file; do
    
#     if [ -f "$file" ]; then
#         # Loop through each key-value pair in the dictionary
#         for sampleID in "${!dict[@]}"; do
#             barcode="${dict[$sampleID]}"
#             # Check if the filename contains the pattern "$sampleID==$barcode"
#             if [[ "$file" == *"$barcode"* && "$file" == *"$pool"* ]]; then
#                 echo "yes"
#             fi
#         done



#         # if [[ "$file" == *"$barcode"* && "$file" == *"$pool"* ]]; then
#         #     echo "File '$file' contains '$pool'"
#         # fi
#     fi
#     #echo $file 
# done









# # Loop through each file in the directory
# ls -1 | while read -r file; do
#     if [ -f "$file" ]; then
#         # Loop through each key-value pair in the dictionary
#         for sampleID in "${!dict[@]}"; do
#             barcode="${dict[$sampleID]}"
#             # Check if the filename contains both $barcode and $sampleID
#             if [[ "$file" == *"$barcode"* && "$file" == *"$sampleID"* ]]; then
#                 echo "File '$file' contains both '$barcode' and '$sampleID'"
#             fi
#         done
#     fi
# done
# if pool # and barcode = pool # and barcode in dictionary, rename the file with the proper information















# Read file line by line and save to dictionary
# while IFS= read -r line; do
#     echo $line | cut -d',' -f3
#     key=$(echo "$line" | cut -d',' -f2)  # Assuming lines are key:value format
#     value=$(echo "$line" | cut -d',' -f1)
#     dict["$key"]+="$value"
# done < "$file"

# # Print the contents of the dictionary (for verification)
# echo "Dictionary Contents:"
# for key in "${!dict[@]}"; do
#     echo "Key: $key, Value: ${dict[$key]}"
# done






