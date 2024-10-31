#!/bin/bash
#SBATCH --job-name=fuzzywuzzy
#SBATCH --mail-user=rossellwood@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=3            # Number of CPU cores per task
#SBATCH --mem=10G
#SBATCH --time=12:00:00
#SBATCH --account=lien.nguyen
#SBATCH --qos=lien.nguyen-b

# conda activate fuzzywuzzy to run script


fastq_dir="/blue/lien.nguyen/rossellwood/amplicon_pool_data_Huong/demux_out/demux_fastq/demux_fastq_perSample"  
#TEST:
# fastq_dir="/blue/lien.nguyen/rossellwood/TEST"

cd $fastq_dir

files=$(ls -1 *.fastq)

for file in $files
do
    echo "filter file: $file"
    #concatenating all lengths into one file
    
    /usr/bin/time -v /blue/lien.nguyen/rossellwood/repeat_expansion_tools/sort_expansion_type/fuzzywuzzy_str_comparison.py -f $file -o /blue/lien.nguyen/rossellwood/repeat_expansion_tools/sort_expansion_type/filtered_full_FASTQ/${file}_filtered.fastq
done