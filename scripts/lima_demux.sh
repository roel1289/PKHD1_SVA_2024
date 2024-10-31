#!/bin/bash
#SBATCH --job-name=lima_demux
#SBATCH --mail-user=rossellwood@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=3            # Number of CPU cores per task
#SBATCH --mem=12G
#SBATCH --time=12:00:00
#SBATCH --account=lien.nguyen

#note: use the first section of code when the input is BAM files, and use the second section of code when the input is FASTQ


##########################
### When demultiplexing using BAM files:
##########################
# bam_dir="/blue/lien.nguyen/rossellwood/amplicon_pool_data_Huong/stc778-huong-nguyen-amplicon" 
# barcodes="/blue/lien.nguyen/rossellwood/amplicon_pool_data_Huong/Fwd_Pacbioampliconsequencing/barcodes.fasta"


# cd $bam_dir

# files=$(ls -1 *.bam)

# for file in $files
# do
#     echo $file
#     /usr/bin/time -v lima $file $barcodes ../demux_out/$file.demux.bam --hifi-preset SYMMETRIC
# done
#############################

##########################################
### When demultiplexing using FASTQ files: 
##########################################

fastq_dir="/blue/lien.nguyen/rossellwood/amplicon_pool_data_Huong/fastq_files" 
barcode_dir="/blue/lien.nguyen/rossellwood/amplicon_pool_data_Huong/Fwd_Pacbioampliconsequencing/barcodes_dir"


#looping through the barcode directory at the same timea s looping through the fastq directory. 
#We will need the matching pool number for each file  

#counter to make sure everything is iterated in the same way:
counter1=0
for fastq in $(ls -1 $fastq_dir/);do
    counter1=$((counter1+1))
    counter2=0
    for barcode_fa in $(ls -1 $barcode_dir/);do
        counter2=$((counter2+1))
        if [[ $counter1 -eq $counter2 ]]; then
        echo "fastq ----> $fastq"
        echo "barcode ---->$barcode_fa"
        /usr/bin/time -v lima -j 3 --split-named  \
            /blue/lien.nguyen/rossellwood/amplicon_pool_data_Huong/fastq_files/$fastq \
            /blue/lien.nguyen/rossellwood/amplicon_pool_data_Huong/Fwd_Pacbioampliconsequencing/barcodes_dir/$barcode_fa \
            /blue/lien.nguyen/rossellwood/amplicon_pool_data_Huong/demux_out/demux_fastq/TEST_FASTQ_dir/$fastq.demux.fastq --hifi-preset SYMMETRIC

        fi
    done
done












# files=$(ls -1 *.fastq)

# for file in $files
# do
#     echo $file
#     /usr/bin/time -v lima -j 3 --split-named $file $barcodes ../demux_out/demux_fastq/TEST_FASTQ_dir/$file.demux.fastq --hifi-preset SYMMETRIC
# done