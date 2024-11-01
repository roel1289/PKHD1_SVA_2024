#!/bin/bash
#SBATCH --job-name=minimap2
#SBATCH --mail-user=rossellwood@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8            # Number of CPU cores per task
#SBATCH --mem=10G
#SBATCH --time=12:00:00
#SBATCH --account=lien.nguyen
#SBATCH --qos=lien.nguyen-b

# conda activate data_processing  (contains minimap2 and samtools)  


# This script aligns FASTQ files to a reference genome
# When running with a new reference genome, all you have to do is change the ref genome file path, and 
# which ref genome you will be aligning to


#paths:
fastq_demux="/blue/lien.nguyen/rossellwood/amplicon_pool_data_Huong/demux_out/demux_fastq/demux_fastq_perSample"
hg19_ref="/blue/lien.nguyen/hg19_ref_minimap/hg19.mmi"
GRCh38_ref="/blue/lien.nguyen/rossellwood/humanGenomeAssembly/Homo_sapiens.GRCh38.dna.primary_assembly.mmi"
T2T_ref="/blue/lien.nguyen/rossellwood/T2T_assembly/T2T.mmi"
bam_files="/blue/lien.nguyen/rossellwood/amplicon_pool_data_Huong/bam_files/bam_files_T2T"
  
# aligning to the ref genome
cd $fastq_demux

files=$(ls -1 *.fastq)

for file in $files
do
    #align file
    echo "mapping: $file"
    /usr/bin/time -v minimap2 -a -t 8 $T2T_ref $file > ../../../bam_files/bam_files_T2T/$file.sam
    #sort and conver to bam
    echo "converting to bam: $file"
    /usr/bin/time -v samtools sort ../../../bam_files/bam_files_T2T/$file.sam -o ../../../bam_files/bam_files_T2T/$file.bam

done


# renaming the files so that they don't have "demux.fastq", and then indexing the BAM file
cd $bam_files
for file in $(ls -1 *.bam); do
    # rename file so it doesn't have an extra ".fastq." in it
    echo "renaming $file"
    mv "$file" "${file//.fastq/}"
done

#indexing file
for file in $(ls -1 *.bam); do
    echo "now indexing BAM: $file"
    /usr/bin/time -v samtools index $file
done




#remove all sam files to save more space:
rm *.sam


##### other samtools commands:
# samtools view -bS 185p95_1230GR.sam | samtools sort -o 185p95_1230GR.sorted.bam
# samtools index 185p95_1230GR.sorted.bam

#can also do it all in one go: bwa mem genome.fa reads.fastq | samtools sort -o myfile_sorted.bam