# PKHD1 SVA Computational Analysis

This repository contains all of the scripts used during the analysis of the PKHD1 SVA. The scripts span from the downloading of data to the final data analysis. 


## Project WorkFlow
#### Data Generation
* Began by downloading raw FASTQ data from Azure blob storage

#### Pre-Processing
* Once raw FASTQ data was downloaded, [lima_demux.sh](https://github.com/roel1289/PKHD1_SVA_2024/blob/main/scripts/lima_demux.sh) was used to demultiplex data
* [primer_csv_to_fasta.py](https://github.com/roel1289/PKHD1_SVA_2024/blob/main/scripts/primer_csv_to_fasta.py) as used to rename the samples
* Next, to verify the proper region of the genome, [minimap2.sh](https://github.com/roel1289/PKHD1_SVA_2024/blob/main/scripts/minimap2.sh) was used to align reads to hg19, GCHr38, and T2T
* Then, using the FASTQ files again, [fuzzywuzzy_str_comparison.py](https://github.com/roel1289/PKHD1_SVA_2024/blob/main/scripts/fuzzywuzzy_str_comparison.py) and [run_fuzzywuzzy.sh](https://github.com/roel1289/PKHD1_SVA_2024/blob/main/scripts/run_fuzzywuzzy.sh) were used to filter out truncated reads, leaving us with only the reads that contain a full SVA


#### Data Analysis
* Once data was available, [allele_call.py](https://github.com/roel1289/PKHD1_SVA_2024/blob/main/scripts/allele_call.py) and [run_allele_call.sh](https://github.com/roel1289/PKHD1_SVA_2024/blob/main/scripts/run_allele_call.sh) were used to call allele types for each sample
* To visualize the data, various R scripts were used, such as [PKHD1_allele_calling.Rmd](https://github.com/roel1289/PKHD1_SVA_2024/blob/main/scripts/PKHD1_allele_calling.Rmd), [PKHD1_stats.Rmd](https://github.com/roel1289/PKHD1_SVA_2024/blob/main/scripts/PKHD1_stats.Rmd), and [PKHD1_Summary.Rmd](https://github.com/roel1289/PKHD1_SVA_2024/blob/main/scripts/PKHD1_Summary.Rmd). [PKHD1_Summary.Rmd](https://github.com/roel1289/PKHD1_SVA_2024/blob/main/scripts/PKHD1_Summary.Rmd) was the script used to make the final report. 


*note:* there are various other scripts in this repository used along the way to learn more about the sequences and SVA retrotransposon. 