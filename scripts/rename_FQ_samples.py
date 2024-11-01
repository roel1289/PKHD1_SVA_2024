#!/usr/bin/env python

import argparse
import re
import os

def get_args():
    parser = argparse.ArgumentParser(description="A program to hold input + output file name")
    parser.add_argument("-f", "--file", help="designates absolute file path to input csv file", type = str)
    parser.add_argument("-d", "--FQ_dir", help="designates absolute file path to the folder containing the FASTQ files", type = str)
    return parser.parse_args()
    
args = get_args()

# This script renames file names using an excel file. This way the file names have meaning to them in terms of repeat length and AD vs Control


# making a set that contains tuples (info for each part)
name_set = set()


with open(args.file) as name_CSV:
    while True:
        line = name_CSV.readline().strip()
        if(line == ""):
            break
        #print(line)
        parts = line.split(',')

        # defining parts of the csv:
        sex = parts[0]
        pool = parts[1]
        barcode = parts[2]
        sampleID = parts[3]
        cohort = parts[4]
        RU1 = parts[5]
        RU2 = parts[6]
        pattern = parts[7]
 
        if (sex, pool, barcode, sampleID, cohort, RU1, RU2, pattern) not in name_set:
            name_set.add(tuple((sex, pool, barcode, sampleID, cohort, RU1, RU2, pattern)))


print(name_set)


#change directory to the folder that has all of the FASTQ files:  
os.chdir(args.FQ_dir)

#print working directory
print(f'FASTQ dir:{os.getcwd()}')

#changing the file names:
for file in os.listdir("."):
    if file.endswith('.fastq'):
        for sex, pool, barcode, sampleID, cohort, RU1, RU2, pattern in name_set:
            #print(pool)
            #print(file)
            if f'pool-{pool}' in file and barcode in file:
                #print(f'this is the pool #: {pool} and file: {file} and barcode: {barcode}')
            
                # renaming files
                new_name = f'pool-{pool}_{sampleID}_{barcode}_{sex}_{cohort}_{RU1}R1_{RU2}R2_p{pattern}.fastq'
                os.rename(file, new_name)
                print(f'Renamed {file} ---> {new_name}')






# ./rename_FQ_samples.py -f /blue/lien.nguyen/rossellwood/amplicon_pool_data_Huong/pool_samples_info/name_PKHD1SVA-pool_mixing-samples_detail.csv -d /blue/lien.nguyen/rossellwood/amplicon_pool_data_Huong/demux_out/demux_fastq/TEST_FASTQ_dir