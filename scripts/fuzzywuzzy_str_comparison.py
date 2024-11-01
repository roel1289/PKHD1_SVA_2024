#!/usr/bin/env python

import argparse
import re
import os
from fuzzywuzzy import fuzz
from fuzzywuzzy import process

def get_args():
    parser = argparse.ArgumentParser(description="A program to hold input + output file name")
    parser.add_argument("-f", "--file", help="designates absolute file path to the input fastq file", type = str)
    parser.add_argument("-o", "--filtered_fq", help="designates absolute file path to the filtered FQ file", type = str)    
    return parser.parse_args()
    
args = get_args()

### conda activate fuzzywuzzy ###

# This script is used to take a raw fastq file, compare reads, and filter out any reads that are abnormal. In my project,
# this file was used to filter our non-full-SVA containing reads

##################################################
# Functions
##################################################

# Function to filter abnormal reads
def filter_abnormal_reads(reference, seq, threshold, header, plus, qual):
    '''This function takes in the reference, threshold, and each read. It compares the sequnce against reference sequence. It only keeps reads above a certain fuzzy 
    score threshold.'''
    
    similarity = fuzz.ratio(reference, seq)
    if similarity >= threshold and ("GGGAGA" in seq or "TCTCCC" in seq):
        filtered_reads.append(plus)
        filtered_fq.write(f"{header}\n{seq}\n{plus}\n{qual}\n")
    return filtered_reads



###############################################
# Globals
###############################################
# Define a similarity threshold (e.g., 70% similarity)
threshold = 77
filtered_reads = []

seqLenDict = dict()

##############################################
# Main
##############################################

###first pass: determine the most abundant read length
with open(args.file, "r") as in_fastq, open(args.filtered_fq, "w") as filtered_fq:
    while True:
        header=in_fastq.readline().strip()
        if header == "":              #break if detects blank new line
            break 
        seq=in_fastq.readline().strip()
        plus=in_fastq.readline().strip()
        qual=in_fastq.readline().strip()

        # any read that is full SVA is above around 2000 bp. Doing this because there are many truncated reads located in the end of the files:
        if len(seq) > 2000 and ("GGGAGAGGGAGA" in seq or "TCTCCCTCTCCC" in seq):
        #determine most abundant sequence length:
            seqLen = (len(seq))
            if seqLen not in seqLenDict:
                seqLenDict[seqLen] = 1
            else:
                seqLenDict[seqLen] += 1
most_abundant_len = max(seqLenDict, key=seqLenDict.get)

###second pass: find the first read with the most abundant seq length
with open(args.file, "r") as in_fastq, open(args.filtered_fq, "w") as filtered_fq:
    while True:
        header=in_fastq.readline().strip()
        if header == "":              #break if detects blank new line
            break 
        seq=in_fastq.readline().strip()
        plus=in_fastq.readline().strip()
        qual=in_fastq.readline().strip()

        #generate reference:
        if len(seq) == most_abundant_len:
            reference = seq
            break

###third pass: use fuzzywuzzy to filter reads and generate new fastq file
with open(args.file, "r") as in_fastq, open(args.filtered_fq, "w") as filtered_fq:
    while True:
        header=in_fastq.readline().strip()
        if header == "":              #break if detects blank new line
            break 
        seq=in_fastq.readline().strip()
        plus=in_fastq.readline().strip()
        qual=in_fastq.readline().strip()
        #filter reads that don't meet fuzzwuzzy threshold
        filter_abnormal_reads(reference, seq, threshold, header, plus, qual)

print(f'{len(filtered_reads)=}')

print(f'{most_abundant_len=}')
print(f'{reference=}')

#  ./fuzzywuzzy_str_comparison.py -f /blue/lien.nguyen/rossellwood/amplicon_pool_data_Huong/demux_out/demux_fastq/demux_fastq_perSample/pool-1_1054_PK1sva_bc1010_Female_Control_19R1_53R2_p1.fastq

