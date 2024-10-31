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

##################################################
# Functions
##################################################

# Function to filter abnormal reads
def filter_abnormal_reads(reference, seq, threshold, header, plus, qual):
    '''This function takes in the reference, threshold, and each read. It compares the sequnce against reference sequence. It only keeps reads above a certain fuzzy 
    score threshold.'''
    
    similarity = fuzz.ratio(reference, seq)
    if similarity >= threshold:
        filtered_reads.append(plus)
        filtered_fq.write(f"{header}\n{seq}\n{plus}\n{qual}\n")
    return filtered_reads

# def average_read_len(seq):
#     '''this function will return the average read length for'''
#     seqLen = (len(seq))
#     return seqLen


# def generate_reference_seq(average_len):
#     '''This function will input the average sequence length, and will extract the first sequence that has this length. This sequence will be the reference 
#     sequence in which other sequences will compare strings against'''

#     return reference

###############################################
# Globals
###############################################
# Define a similarity threshold (e.g., 70% similarity)
threshold = 60
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
        if len(seq) > 1500:
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
# print(reference)

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





















# #!/usr/bin/env python

# import argparse
# import re
# import os
# from fuzzywuzzy import fuzz
# from fuzzywuzzy import process

# def get_args():
#     parser = argparse.ArgumentParser(description="A program to hold input + output file name")
#     parser.add_argument("-f", "--file", help="designates absolute file path to the input fastq file", type = str)
#     parser.add_argument("-o", "--filtered_fq", help="designates absolute file path to the filtered FQ file", type = str)    
#     return parser.parse_args()
    
# args = get_args()

# ### conda activate fuzzywuzzy ###

# ##################################################
# # Functions
# ##################################################

# # Function to filter abnormal reads
# def filter_abnormal_reads(reference, seq, threshold, header, plus, qual):
#     '''This function takes in the reference, threshold, and each read. It compares the sequnce against reference sequence. It only keeps reads above a certain fuzzy 
#     score threshold.'''
    
#     similarity = fuzz.ratio(reference, seq)
#     if similarity >= threshold:
#         filtered_reads.append(plus)
#         filtered_fq.write(f"{header}\n{seq}\n{plus}\n{qual}\n")
#     return filtered_reads

# # def average_read_len(seq):
# #     '''this function will return the average read length for'''
# #     seqLen = (len(seq))
# #     return seqLen


# # def generate_reference_seq(average_len):
# #     '''This function will input the average sequence length, and will extract the first sequence that has this length. This sequence will be the reference 
# #     sequence in which other sequences will compare strings against'''

# #     return reference

# ###############################################
# # Globals
# ###############################################
# # Define a similarity threshold (e.g., 70% similarity)
# threshold = 60
# filtered_reads = []

# nonexp_reference = ""
# exp_reference = ""

# nonexp_seqLenDict = dict()
# exp_seqLenDict = dict()

# ##############################################
# # Main
# ##############################################

# ###first pass: determine the most abundant read length, split up by allele
# with open(args.file, "r") as in_fastq, open(args.filtered_fq, "w") as filtered_fq:
#     print(in_fastq)
#     while True:
#         header=in_fastq.readline().strip()
#         if header == "":              #break if detects blank new line
#             break 
#         seq=in_fastq.readline().strip()
#         plus=in_fastq.readline().strip()
#         qual=in_fastq.readline().strip()

#         # ReGex to first separate Hexamer expansion or nonexpansion as well as interruption cases
#         match = re.search("([ATGC]*GTGAGGCGTAG)([ATGC]*CTCTGTGCAATCGGAGTAGAGG)", seq) #hexamer
#         match2 = re.search("(CCTCTACTCCGATTGCACAGAG[ATGC]*)(CTACGCCTCAC[AGTC]*)", seq) #reverse complement of first regex statement hexamer

#         if match:
#             hexamer = match.group(2)
#             VNTR = match.group(1) 

#             if len(hexamer) < 403:
#             #determine most abundant sequence length:
#                 seqLen = (len(seq))
#                 if seqLen not in nonexp_seqLenDict:
#                     nonexp_seqLenDict[seqLen] = 1
#                 else:
#                     nonexp_seqLenDict[seqLen] += 1
#             if len(hexamer) >= 403:
#                 seqLen = (len(seq))
#                 if seqLen not in exp_seqLenDict:
#                     exp_seqLenDict[seqLen] = 1
#                 else:
#                     exp_seqLenDict[seqLen] += 1
# nonexp_most_abundant_len = max(nonexp_seqLenDict, key=nonexp_seqLenDict.get)
# exp_most_abundant_len = max(exp_seqLenDict, key=exp_seqLenDict.get)

# print(f'{nonexp_most_abundant_len=}')
# print(f'{exp_most_abundant_len=}')

# ###second pass: find the first read with the most abundant seq length
# with open(args.file, "r") as in_fastq, open(args.filtered_fq, "w") as filtered_fq:
#     while True:
#         header=in_fastq.readline().strip()
#         if header == "":              #break if detects blank new line
#             break 
#         seq=in_fastq.readline().strip()
#         plus=in_fastq.readline().strip()
#         qual=in_fastq.readline().strip()

#         #generate reference for nonexp:
#         if len(seq) == nonexp_most_abundant_len:
#             nonexp_reference = seq
#             continue
#         #generate reference for exp:
#         if len(seq) == exp_most_abundant_len:
#             exp_reference = seq
#             break
#     print("________")
#     print(f'{nonexp_reference=}')
#     print("__________")
#     print(exp_reference)
#     print("__________")

# ###third pass: use fuzzywuzzy to filter reads and generate new fastq file
# with open(args.file, "r") as in_fastq, open(args.filtered_fq, "w") as filtered_fq:
#     while True:
#         header=in_fastq.readline().strip()
#         if header == "":              #break if detects blank new line
#             break 
#         seq=in_fastq.readline().strip()
#         plus=in_fastq.readline().strip()
#         qual=in_fastq.readline().strip()
#         #filter reads that don't meet fuzzwuzzy threshold
#         filter_abnormal_reads(reference, seq, threshold, header, plus, qual)

# print(len(filtered_reads))

# print(most_abundant_len)
# print(reference)

# #  ./fuzzywuzzy_str_comparison.py -f /blue/lien.nguyen/rossellwood/amplicon_pool_data_Huong/demux_out/demux_fastq/demux_fastq_perSample/pool-1_1054_PK1sva_bc1010_Female_Control_19R1_53R2_p1.fastq










