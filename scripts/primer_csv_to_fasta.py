#!/usr/bin/env python

import argparse
import re

def get_args():
    parser = argparse.ArgumentParser(description="A program to hold input + output file name")
    parser.add_argument("-f", "--file", help="designates absolute file path to input csv file", type = str)
    parser.add_argument("-o", "--outFasta", help="designates absolute file path to output fasta file", type = str)
    return parser.parse_args()
    
args = get_args()

# This script takes primers from a csv and turns them into FASTA format

seen_headers = set()

with open(args.file) as csv_file, open(args.outFasta, "w") as fasta:
    while True:
        line = csv_file.readline().strip()
        if(line == ""):
            break
        if line.startswith("Pool"):
             fasta.write(f'{line}\n')
        
        if line.startswith("PK1"):
           # break into parts then print out each part in specific manner
            parts = line.split(',')
            

            header = parts[0]
            seq = parts[1]

            #only keep the unique headers
            #if header not in seen_headers:
            seen_headers.add(header)
            fasta.write(f'>{header}\n{seq}\n')
            
###
# to run: ./primer_csv_to_fasta.py -f ../Fwd_Pacbioampliconsequencing/PKHD1SVA-amplicon\ sequencing-submitting.csv -o ../Fwd_Pacbioampliconsequencing/barcodes.fasta
###