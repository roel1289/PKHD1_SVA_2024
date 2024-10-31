#!/usr/bin/env python

import argparse
import re
import os

def get_args():
    parser = argparse.ArgumentParser(description="A program to hold input + output file name")
    parser.add_argument("-f", "--file", help="designates absolute file path to mTR text file", type = str)
    parser.add_argument("-o", "--out", help="designates absolute file path to filtered mTR text output", type = str)
    return parser.parse_args()
    
args = get_args()

#dict. Key = sequence, value is the amount of times this sequence is seen
seq_dict = dict()


with open(args.file, "r") as mTR_in:
    while True:
        line=mTR_in.readline().strip()
        if line == "":              #break if detects blank new line
            break 
        
        parts = line.split('\t')
        read_id = str(parts[0])
        read_len = int(parts[1])
        start= int(parts[2])
        end = int(parts[3])
        length = int(parts[4])
        unit_length = int(parts[5])
        unit_occurrences = int(parts[6])
        ratio = float(parts[8])
        seq = str(parts[12])

        #print(start, end, length)
        #print(f'verify ---> {end - start + 1} for {read_id} occur {unit_occurrences=}')
        #print(f'{len(seq)} ---> {unit_occurrences=}')


        # if float(ratio) > float(0.90):
        #     print("test")

################################################ Get VNTR length:
        #focus on longest VNTR region 
        # if  unit_length > 10 and length > 800:
        #     #print(f'{start=}, {end=}, {length=}, {unit_occurrences=}, {seq=}')
        #     #print(f'{read_id}, {read_len=},{length=}, and {seq=}')
        #     print(line)
#################################################

        # find out which sequences are most common:
        if unit_length > 10:
            if start in seq_dict:
                seq_dict[start] += 1
            else:
                seq_dict[start] = 1
        


        if start == 764:
            print(start, end)

        # if seq_dict[seq] > 6000:
        #     print(seq_dict[seq])
#confirm the most common repeat sequence
for item in seq_dict:
    if seq_dict[item] > 6340:
        x = item, seq_dict[item]
        #print(seq_dict[item])
        print(f'{x=}')

# print(seq_dict)