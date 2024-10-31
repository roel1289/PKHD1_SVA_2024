#!/usr/bin/env python

import argparse
import re
import os
from fuzzywuzzy import fuzz
from fuzzywuzzy import process
from fuzzysearch import find_near_matches

def get_args():
    parser = argparse.ArgumentParser(description="A program to hold input + output file name")
    parser.add_argument("-f", "--in_file", help="designates absolute file path to the input file that contains the FASTQ files", type = str)
    parser.add_argument("-o", "--out_table", help="designates absolute file path to the output allele table. It will be a .txt", type = str)
    
    return parser.parse_args()
    
args = get_args()

#hexamer-hexvar
SHex_hexv0 = 0
SHex_hexv1 = 0
SHex_hexv2 = 0
SHex_hexv3 = 0
SHex_hexv4 = 0
SHex_hexv5 = 0
Shex_hexv6 = 0
LHex_hexv0 = 0
LHex_hexv1 = 0
LHex_hexv2 = 0
LHex_hexv3 = 0
LHex_hexv4 = 0
LHex_hexv5 = 0
Lhex_hexv6 = 0


SHex_total_ratio = 0
LHex_total_ratio = 0
total_reads_SHex = 0
total_reads_LHex = 0
    
with open(args.in_file, "r") as in_fastq:
    while True:
        header=in_fastq.readline().strip()
        if header == "":              #break if detects blank new line
            break 
        seq=in_fastq.readline().strip()
        plus=in_fastq.readline().strip()
        qual=in_fastq.readline().strip()         

        ### REGEX ###
        G4AGA_match = re.finditer("GGGGAGA", seq)
        GGGAGA_match = re.finditer("AGGGAGA", seq)
        G4AGA_rev_match = re.search("TCTCCCC", seq)
        GGGAGA_rev_match = re.search("CTCTCCC", seq)


        G4AGA_count = 0
        GGGAGA_count = 0
        Shex_ratio = 0
        Lhex_ratio = 0

        # ReGex to first separate Hexamer expansion or nonexpansion as well as interruption cases
        match = re.search("([ATGC]*GTGAGGCGTAG)([ATGC]*CTCTGTGCAATCGGAGTAGAGG)", seq) 
        match2 = re.search("(CCTCTACTCCGATTGCACAGAG[ATGC]*)(CTACGCCTCAC[AGTC]*)", seq) #reverse complement of first regex statement   

        #ReGex for Hexamer variants
        hex_var0_match = re.search("GGGAGAGGGAGAGGGAGA", seq)
        hex_var1_match = re.search("[AG]+[C|G|A]TTTG[G]+AG[C]+TCTC[T]+AA[A]*GG", seq)
        hex_var2_match = re.search("[AG]+([C|G]TTTGGAGCCTCTCTAAAAGGGCTCTGTACC)[AG]+(GGGAGA){3,}[GA]+([C|G]TTTGGAGCCTCTCTAAAAGGGCTCTGTACC)", seq)
        hex_var3_match = re.search("[AG]+[C|G]TTTGGAGCCTCTCTAAAAGGGCTCTGTACC[AG]+(G{5}A){2,}[GA]+([C|G]TTTGGAGCCTCTCTAAAAGGGCTCTGTACC)", seq) #use se
        hex_var4_match = re.search("[AG]+CGGAGACGGAGA[AG]+", seq)
        hex_var5_match = re.search("GGGGAGAGGGGAGA[GA]+[C|G][AG]+[C|G]TTTGGAGCCTCTCTAA[A]*GG", seq)
        hex_var6_match = re.search("GGGGAGAGGGGAGAGGGGAGA", seq)  

        
        if match:
            hexamer = match.group(2)
            VNTR = match.group(1)
            
            #short hexamer
            if len(hexamer) < 403:
                if hex_var2_match:
                    SHex_hexv2+=1
                elif hex_var3_match:
                    SHex_hexv3+=1
                elif hex_var4_match:
                    SHex_hexv4+=1
                elif hex_var5_match:
                    SHex_hexv5+=1
                    # for m in G4AGA_match:
                    #     print(m)
                    #     G4AGA_count+=1
                    # for m in GGGAGA_match:
                    #     GGGAGA_count+=1
                    #     print(m)
                    # print(G4AGA_count+GGGAGA_count)
                elif hex_var6_match:
                    Shex_hexv6+=1
                elif hex_var1_match:
                    SHex_hexv1+=1
                elif hex_var0_match:
                    SHex_hexv0+=1

                # print((hexamer))

            #long hexamer
            if len(hexamer) >= 403:
                if hex_var2_match:
                    LHex_hexv2+=1
                elif hex_var3_match:
                    LHex_hexv3+=1
                elif hex_var4_match:
                    LHex_hexv4+=1
                elif hex_var5_match:
                    LHex_hexv5+=1
                    # for m in G4AGA_match:
                    #     # print(m)
                    #     G4AGA_count+=1
                    # for m in GGGAGA_match:
                    #     GGGAGA_count+=1
                    #     # print(m)
                    # print(G4AGA_count+GGGAGA_count)


                    print((len(hexamer)- 229)/6)
                elif hex_var6_match:
                    Lhex_hexv6+=1
                elif hex_var1_match:
                    LHex_hexv1+=1
                elif hex_var0_match:
                    LHex_hexv0+=1

                # print(len(hexamer))