#!/usr/bin/env python

import argparse
import re
import os
from fuzzywuzzy import fuzz
from fuzzywuzzy import process
from fuzzysearch import find_near_matches

def get_args():
    parser = argparse.ArgumentParser(description="A program to hold input + output file name")
    parser.add_argument("-d", "--in_dir", help="designates absolute file path to the input directory that contains the FASTQ files", type = str)
    parser.add_argument("-o", "--out_table", help="designates absolute file path to the output allele table. It will be a .txt", type = str)
    
    return parser.parse_args()
    
args = get_args()


# This script will call alleles. It will begin by calling alleles based off of hexamer and VNTR lengths, inwhich there are 4 combinations
# If needed, I will dive more deeply distinguishing between interruption types
# This file will input a FASTQ file and will output the 2 most abundant allele types. Calling based on have the largest two counts. 
# In future versions of this file, I will try and make it more concise and object oriented. 


####################################
### MAIN ###
######################################

os.chdir(args.in_dir)
print(f'FASTQ dir:{os.getcwd()}')
# #name directory

for file in os.listdir("."):
    if file.endswith('.fastq'):
        
        
        print(file)
        #grouping files into +/-, +/+, and -/-:
        filename_match = re.search("([0-9]+)R1_([0-9]+)R2", file)
        if filename_match:
            R1 = int(filename_match.group(1))
            R2 = int(filename_match.group(2))

            print(f'{R1}----------{R2}')
            if R1 < 29 and R2 >= 29:
                print(file)



                print(f'Calling +/- allele for: {file}')
                
                with open(file, "r") as in_fastq:
                    #hexamer counts
                    short_hex = 0
                    long_hex = 0
                    short_VNTR = 0
                    long_VNTR = 0

                    #hexamer-VNTR counts
                    S_short_VNTR = 0
                    S_long_VNTR = 0
                    L_short_VNTR = 0
                    L_long_VNTR = 0

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

                    #hexamer-VNTR-hexvar counts
                    SHex_SVNTR_hexv0 = 0
                    SHex_SVNTR_hexv1 = 0
                    SHex_SVNTR_hexv2 = 0
                    SHex_SVNTR_hexv3 = 0
                    SHex_SVNTR_hexv4 = 0
                    SHex_SVNTR_hexv5 = 0
                    SHex_SVNTR_hexv6 = 0
                    SHex_LVNTR_hexv0 = 0
                    SHex_LVNTR_hexv1 = 0
                    SHex_LVNTR_hexv2 = 0
                    SHex_LVNTR_hexv3 = 0
                    SHex_LVNTR_hexv4 = 0
                    SHex_LVNTR_hexv5 = 0
                    SHex_LVNTR_hexv6 = 0

                    LHex_SVNTR_hexv0 = 0
                    LHex_SVNTR_hexv1 = 0
                    LHex_SVNTR_hexv2 = 0
                    LHex_SVNTR_hexv3 = 0
                    LHex_SVNTR_hexv4 = 0
                    LHex_SVNTR_hexv5 = 0
                    LHex_SVNTR_hexv6 = 0
                    LHex_LVNTR_hexv0 = 0
                    LHex_LVNTR_hexv1 = 0
                    LHex_LVNTR_hexv2 = 0
                    LHex_LVNTR_hexv3 = 0
                    LHex_LVNTR_hexv4 = 0
                    LHex_LVNTR_hexv5 = 0
                    LHex_LVNTR_hexv6 = 0

                    total_reads = 0
                    
                    while True:
                        header=in_fastq.readline().strip()
                        if header == "":              #break if detects blank new line
                            break 
                        seq=in_fastq.readline().strip()
                        plus=in_fastq.readline().strip()
                        qual=in_fastq.readline().strip()                       

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

                        #ReGex for hexamer variants that are REVERSE COMPLEMENT ###
                        rev_hex_var0_match = re.search("TCTCCCTCTCCCTCTCCC", seq)
                        rev_hex_var1_match = re.search("GGTACAGAGCCCTTTTAGAGAGGCTCCAAAGCTC[CT]+", seq)
                        rev_hex_var2_match = re.search("GGTACAGAGCCCTTTTAGAGAGGCTCCAAAGCTC[CT]+(TCTCCC){3,}[CT]+GGTACAGAGCCCTTTTAGAGAGGCTCCAAAGCTC[TC]+", seq)
                        rev_hex_var3_match = re.search("GGTACAGAGCCCTTTTAGAGAGGCTCCAAAGCTC[CT]+(TC{5}){2,}[CT]+GGTACAGAGCCCTTTTAGAGAGGCTCCAAAGCTC[TC]+", seq) #use se
                        rev_hex_var4_match = re.search("[TC]+TCTCCGTCTCCG[TC]+", seq)
                        rev_hex_var5_match = re.search("GGTACAGAGCCCTTTTAGAGAGGCTCCAAAGCTC[CT]+[TCTCCCCTCTCCCC]", seq)
                        rev_hex_var6_match = re.search("TCTCCCCTCTCCCCTCTCCCC", seq)

                        # fuzzy string matches of the different cases:
                        hex_var4 = "GGGAGACGGAGACGGAGAGGGAGA"
                        find_hexvar4 = find_near_matches(hex_var4, seq, max_l_dist = 1)
                        find_hexvar4_start = [match.start for match in find_hexvar4]
                        
                        
                        if match:
                            hexamer = match.group(2)
                            VNTR = match.group(1)
                            
                            #short hexamer
                            if len(hexamer) < 403:
                                short_hex += 1
                                if hex_var2_match:
                                    SHex_hexv2+=1
                                elif hex_var3_match:
                                    SHex_hexv3+=1
                                    print(f'hex3 seq; ---------------:\n {seq}')
                                elif hex_var4_match:
                                    SHex_hexv4+=1
                                elif hex_var5_match:
                                    SHex_hexv5+=1
                                    # print(f'hexv5: -----:\n {seq}')
                                elif hex_var6_match:
                                    Shex_hexv6+=1
                                    # print(f'hexv6: -----:\n {seq}')
                                elif hex_var1_match:
                                    SHex_hexv1+=1
                                elif hex_var0_match:
                                    SHex_hexv0+=1
                                    # print(f'hexv0: -----:\n {seq}')
                                    

                                if len(VNTR) < 2500: #short vntr
                                    S_short_VNTR +=1
                                    short_VNTR+=1
                                    
                                    if hex_var2_match:
                                        SHex_SVNTR_hexv2+=1
                                    elif hex_var3_match:
                                        SHex_SVNTR_hexv3+=1
                                    elif hex_var4_match:
                                        SHex_SVNTR_hexv4+=1
                                    elif hex_var5_match:
                                        SHex_SVNTR_hexv5+=1
                                    elif hex_var6_match:
                                        SHex_SVNTR_hexv6+=1
                                    elif hex_var1_match:
                                        SHex_SVNTR_hexv1+=1
                                    elif hex_var0_match:
                                        SHex_SVNTR_hexv0+=1
                                if len(VNTR) >= 2500: #long vntr
                                    S_long_VNTR +=1
                                    long_VNTR+=1   
                                    
                                    if hex_var2_match:
                                        SHex_LVNTR_hexv2+=1
                                    elif hex_var3_match:
                                        SHex_LVNTR_hexv3+=1
                                    elif hex_var4_match:
                                        SHex_LVNTR_hexv4+=1
                                    elif hex_var5_match:
                                        SHex_LVNTR_hexv5+=1
                                    elif hex_var6_match:
                                        SHex_LVNTR_hexv6+=1
                                    elif hex_var1_match:
                                        SHex_LVNTR_hexv1+=1 
                                    elif hex_var0_match:
                                        SHex_LVNTR_hexv0+=1  
                                

                            #long hexamer
                            if len(hexamer) >= 403:
                                long_hex += 1
                                if hex_var2_match:
                                    LHex_hexv2+=1
                                elif hex_var3_match:
                                    LHex_hexv3+=1
                                elif hex_var4_match:
                                    LHex_hexv4+=1
                                elif hex_var5_match:
                                    LHex_hexv5+=1
                                elif hex_var6_match:
                                    Lhex_hexv6+=1
                                elif hex_var1_match:
                                    LHex_hexv1+=1
                                elif hex_var0_match:
                                    LHex_hexv0+=1

                                if len(VNTR) < 2500: #short vntr
                                    L_short_VNTR +=1
                                    short_VNTR+=1
                                    
                                    if hex_var2_match:
                                        LHex_SVNTR_hexv2+=1
                                    elif hex_var3_match:
                                        LHex_SVNTR_hexv3+=1
                                    elif hex_var4_match:
                                        LHex_SVNTR_hexv4+=1
                                    elif hex_var5_match:
                                        LHex_SVNTR_hexv5+=1
                                    elif hex_var6_match:
                                        LHex_SVNTR_hexv6+=1
                                    elif hex_var1_match:
                                        LHex_SVNTR_hexv1+=1
                                    elif hex_var0_match:
                                        LHex_SVNTR_hexv0+=1

                                if len(VNTR) >= 2500: #long vntr
                                    L_long_VNTR +=1
                                    long_VNTR+=1
                                    
                                    if hex_var2_match:
                                        LHex_LVNTR_hexv2+=1
                                    elif hex_var3_match:
                                        LHex_LVNTR_hexv3+=1
                                    elif hex_var4_match:
                                        LHex_LVNTR_hexv4+=1
                                    elif hex_var5_match:
                                        LHex_LVNTR_hexv5+=1
                                    elif hex_var6_match:
                                        LHex_LVNTR_hexv6+=1
                                    elif hex_var1_match:
                                        LHex_LVNTR_hexv1+=1
                                    elif hex_var0_match:
                                        LHex_LVNTR_hexv0+=1


                        elif match2: #reverse complement:
                            hexamer = match2.group(1)
                            VNTR = match2.group(2)
                            #short hexamer
                            if len(hexamer) < 403:
                                short_hex += 1
                                if rev_hex_var2_match:
                                    SHex_hexv2+=1
                                elif rev_hex_var3_match:
                                    SHex_hexv3+=1
                                elif rev_hex_var4_match:
                                    SHex_hexv4+=1
                                elif rev_hex_var5_match:
                                    SHex_hexv5+=1
                                elif rev_hex_var6_match:
                                    Shex_hexv6+=1
                                elif rev_hex_var1_match:
                                    SHex_hexv1+=1
                                elif rev_hex_var0_match:
                                    SHex_hexv0+=1
                                
                                if len(VNTR) < 2500:
                                        S_short_VNTR +=1
                                        short_VNTR+=1
                                        
                                        if rev_hex_var2_match:
                                            SHex_SVNTR_hexv2+=1
                                        elif rev_hex_var3_match:
                                            SHex_SVNTR_hexv3+=1
                                        elif rev_hex_var4_match:
                                            SHex_SVNTR_hexv4+=1
                                        elif rev_hex_var5_match:
                                            SHex_SVNTR_hexv5+=1
                                        elif rev_hex_var6_match:
                                            SHex_SVNTR_hexv6+=1
                                        elif rev_hex_var1_match:
                                            SHex_SVNTR_hexv1+=1
                                        elif rev_hex_var0_match:
                                            SHex_SVNTR_hexv0+=1
                                if len(VNTR) >= 2500:
                                        S_long_VNTR +=1
                                        long_VNTR+=1
                                        
                                        if rev_hex_var2_match:
                                            SHex_LVNTR_hexv2+=1
                                        elif rev_hex_var3_match:
                                            SHex_LVNTR_hexv3+=1                                        
                                        elif rev_hex_var4_match:
                                            SHex_LVNTR_hexv4+=1
                                        elif rev_hex_var5_match:
                                            SHex_LVNTR_hexv5+=1
                                        elif rev_hex_var6_match:
                                            SHex_LVNTR_hexv6+=1
                                        elif rev_hex_var1_match:
                                            SHex_LVNTR_hexv1+=1
                                        elif rev_hex_var0_match:
                                            SHex_LVNTR_hexv0+=1
         
                            #long hexamer
                            if len(hexamer) >= 403:
                                long_hex += 1
                                if rev_hex_var2_match:
                                    LHex_hexv2+=1
                                elif rev_hex_var3_match:
                                    LHex_hexv3+=1
                                elif rev_hex_var4_match:
                                    LHex_hexv4+=1
                                elif rev_hex_var5_match:
                                    LHex_hexv5+=1
                                elif rev_hex_var6_match:
                                    Lhex_hexv6+=1
                                elif rev_hex_var1_match:
                                    LHex_hexv1+=1
                                elif rev_hex_var0_match:
                                    LHex_hexv0+=1
                                if len(VNTR) < 2500:
                                    L_short_VNTR +=1
                                    short_VNTR+=1
                                    
                                    if rev_hex_var2_match:
                                        LHex_SVNTR_hexv2+=1
                                    elif rev_hex_var3_match:
                                        LHex_SVNTR_hexv3+=1
                                    elif rev_hex_var4_match:
                                        LHex_SVNTR_hexv4+=1
                                    elif rev_hex_var5_match:
                                        LHex_SVNTR_hexv5+=1
                                    elif rev_hex_var6_match:
                                        LHex_SVNTR_hexv6+=1
                                    elif rev_hex_var1_match:
                                        LHex_SVNTR_hexv1+=1
                                    elif rev_hex_var0_match:
                                        LHex_SVNTR_hexv0+=1

                                if len(VNTR) >= 2500:
                                    L_long_VNTR +=1
                                    long_VNTR+=1
                                    
                                    if rev_hex_var2_match:
                                        LHex_LVNTR_hexv2+=1
                                    elif rev_hex_var3_match:
                                        LHex_LVNTR_hexv3+=1
                                    elif rev_hex_var4_match:
                                        LHex_LVNTR_hexv4+=1
                                    elif rev_hex_var5_match:
                                        LHex_LVNTR_hexv5+=1
                                    elif rev_hex_var6_match:
                                        LHex_LVNTR_hexv6+=1
                                    elif rev_hex_var1_match:
                                        LHex_LVNTR_hexv1+=1
                                    elif rev_hex_var0_match:
                                        LHex_LVNTR_hexv0+=1
                                    
                            
                        else: # use approximate string matching to determine lengths
                            primer1seq = "ACGAAAACCAGTGAGGCGTAG"
                            primer1 = find_near_matches(primer1seq,seq, max_l_dist=2)
                            # print(primer1)
                            start_position1 = [match.start for match in primer1]
                            if start_position1: #check if the position exists
                                start_position1 = start_position1[0]
                                VNTR = start_position1
                            
                            primer2seq = "CTCTGTGCAATCGGAGTAGAGG"
                            primer2 = find_near_matches(primer2seq,seq, max_l_dist=2)
                            # print(primer2)
                            start_position2 = [match.end for match in primer2]
                            if start_position2: #check if the position exists
                                start_position2 = start_position2[0]
                            # print(f'{start_position2=}')

                            hex_var1 = ""
                            hex_var2 = ""
                            hex_var3 = "GGGAGAGCTTTGGAGCCTCTCTAAAAGGGCTCTGTACCGGAGAGGGGAGAGGGGGAGGGGGAGGGGAGGGGGAGGGGGAGGGGGAGCTTTGGAGCCTCTCTAAAAGGGCTCTGTACC"
                            
                            hex_var4 = "GGGAGACGGAGACGGAGAGGGAGA"
                            find_hexvar4 = find_near_matches(hex_var4, seq, max_l_dist = 1)
                            find_hexvar4_start = [match.start for match in find_hexvar4]

                            if start_position1 and start_position2: 
                                # print(f'{start_position2 - start_position1}')
                            # print(f'Length between primers:{start_position2[0]-start_position1[0]}')
                                if (start_position2 - start_position1) < 403: #short hexamer
                                    short_hex += 1
                                    if hex_var2_match:
                                        SHex_hexv2+=1
                                    elif hex_var3_match:
                                        SHex_hexv3+=1
                                    elif hex_var4_match:
                                        SHex_hexv4+=1
                                    elif hex_var5_match:
                                        SHex_hexv5+=1
                                    elif hex_var6_match:
                                        Shex_hexv6+=1
                                    elif hex_var1_match:
                                        SHex_hexv1+=1
                                    elif hex_var0_match:
                                        SHex_hexv0+=1

                                    if (VNTR) < 2500:
                                        S_short_VNTR +=1
                                        short_VNTR+=1
                                        
                                        if hex_var2_match:
                                            SHex_SVNTR_hexv2+=1
                                        elif hex_var3_match:
                                            SHex_SVNTR_hexv3+=1
                                        elif hex_var4_match:
                                            SHex_SVNTR_hexv4+=1
                                        elif hex_var5_match:
                                            SHex_SVNTR_hexv5+=1
                                        elif hex_var6_match:
                                            SHex_SVNTR_hexv6+=1
                                        elif hex_var1_match:
                                            SHex_SVNTR_hexv1+=1
                                        elif hex_var0_match:
                                            SHex_SVNTR_hexv0+=1

                                    if (VNTR) >= 2500:
                                        S_long_VNTR +=1
                                        long_VNTR+=1
                                        
                                        if hex_var2_match:
                                            SHex_LVNTR_hexv2+=1
                                        elif hex_var3_match:
                                            SHex_LVNTR_hexv3+=1
                                        elif hex_var4_match:
                                            SHex_LVNTR_hexv4+=1
                                        elif hex_var5_match:
                                            SHex_LVNTR_hexv5+=1
                                        elif hex_var6_match:
                                            SHex_LVNTR_hexv6+=1
                                        elif hex_var1_match:
                                            SHex_LVNTR_hexv1+=1
                                        elif hex_var0_match:
                                            SHex_LVNTR_hexv0+=1
                                else: #long hexamer
                                    long_hex += 1
                                    if hex_var2_match:
                                        LHex_hexv2+=1
                                    elif hex_var3_match:
                                        LHex_hexv3+=1
                                    elif hex_var4_match:
                                        LHex_hexv4+=1
                                    elif hex_var5_match:
                                        LHex_hexv5+=1
                                    elif hex_var6_match:
                                        Lhex_hexv6+=1
                                    elif hex_var1_match:
                                        LHex_hexv1+=1
                                    elif hex_var0_match:
                                        LHex_hexv0+=1

                                    if (VNTR) < 2500:
                                        L_short_VNTR +=1
                                        short_VNTR+=1
                                        
                                        if hex_var2_match:
                                            LHex_SVNTR_hexv2+=1
                                        elif hex_var3_match:
                                            LHex_SVNTR_hexv3+=1
                                        elif hex_var4_match:
                                            LHex_SVNTR_hexv4+=1
                                        elif hex_var5_match:
                                            LHex_SVNTR_hexv5+=1
                                        elif hex_var6_match:
                                            LHex_SVNTR_hexv6+=1
                                        elif hex_var1_match:
                                            LHex_SVNTR_hexv1+=1
                                        elif hex_var0_match:
                                            LHex_SVNTR_hexv0+=1

                                    if (VNTR) >= 2500:
                                        L_long_VNTR +=1 
                                        long_VNTR+=1
                                        
                                        if hex_var2_match:
                                            LHex_LVNTR_hexv2+=1
                                        elif hex_var3_match:
                                            LHex_LVNTR_hexv3+=1
                                        elif hex_var4_match:
                                            LHex_LVNTR_hexv4+=1
                                        elif hex_var5_match:
                                            LHex_LVNTR_hexv5+=1
                                        elif hex_var6_match:
                                            LHex_LVNTR_hexv6+=1
                                        elif hex_var1_match:
                                            LHex_LVNTR_hexv1+=1
                                        elif hex_var0_match:
                                            LHex_LVNTR_hexv0+=1
                                        
                        total_reads += 1

#########################################################################################

            elif R1 >= 29 and R2 >= 29: ### in +/+ cases
                print("+/+", file)
                print(f'Calling +/+ allele for: {file}')
                
                with open(file, "r") as in_fastq:
                    #hexamer counts
                    short_hex = 0
                    long_hex = 0
                    short_VNTR = 0
                    long_VNTR = 0

                    #hexamer-VNTR counts
                    S_short_VNTR = 0
                    S_long_VNTR = 0
                    L_short_VNTR = 0
                    L_long_VNTR = 0

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

                    #hexamer-VNTR-hexvar counts
                    SHex_SVNTR_hexv0 = 0
                    SHex_SVNTR_hexv1 = 0
                    SHex_SVNTR_hexv2 = 0
                    SHex_SVNTR_hexv3 = 0
                    SHex_SVNTR_hexv4 = 0
                    SHex_SVNTR_hexv5 = 0
                    SHex_SVNTR_hexv6 = 0
                    SHex_LVNTR_hexv0 = 0
                    SHex_LVNTR_hexv1 = 0
                    SHex_LVNTR_hexv2 = 0
                    SHex_LVNTR_hexv3 = 0
                    SHex_LVNTR_hexv4 = 0
                    SHex_LVNTR_hexv5 = 0
                    SHex_LVNTR_hexv6 = 0

                    LHex_SVNTR_hexv0 = 0
                    LHex_SVNTR_hexv1 = 0
                    LHex_SVNTR_hexv2 = 0
                    LHex_SVNTR_hexv3 = 0
                    LHex_SVNTR_hexv4 = 0
                    LHex_SVNTR_hexv5 = 0
                    LHex_SVNTR_hexv6 = 0
                    LHex_LVNTR_hexv0 = 0
                    LHex_LVNTR_hexv1 = 0
                    LHex_LVNTR_hexv2 = 0
                    LHex_LVNTR_hexv3 = 0
                    LHex_LVNTR_hexv4 = 0
                    LHex_LVNTR_hexv5 = 0
                    LHex_LVNTR_hexv6 = 0

                    total_reads = 0
                    
                    while True:
                        header=in_fastq.readline().strip()
                        if header == "":              #break if detects blank new line
                            break 
                        seq=in_fastq.readline().strip()
                        plus=in_fastq.readline().strip()
                        qual=in_fastq.readline().strip()

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

                        #ReGex for hexamer variants that are REVERSE COMPLEMENT ###
                        rev_hex_var0_match = re.search("TCTCCCTCTCCCTCTCCC", seq)
                        rev_hex_var1_match = re.search("GGTACAGAGCCCTTTTAGAGAGGCTCCAAAGCTC[CT]+", seq)
                        rev_hex_var2_match = re.search("GGTACAGAGCCCTTTTAGAGAGGCTCCAAAGCTC[CT]+(TCTCCC){3,}[CT]+GGTACAGAGCCCTTTTAGAGAGGCTCCAAAGCTC[TC]+", seq)
                        rev_hex_var3_match = re.search("GGTACAGAGCCCTTTTAGAGAGGCTCCAAAGCTC[CT]+(TC{5}){2,}[CT]+GGTACAGAGCCCTTTTAGAGAGGCTCCAAAGCTC[TC]+", seq) #use se
                        rev_hex_var4_match = re.search("[TC]+TCTCCGTCTCCG[TC]+", seq)
                        rev_hex_var5_match = re.search("GGTACAGAGCCCTTTTAGAGAGGCTCCAAAGCTC[CT]+[TCTCCCCTCTCCCC]", seq)
                        rev_hex_var6_match = re.search("TCTCCCCTCTCCCCTCTCCCC", seq)

                        # fuzzy string matches of the different cases:
                        hex_var4 = "GGGAGACGGAGACGGAGAGGGAGA"
                        find_hexvar4 = find_near_matches(hex_var4, seq, max_l_dist = 1)
                        find_hexvar4_start = [match.start for match in find_hexvar4]
                        
                        
                        if match:
                            hexamer = match.group(2)
                            VNTR = match.group(1)

                            # print(f'hexamer length per read:{len(hexamer)}')

                            #long hexamer
                            if len(hexamer) >= 380:
                                long_hex += 1
                                if hex_var2_match:
                                    LHex_hexv2+=1
                                elif hex_var3_match:
                                    LHex_hexv3+=1
                                elif hex_var4_match:
                                    LHex_hexv4+=1
                                elif hex_var5_match:
                                    LHex_hexv5+=1
                                elif hex_var6_match:
                                    Lhex_hexv6+=1
                                elif hex_var1_match:
                                    LHex_hexv1+=1
                                elif hex_var0_match:
                                    LHex_hexv0+=1
                                # print(f'This is the hexamer length: {len(hexamer)}')
                                
                                if len(VNTR) < 2500: #short vntr
                                    L_short_VNTR +=1
                                    short_VNTR+=1
                                    
                                    if hex_var2_match:
                                        LHex_SVNTR_hexv2+=1
                                    elif hex_var3_match:
                                        LHex_SVNTR_hexv3+=1
                                    elif hex_var4_match:
                                        LHex_SVNTR_hexv4+=1
                                    elif hex_var5_match:
                                        LHex_SVNTR_hexv5+=1
                                    elif hex_var6_match:
                                        LHex_SVNTR_hexv6+=1
                                    elif hex_var1_match:
                                        LHex_SVNTR_hexv1+=1
                                    elif hex_var0_match:
                                        LHex_SVNTR_hexv0+=1

                                if len(VNTR) >= 2500: #long vntr
                                    L_long_VNTR +=1
                                    long_VNTR+=1
                                    
                                    if hex_var2_match:
                                        LHex_LVNTR_hexv2+=1
                                    elif hex_var3_match:
                                        LHex_LVNTR_hexv3+=1
                                    elif hex_var4_match:
                                        LHex_LVNTR_hexv4+=1
                                    elif hex_var5_match:
                                        LHex_LVNTR_hexv5+=1
                                    elif hex_var6_match:
                                        LHex_LVNTR_hexv6+=1
                                    elif hex_var1_match:
                                        LHex_LVNTR_hexv1+=1
                                    elif hex_var0_match:
                                        LHex_LVNTR_hexv0+=1


                        elif match2: #reverse complement:
                            hexamer = match2.group(1)
                            VNTR = match2.group(2)
                            # print("match 2 : rev comp found")
                            # print(f'seq: {seq}')
          
                            #long hexamer
                            # print(f'Length of the hexamer: {len(hexamer)}')
                            if len(hexamer) >= 380:
                                long_hex += 1
                                if rev_hex_var2_match:
                                    LHex_hexv2+=1
                                elif rev_hex_var3_match:
                                    LHex_hexv3+=1
                                elif rev_hex_var4_match:
                                    LHex_hexv4+=1
                                elif rev_hex_var5_match:
                                    LHex_hexv5+=1
                                elif rev_hex_var6_match:
                                    Lhex_hexv6+=1
                                elif rev_hex_var1_match:
                                    LHex_hexv1+=1
                                elif rev_hex_var0_match:
                                    LHex_hexv0+=1
                                # print("shortvntr")
                                if len(VNTR) < 2500:
                                    L_short_VNTR +=1
                                    short_VNTR+=1
                                    # print("shortvntr")
                                    
                                    if rev_hex_var2_match:
                                        LHex_SVNTR_hexv2+=1
                                        # print("c2vs")
                                    elif rev_hex_var3_match:
                                        LHex_SVNTR_hexv3+=1
                                        # print("c3vs")
                                    elif rev_hex_var4_match:
                                        LHex_SVNTR_hexv4+=1
                                        # print("v4vs")
                                    elif rev_hex_var5_match:
                                        LHex_SVNTR_hexv5+=1
                                    elif rev_hex_var6_match:
                                        LHex_SVNTR_hexv6+=1
                                    elif rev_hex_var1_match:
                                        LHex_SVNTR_hexv1+=1
                                        # print("case1shortvntr")
                                    elif rev_hex_var0_match:
                                        LHex_SVNTR_hexv0+=1

                                if len(VNTR) >= 2500:
                                    L_long_VNTR +=1
                                    long_VNTR+=1
                                    # print("longvntr")
                                    
                                    if rev_hex_var2_match:
                                        LHex_LVNTR_hexv2+=1
                                        # print("v2vL")
                                    elif rev_hex_var3_match:
                                        LHex_LVNTR_hexv3+=1
                                        # print("v3vL")
                                    elif rev_hex_var4_match:
                                        LHex_LVNTR_hexv4+=1
                                        # print("v4vL")
                                    elif rev_hex_var5_match:
                                        LHex_LVNTR_hexv5+=1
                                    elif rev_hex_var6_match:
                                        LHex_LVNTR_hexv6+=1
                                    elif rev_hex_var1_match:
                                        LHex_LVNTR_hexv1+=1
                                        # print("case1longvntr")
                                    elif rev_hex_var0_match:
                                        LHex_LVNTR_hexv0+=1
                                    
                            


                        else: # use approximate string matching to determine lengths
                            primer1seq = "ACGAAAACCAGTGAGGCGTAG"
                            primer1 = find_near_matches(primer1seq,seq, max_l_dist=2)
                            # print(primer1)
                            start_position1 = [match.start for match in primer1]
                            if start_position1: #check if the position exists
                                start_position1 = start_position1[0]
                                VNTR = start_position1
                            
                            primer2seq = "CTCTGTGCAATCGGAGTAGAGG"
                            primer2 = find_near_matches(primer2seq,seq, max_l_dist=2)
                            # print(primer2)
                            start_position2 = [match.end for match in primer2]
                            if start_position2: #check if the position exists
                                start_position2 = start_position2[0]
                            # print(f'{start_position2=}')

                            hex_var1 = ""
                            hex_var2 = ""
                            hex_var3 = "GGGAGAGCTTTGGAGCCTCTCTAAAAGGGCTCTGTACCGGAGAGGGGAGAGGGGGAGGGGGAGGGGAGGGGGAGGGGGAGGGGGAGCTTTGGAGCCTCTCTAAAAGGGCTCTGTACC"
                            
                            hex_var4 = "GGGAGACGGAGACGGAGAGGGAGA"
                            find_hexvar4 = find_near_matches(hex_var4, seq, max_l_dist = 1)
                            find_hexvar4_start = [match.start for match in find_hexvar4]

                            if start_position1 and start_position2:
                                # print(f'{start_position2 - start_position1}')
                            # print(f'Length between primers:{start_position2[0]-start_position1[0]}')
                                
                                # long hexamer
                                if (start_position2 - start_position1) >= 380:
                                    long_hex += 1
                                    if hex_var2_match:
                                        LHex_hexv2+=1
                                    elif hex_var3_match:
                                        LHex_hexv3+=1
                                    elif hex_var4_match:
                                        LHex_hexv4+=1
                                    elif hex_var5_match:
                                        LHex_hexv5+=1
                                    elif hex_var6_match:
                                        Lhex_hexv6+=1
                                    elif hex_var1_match:
                                        LHex_hexv1+=1
                                    elif hex_var0_match:
                                        LHex_hexv0+=1

                                    if (VNTR) < 2500:
                                        L_short_VNTR +=1
                                        short_VNTR+=1
                                        if hex_var2_match:
                                            LHex_SVNTR_hexv2+=1
                                        elif hex_var3_match:
                                            LHex_SVNTR_hexv3+=1
                                        elif hex_var4_match:
                                            LHex_SVNTR_hexv4+=1
                                        elif hex_var5_match:
                                            LHex_SVNTR_hexv5+=1
                                        elif hex_var6_match:
                                            LHex_SVNTR_hexv6+=1
                                        elif hex_var1_match:
                                            LHex_SVNTR_hexv1+=1
                                        elif hex_var0_match:
                                            LHex_SVNTR_hexv0+=1

                                    if (VNTR) >= 2500:
                                        L_long_VNTR +=1 
                                        long_VNTR+=1
                                        
                                        if hex_var2_match:
                                            LHex_LVNTR_hexv2+=1
                                        elif hex_var3_match:
                                            LHex_LVNTR_hexv3+=1
                                        elif hex_var4_match:
                                            LHex_LVNTR_hexv4+=1
                                        elif hex_var5_match:
                                            LHex_LVNTR_hexv5+=1
                                        elif hex_var6_match:
                                            LHex_LVNTR_hexv6+=1
                                        elif hex_var1_match:
                                            LHex_LVNTR_hexv1+=1
                                        elif hex_var0_match:
                                            LHex_LVNTR_hexv0+=1
                                        
                        total_reads += 1

###############################################################################

            elif R1 < 29 and R2 < 29: ### in -/- cases
                print("-/-", file)
                print(f'Calling -/- allele for: {file}')
                
                with open(file, "r") as in_fastq:
                    #hexamer counts
                    short_hex = 0
                    long_hex = 0
                    short_VNTR = 0
                    long_VNTR = 0

                    #hexamer-VNTR counts
                    S_short_VNTR = 0
                    S_long_VNTR = 0
                    L_short_VNTR = 0
                    L_long_VNTR = 0

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

                    #hexamer-VNTR-hexvar counts
                    SHex_SVNTR_hexv0 = 0
                    SHex_SVNTR_hexv1 = 0
                    SHex_SVNTR_hexv2 = 0
                    SHex_SVNTR_hexv3 = 0
                    SHex_SVNTR_hexv4 = 0
                    SHex_SVNTR_hexv5 = 0
                    SHex_SVNTR_hexv6 = 0
                    SHex_LVNTR_hexv0 = 0
                    SHex_LVNTR_hexv1 = 0
                    SHex_LVNTR_hexv2 = 0
                    SHex_LVNTR_hexv3 = 0
                    SHex_LVNTR_hexv4 = 0
                    SHex_LVNTR_hexv5 = 0
                    SHex_LVNTR_hexv6 = 0

                    LHex_SVNTR_hexv0 = 0
                    LHex_SVNTR_hexv1 = 0
                    LHex_SVNTR_hexv2 = 0
                    LHex_SVNTR_hexv3 = 0
                    LHex_SVNTR_hexv4 = 0
                    LHex_SVNTR_hexv5 = 0
                    LHex_SVNTR_hexv6 = 0
                    LHex_LVNTR_hexv0 = 0
                    LHex_LVNTR_hexv1 = 0
                    LHex_LVNTR_hexv2 = 0
                    LHex_LVNTR_hexv3 = 0
                    LHex_LVNTR_hexv4 = 0
                    LHex_LVNTR_hexv5 = 0
                    LHex_LVNTR_hexv6 = 0

                    total_reads = 0
                    
                    while True:
                        header=in_fastq.readline().strip()
                        if header == "":              #break if detects blank new line
                            break 
                        seq=in_fastq.readline().strip()
                        plus=in_fastq.readline().strip()
                        qual=in_fastq.readline().strip()

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

                        #ReGex for hexamer variants that are REVERSE COMPLEMENT ###
                        rev_hex_var0_match = re.search("TCTCCCTCTCCCTCTCCC", seq)
                        rev_hex_var1_match = re.search("GGTACAGAGCCCTTTTAGAGAGGCTCCAAAGCTC[CT]+", seq)
                        rev_hex_var2_match = re.search("GGTACAGAGCCCTTTTAGAGAGGCTCCAAAGCTC[CT]+(TCTCCC){3,}[CT]+GGTACAGAGCCCTTTTAGAGAGGCTCCAAAGCTC[TC]+", seq)
                        rev_hex_var3_match = re.search("GGTACAGAGCCCTTTTAGAGAGGCTCCAAAGCTC[CT]+(TC{5}){2,}[CT]+GGTACAGAGCCCTTTTAGAGAGGCTCCAAAGCTC[TC]+", seq) #use se
                        rev_hex_var4_match = re.search("[TC]+TCTCCGTCTCCG[TC]+", seq)
                        rev_hex_var5_match = re.search("GGTACAGAGCCCTTTTAGAGAGGCTCCAAAGCTC[CT]+[TCTCCCCTCTCCCC]", seq)
                        rev_hex_var6_match = re.search("TCTCCCCTCTCCCCTCTCCCC", seq)


                        # fuzzy string matches of the different cases:
                        hex_var4 = "GGGAGACGGAGACGGAGAGGGAGA"
                        find_hexvar4 = find_near_matches(hex_var4, seq, max_l_dist = 1)
                        find_hexvar4_start = [match.start for match in find_hexvar4]
                        
                        
                        if match:
                            hexamer = match.group(2)
                            VNTR = match.group(1)
                            
                    
                            #short hexamer
                            if len(hexamer) < 403:
                                short_hex += 1
                                if hex_var2_match:
                                    SHex_hexv2+=1
                                elif hex_var3_match:
                                    SHex_hexv3+=1
                                elif hex_var4_match:
                                    SHex_hexv4+=1
                                elif hex_var5_match:
                                    SHex_hexv5+=1
                                elif hex_var6_match:
                                    Shex_hexv6+=1
                                elif hex_var1_match:
                                    SHex_hexv1+=1
                                elif hex_var0_match:
                                    SHex_hexv0+=1

                                if len(VNTR) < 2500: #short vntr
                                        S_short_VNTR +=1
                                        short_VNTR+=1
                                        
                                        if hex_var2_match:
                                            SHex_SVNTR_hexv2+=1
                                        elif hex_var3_match:
                                            SHex_SVNTR_hexv3+=1
                                        elif hex_var4_match:
                                            SHex_SVNTR_hexv4+=1
                                        elif hex_var5_match:
                                            SHex_SVNTR_hexv5+=1
                                        elif hex_var6_match:
                                            SHex_SVNTR_hexv6+=1
                                        elif hex_var1_match:
                                            SHex_SVNTR_hexv1+=1
                                        elif hex_var0_match:
                                            SHex_SVNTR_hexv0+=1

                                if len(VNTR) >= 2500: #long vntr
                                        S_long_VNTR +=1
                                        long_VNTR+=1   
                                        
                                        if hex_var2_match:
                                            SHex_LVNTR_hexv2+=1
                                        elif hex_var3_match:
                                            SHex_LVNTR_hexv3+=1
                                        elif hex_var4_match:
                                            SHex_LVNTR_hexv4+=1  
                                        elif hex_var5_match:
                                            SHex_LVNTR_hexv5+=1
                                        elif hex_var6_match:
                                            SHex_LVNTR_hexv6+=1
                                        elif hex_var1_match:
                                            SHex_LVNTR_hexv1+=1 
                                        elif hex_var0_match:
                                            SHex_LVNTR_hexv0+=1
                            

                        elif match2: #reverse complement:
                            hexamer = match2.group(1)
                            VNTR = match2.group(2)
                            #short hexamer
                            if len(hexamer) < 403:
                                short_hex += 1
                                if rev_hex_var2_match:
                                    SHex_hexv2+=1
                                elif rev_hex_var3_match:
                                    SHex_hexv3+=1
                                elif rev_hex_var4_match:
                                    SHex_hexv4+=1
                                elif rev_hex_var5_match:
                                    SHex_hexv5+=1
                                elif rev_hex_var6_match:
                                    Shex_hexv6+=1
                                elif rev_hex_var1_match:
                                    SHex_hexv1+=1
                                elif rev_hex_var0_match:
                                    SHex_hexv0+=1
                                if len(VNTR) < 2500:
                                        S_short_VNTR +=1
                                        short_VNTR+=1
                                        
                                        if rev_hex_var2_match:
                                            SHex_SVNTR_hexv2+=1
                                        elif rev_hex_var3_match:
                                            SHex_SVNTR_hexv3+=1
                                        elif rev_hex_var4_match:
                                            SHex_SVNTR_hexv4+=1
                                        elif rev_hex_var5_match:
                                            SHex_SVNTR_hexv5+=1
                                        elif rev_hex_var6_match:
                                            SHex_SVNTR_hexv6+=1
                                        elif rev_hex_var1_match:
                                            SHex_SVNTR_hexv1+=1
                                        elif rev_hex_var0_match:
                                            SHex_SVNTR_hexv0+=1

                                if len(VNTR) >= 2500:
                                        S_long_VNTR +=1
                                        long_VNTR+=1
                                        
                                        if rev_hex_var2_match:
                                            SHex_LVNTR_hexv2+=1
                                        elif rev_hex_var3_match:
                                            SHex_LVNTR_hexv3+=1
                                        elif rev_hex_var4_match:
                                            SHex_LVNTR_hexv4+=1
                                        elif rev_hex_var5_match:
                                            SHex_LVNTR_hexv5+=1
                                        elif rev_hex_var6_match:
                                            SHex_LVNTR_hexv6+=1
                                        elif rev_hex_var1_match:
                                            SHex_LVNTR_hexv1+=1
                                        elif rev_hex_var0_match:
                                            SHex_LVNTR_hexv0+=1
         
                            
                        else: # use approximate string matching to determine lengths
                            primer1seq = "ACGAAAACCAGTGAGGCGTAG"
                            primer1 = find_near_matches(primer1seq,seq, max_l_dist=2)
                            # print(primer1)
                            start_position1 = [match.start for match in primer1]
                            if start_position1: #check if the position exists
                                start_position1 = start_position1[0]
                                VNTR = start_position1
                            
                            primer2seq = "CTCTGTGCAATCGGAGTAGAGG"
                            primer2 = find_near_matches(primer2seq,seq, max_l_dist=2)
                            # print(primer2)
                            start_position2 = [match.end for match in primer2]
                            if start_position2: #check if the position exists
                                start_position2 = start_position2[0]
                            # print(f'{start_position2=}')

                            hex_var1 = ""
                            hex_var2 = ""
                            hex_var3 = "GGGAGAGCTTTGGAGCCTCTCTAAAAGGGCTCTGTACCGGAGAGGGGAGAGGGGGAGGGGGAGGGGAGGGGGAGGGGGAGGGGGAGCTTTGGAGCCTCTCTAAAAGGGCTCTGTACC"
                            
                            hex_var4 = "GGGAGACGGAGACGGAGAGGGAGA"
                            find_hexvar4 = find_near_matches(hex_var4, seq, max_l_dist = 1)
                            find_hexvar4_start = [match.start for match in find_hexvar4]

                            if start_position1 and start_position2:
                                # print(f'{start_position2 - start_position1}')
                            # print(f'Length between primers:{start_position2[0]-start_position1[0]}')
                                if (start_position2 - start_position1) < 403:
                                    short_hex += 1
                                    if hex_var2_match:
                                        SHex_hexv2+=1
                                    elif hex_var3_match:
                                        SHex_hexv3+=1
                                    elif hex_var4_match:
                                        SHex_hexv4+=1
                                    elif hex_var5_match:
                                        SHex_hexv5+=1
                                    elif hex_var6_match:
                                        Shex_hexv6+=1
                                    elif hex_var1_match:
                                        SHex_hexv1+=1
                                    elif hex_var0_match:
                                        SHex_hexv0+=1

                                    if (VNTR) < 2500:
                                        S_short_VNTR +=1
                                        short_VNTR+=1
                                        
                                        if hex_var2_match:
                                            SHex_SVNTR_hexv2+=1
                                        elif hex_var3_match:
                                            SHex_SVNTR_hexv3+=1
                                        elif hex_var4_match:
                                            SHex_SVNTR_hexv4+=1
                                        elif hex_var5_match:
                                            SHex_SVNTR_hexv5+=1
                                        elif hex_var6_match:
                                            SHex_SVNTR_hexv6+=1
                                        elif hex_var1_match:
                                            SHex_SVNTR_hexv1+=1
                                        elif hex_var0_match:
                                            SHex_SVNTR_hexv0+=1

                                    if (VNTR) >= 2500:
                                        S_long_VNTR +=1
                                        long_VNTR+=1
                                        
                                        if hex_var2_match:
                                            SHex_LVNTR_hexv2+=1
                                        elif hex_var3_match:
                                            SHex_LVNTR_hexv3+=1
                                        elif hex_var4_match:
                                            SHex_LVNTR_hexv4+=1
                                        elif hex_var5_match:
                                            SHex_LVNTR_hexv5+=1
                                        elif hex_var6_match:
                                            SHex_LVNTR_hexv6+=1
                                        elif hex_var1_match:
                                            SHex_LVNTR_hexv1+=1
                                        elif hex_var0_match:
                                            SHex_LVNTR_hexv0+=1
                                        
                        total_reads += 1
                    
        print(f'{total_reads=}')
        
############################################################################################################        

        # Determine the most abundant VNTR type for short hexamers
        most_abundant_shortHex_VNTR = 'short_VNTR' if S_short_VNTR >= S_long_VNTR else 'long_VNTR'
        most_abundant_shortHex_VNTR_count = max(S_short_VNTR, S_long_VNTR)
        second_abundant_short_VNTR = 'short_VNTR' if S_short_VNTR < S_long_VNTR else 'long_VNTR'
        second_abundant_short_VNTR_count = min(S_short_VNTR, S_long_VNTR)
        
        # Determine the most abundant VNTR type for long hexamers
        most_abundant_longHex_VNTR = 'short_VNTR' if L_short_VNTR >= L_long_VNTR else 'long_VNTR'
        most_abundant_longHex_VNTR_count = max(L_short_VNTR, L_long_VNTR)
        second_abundant_long_VNTR = 'long_VNTR' if L_short_VNTR < L_long_VNTR else 'long_VNTR'
        second_abundant_long_VNTR_count = min(S_short_VNTR, S_long_VNTR)

        print(f'{most_abundant_shortHex_VNTR_count=}')
        print(f'{most_abundant_longHex_VNTR_count=}')


        # Determine most abundant hex-hexvar for short hexamers:
        shorthex_hexvar_counts = {
            'SHex_hexv0': SHex_hexv0,
            'SHex_hexv1': SHex_hexv1,
            'SHex_hexv2': SHex_hexv2,
            'SHex_hexv3': SHex_hexv3,
            'SHex_hexv4': SHex_hexv4,
            'SHex_hexv5': SHex_hexv5,
            'SHex_hexv6': Shex_hexv6
        }
        sorted_shorthexhexvar_counts = sorted(shorthex_hexvar_counts.items(), key=lambda item: item[1], reverse=True)
        most_abundantshorthexhexvar = sorted_shorthexhexvar_counts[0]
        most_abundant_typeshorthexhexvar, most_abundant_countshorthexhexvar = most_abundantshorthexhexvar
        print(f"Most abundant short hex vntr full: {most_abundant_typeshorthexhexvar} with count: {most_abundant_countshorthexhexvar}")
        second_abundantshorthexhexvar = sorted_shorthexhexvar_counts[1]
        second_abundant_typeshorthexhexvar, second_abundant_countshorthexhexvar = second_abundantshorthexhexvar
        print(f"SECOND most abundant short hex vntr full: {second_abundant_typeshorthexhexvar} with count: {second_abundant_countshorthexhexvar}")

        # Determine most abundant hex-hexvar for long hexamers:
        longhex_hexvar_counts = {
            'LHex_hexv0': LHex_hexv0,
            'LHex_hexv1': LHex_hexv1,
            'LHex_hexv2': LHex_hexv2,
            'LHex_hexv3': LHex_hexv3,
            'LHex_hexv4': LHex_hexv4,
            'LHex_hexv5': LHex_hexv5,
            'LHex_hexv6': Lhex_hexv6
        }
        sorted_longhexhexvar_counts = sorted(longhex_hexvar_counts.items(), key=lambda item: item[1], reverse=True)
        most_abundantlonghexhexvar = sorted_longhexhexvar_counts[0]
        most_abundant_typelonghexhexvar, most_abundant_countlonghexhexvar = most_abundantlonghexhexvar
        print(f"Most abundant long hex vntr full: {most_abundant_typelonghexhexvar} with count: {most_abundant_countlonghexhexvar}")
        second_abundantlonghexhexvar = sorted_longhexhexvar_counts[1]
        second_abundant_typelonghexhexvar, second_abundant_countlonghexhexvar = second_abundantlonghexhexvar
        print(f"SECOND most abundant long hex vntr full: {second_abundant_typelonghexhexvar} with count: {second_abundant_countlonghexhexvar}")


        # Determine most abundant hex-VNTR-hexvar for short hexamers:
        shorthexVNTRhexvar_counts = {
            'SHex_SVNTR_hexv0': SHex_LVNTR_hexv0,
            'SHex_SVNTR_hexv1': SHex_SVNTR_hexv1,
            'SHex_SVNTR_hexv2': SHex_SVNTR_hexv2,
            'SHex_SVNTR_hexv3': SHex_SVNTR_hexv3,
            'SHex_SVNTR_hexv4': SHex_SVNTR_hexv4,
            'SHex_SVNTR_hexv5': SHex_SVNTR_hexv5,
            'SHex_SVNTR_hexv6': SHex_SVNTR_hexv6,
            'SHex_LVNTR_hexv0': SHex_LVNTR_hexv0,
            'SHex_LVNTR_hexv1': SHex_LVNTR_hexv1,
            'SHex_LVNTR_hexv2': SHex_LVNTR_hexv2,
            'SHex_LVNTR_hexv3': SHex_LVNTR_hexv3,
            'SHex_LVNTR_hexv4': SHex_LVNTR_hexv4,
            'SHex_LVNTR_hexv5': SHex_LVNTR_hexv5,
            'SHex_LVNTR_hexv6': SHex_LVNTR_hexv6
        } 
        sorted_shorthexVNTRhexvar_counts = sorted(shorthexVNTRhexvar_counts.items(), key=lambda item: item[1], reverse=True)
        most_abundantshorthex = sorted_shorthexVNTRhexvar_counts[0]
        most_abundant_typeshort, most_abundant_countshort = most_abundantshorthex
        print(f"Most abundant short hex vntr full: {most_abundant_typeshort} with count: {most_abundant_countshort}")
        second_abundantshorthex = sorted_shorthexVNTRhexvar_counts[1]
        second_abundant_typeshort, second_abundant_countshort = second_abundantshorthex
        print(f"SECOND most abundant short hex vntr full: {second_abundant_typeshort} with count: {second_abundant_countshort}")

        # Determine most abundant hex-VNTR-hexvar for long hexamers:
        longhexVNTRhexvar_counts = {
            'LHex_SVNTR_hexv0': LHex_SVNTR_hexv0,
            'LHex_SVNTR_hexv1': LHex_SVNTR_hexv1,
            'LHex_SVNTR_hexv2': LHex_SVNTR_hexv2,
            'LHex_SVNTR_hexv3': LHex_SVNTR_hexv3,
            'LHex_SVNTR_hexv4': LHex_SVNTR_hexv4,
            'LHex_SVNTR_hexv5': LHex_SVNTR_hexv5,
            'LHex_SVNTR_hexv6': LHex_SVNTR_hexv6,
            'LHex_LVNTR_hexv0': LHex_LVNTR_hexv0,
            'LHex_LVNTR_hexv1': LHex_LVNTR_hexv1,
            'LHex_LVNTR_hexv2': LHex_LVNTR_hexv2,
            'LHex_LVNTR_hexv3': LHex_LVNTR_hexv3,
            'LHex_LVNTR_hexv4': LHex_LVNTR_hexv4,
            'LHex_LVNTR_hexv5': LHex_LVNTR_hexv5,
            'LHex_LVNTR_hexv6': LHex_LVNTR_hexv6
        } 
        sorted_longhexVNTRhexvar_counts = sorted(longhexVNTRhexvar_counts.items(), key=lambda item: item[1], reverse=True)
        most_abundantlonghex = sorted_longhexVNTRhexvar_counts[0]
        most_abundant_typelong, most_abundant_countlong = most_abundantlonghex
        print(f"Most abundant long hex vntr full: {most_abundant_typelong} with count: {most_abundant_countlong}")
        second_abundantlonghex = sorted_longhexVNTRhexvar_counts[1]
        second_abundant_typelong, second_abundant_countlong = second_abundantlonghex
        print(f"SECOND most abundant LONG hex vntrfull: {second_abundant_typelong} with count: {second_abundant_countlong}")
        

        ####
        # Calling allele
        ####  
        with open(args.out_table, "a") as out_table:
               
            if filename_match:
                R1 = int(filename_match.group(1))
                R2 = int(filename_match.group(2))

            
            #+/-
            print(f'{R1}----------{R2}')
            if R1 < 29 and R2 >= 29:

                alleles = [('ShortHex', short_hex, (short_hex/total_reads), most_abundant_shortHex_VNTR, most_abundant_shortHex_VNTR_count, most_abundant_typeshort, most_abundant_countshort, second_abundant_typeshort, second_abundant_countshort, most_abundant_typeshorthexhexvar, most_abundant_countshorthexhexvar, total_reads),
                        ('LongHex',long_hex, (long_hex/total_reads), most_abundant_longHex_VNTR, most_abundant_longHex_VNTR_count, most_abundant_typelonghexhexvar, most_abundant_countlonghexhexvar, most_abundant_typelong, most_abundant_countlong, second_abundant_typelong, second_abundant_countlong ) ]
                top_two_alleles = alleles[:2]
                                
                out_table.write(f'{file}\t{top_two_alleles[0][11]}\t{top_two_alleles[0][0]}\t{top_two_alleles[0][1]}\t{top_two_alleles[0][2]}\t{top_two_alleles[0][3]}\t{top_two_alleles[0][4]}\t{top_two_alleles[0][5]}\t{top_two_alleles[0][6]}\t{top_two_alleles[0][9]}\t{top_two_alleles[0][10]}\t{top_two_alleles[1][0]}\t{top_two_alleles[1][1]}\t{top_two_alleles[1][2]}\t{top_two_alleles[1][3]}\t{top_two_alleles[1][4]}\t{top_two_alleles[1][5]}\t{top_two_alleles[1][6]}\t{top_two_alleles[1][7]}\t{top_two_alleles[1][8]}\t{top_two_alleles[0][7]}\t{top_two_alleles[0][8]}\t{top_two_alleles[1][9]}\t{top_two_alleles[1][10]}\n')

            #+/+
            elif R1 >= 29 and R2 >= 29:
                alleles = [('LongHex',long_hex, (long_hex/total_reads), most_abundant_longHex_VNTR, most_abundant_longHex_VNTR_count, most_abundant_typelong, most_abundant_countlong, most_abundant_typelonghexhexvar, most_abundant_countlonghexhexvar,total_reads),
                        ('LongHex',long_hex, (long_hex/total_reads), most_abundant_longHex_VNTR, most_abundant_longHex_VNTR_count, most_abundant_typelonghexhexvar, most_abundant_countlonghexhexvar, second_abundant_typelong, second_abundant_countlong) ]
                top_two_alleles = alleles[:2]

                out_table.write(f'{file}\t{top_two_alleles[0][9]}\t{top_two_alleles[0][0]}\t{top_two_alleles[0][1]}\t{top_two_alleles[0][2]}\t{top_two_alleles[0][3]}\t{top_two_alleles[0][4]}\t{top_two_alleles[0][5]}\t{top_two_alleles[0][6]}\t{top_two_alleles[0][7]}\t{top_two_alleles[0][8]}\t{top_two_alleles[1][0]}\t{top_two_alleles[1][1]}\t{top_two_alleles[1][2]}\t{top_two_alleles[1][3]}\t{top_two_alleles[1][4]}\t{top_two_alleles[1][5]}\t{top_two_alleles[1][6]}\t{top_two_alleles[1][7]}\t{top_two_alleles[1][8]}\tNA\tNA\tNA\tNA\n')

            #-/-
            elif R1 < 29 and R2 < 29:
                alleles = [('ShortHex', short_hex, (short_hex/total_reads), most_abundant_shortHex_VNTR, most_abundant_shortHex_VNTR_count, most_abundant_typeshort, most_abundant_countshort, most_abundant_typeshorthexhexvar, most_abundant_countshorthexhexvar, total_reads),
                        ('ShortHex', short_hex, (short_hex/total_reads), most_abundant_shortHex_VNTR, most_abundant_shortHex_VNTR_count, most_abundant_typeshorthexhexvar, most_abundant_countshorthexhexvar, second_abundant_typeshort, second_abundant_countlong ) ]
                top_two_alleles = alleles[:2]

                out_table.write(f'{file}\t{top_two_alleles[0][9]}\t{top_two_alleles[0][0]}\t{top_two_alleles[0][1]}\t{top_two_alleles[0][2]}\t{top_two_alleles[0][3]}\t{top_two_alleles[0][4]}\t{top_two_alleles[0][5]}\t{top_two_alleles[0][6]}\t{top_two_alleles[0][7]}\t{top_two_alleles[0][8]}\t{top_two_alleles[1][0]}\t{top_two_alleles[1][1]}\t{top_two_alleles[1][2]}\t{top_two_alleles[1][3]}\t{top_two_alleles[1][4]}\t{top_two_alleles[1][5]}\t{top_two_alleles[1][6]}\t{top_two_alleles[1][7]}\t{top_two_alleles[1][8]}\tNA\tNA\tNA\tNA\n')


# ./allele_call.py -f /blue/lien.nguyen/rossellwood/repeat_expansion_tools/sort_expansion_type/filtered_full_alleles_FASTQ/pool-1_1054.. -o alleleCall.txt