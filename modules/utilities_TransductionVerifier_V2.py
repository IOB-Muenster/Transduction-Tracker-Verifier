#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import regex




#
#------------------------------------------------------------------------------
# Function: process potential Alu transductions and convert them to BED format
#------------------------------------------------------------------------------
#

def process_transduction_file(input_file):
    
    raw_file = ""
    #~~~ file containing the BED format of potential transductions
    transd_bed_format = ""
    
    with open(input_file) as f:
        #~~~ store the header (first line)
        header = f.readline().strip()        
        
        for raw_line in f:
            raw_file += raw_line
            line = raw_line.strip().split("\t")
            strand = line[1]
            alu_type = line[2]
            
            #~~~ splitting the coordinate e.g., "chr1:12-25"
            # to remove any empty strings that may result from the split
            chrom, beg, end = tuple(item for item in regex.split(r"\W", line[3]) if item)
            
            #~~~ combining the required info per line
            #turning the coordinate into a zero-based index
            new_line = "\t".join([chrom, str(int(beg) - 1), end, alu_type, ".", strand])
            transd_bed_format += new_line + "\n"
            
   
    return header, transd_bed_format, raw_file



#
#------------------------------------------------------------------------------
# Function: read a RepeatMasker file and convert its contents to BED format.
#------------------------------------------------------------------------------
#    

def process_ReapetMasker_file_to_bed(input_file):
    RM_bed_format = ""
    
    te_to_search = ["LINE", "Retroposon/SVA", "SINE/Alu"]
    
    for line in input_file:
        line = line.strip()
        line = regex.split(r"\s+", line)
        if line[0].isdigit() and any(line[10].startswith(te) for te in te_to_search):
            chrom = line[4]
            beg = line[5]
            end = line[6]
            strand = line[8].replace("C", "-")
            repeat_type = line[9]
            new_line = "\t".join([chrom, str(int(beg) - 1), end, repeat_type, ".", strand])
            RM_bed_format += new_line + "\n"
        
    #~~~ return the bed format of RM out
    return RM_bed_format


#
#------------------------------------------------------------------------------
# Function: summarize the contents of a segdup file.
#------------------------------------------------------------------------------
#    

def process_SD_file(input_file):
    SD_bed_format = ""
    
    for line in input_file:
        if not line.startswith("#"):
            line = line.strip()
            line = regex.split(r"\s+", line)
            chrom = line[0]
            beg = str(int(line[1]))
            end = line[2]
            strand = line[5] 
            new_line = "\t".join([chrom, beg, end, ".", ".", strand])
            SD_bed_format += new_line + "\n"
            
    return SD_bed_format




#
#------------------------------------------------------------------------------
# Function: custome sort key for chromosomes
#------------------------------------------------------------------------------
# 

def custome_sort_key(line: str) -> tuple:
    chrom_order = {
        "chr1": 1,
        "chr2": 2,
        "chr2A": 3,
        "chr2B": 4,
        "chr3": 5,
        "chr4": 6,
        "chr5": 7,
        "chr6": 8,
        "chr7": 9,
        "chr8": 10,
        "chr9": 11,
        "chr10": 12,
        "chr11": 13,
        "chr12": 14,
        "chr13": 15,
        "chr14": 16,
        "chr15": 17,
        "chr16": 18,
        "chr17": 19,
        "chr18": 20,
        "chr19": 21,
        "chr20": 22,
        "chr21": 23,
        "chr22": 24,
        "chrX": 25,
        "chrY": 26,
        "chrM": 27
        }

    line = line.split("\t")
    chrom = regex.split(r"[:-]", line[0])[0]
    beg = regex.split(r"[:-]", line[0])[1]
    return chrom_order.get(chrom, 999), chrom, int(beg)