#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import colorama
import time
import sys
import gzip
import regex
import os



# 
# ----------------------------------------------------------------------------
# Function: read a RepeatMasker output file 
# ----------------------------------------------------------------------------
# 

def read_repeatmasker_file(rm_file, alu_type):
   
    #dict: keys=chromosomes  values= a list of Alus and associated info per chromosome
    data = {}
    
    # determine how to open the file based on its extension
    open_func = gzip.open if rm_file.suffix == ".gz" else open 
    
    #assign required data to variables
    with  open_func(rm_file, "rt") as fhandle:
        for line in fhandle:
            line = line.strip()
            line = regex.split(r"\s+", line)
            if line[0].isdigit():
                chrom = line[4]
                start = int(line[5]) 
                end = int(line[6])
                strand = line[8].replace('C', '-') #replacing minus strand with standard notation
                te_type = line[9]
                
                # verify if the TE category is a potential active Alu and also a FL element
                if alu_type in te_type:
                    if (strand == "+" and int(line[11]) <= 4 and int(line[12]) >= 267):                             
                        te_id = data.get(chrom, [])
                        te_id.append((chrom, start, end, strand, te_type))
                        data[chrom] = te_id
                    if (strand == "-" and int(line[12]) >= 267 and int(line[13]) <= 4):
                        te_id = data.get(chrom, [])
                        te_id.append((chrom, start, end, strand, te_type))
                        data[chrom] = te_id
    yield data

    
# 
# ----------------------------------------------------------------------------
# Function: read a repeat gff file 
# ----------------------------------------------------------------------------
# 

def read_gff_file(gff_file, alu_type):
    
    data = {}
    
    
    # determine how to open the file based on its extension
    open_func = gzip.open if gff_file.endswith(".gz") else open
    
    
    #assign required data to variables
    with  open_func(gff_file, "rt") as fhandle:
        for line in fhandle:
            line = line.strip().split('\t')
            chrom = line[0]
            start = int(line[3]) 
            end = int(line[4])
            strand = line[6] 
            te_type = line[8]
            
            if alu_type in te_type:
                te_inf = data.get(chrom, [])
                te_inf.append((chrom, start, end, strand, te_type))
                data[chrom] = te_inf
    yield data
   


# 
# ----------------------------------------------------------------------------
# Function: read a repeat annotation (RM, bed format)
# ----------------------------------------------------------------------------
# 

def read_bed_file(bed_file, alu_type):
    
    data = {}
    
    
    # determine how to open the file based on its extension
    open_func = gzip.open if bed_file.endswith(".gz") else open
    
    #assign required data to variables
    with  open_func(bed_file, "rt") as fhandle:
        for line in fhandle:
            line = line.strip().split('\t')
            chrom = line[0]
            start = int(line[1]) + 1
            end = int(line[2])
            strand = line[5] 
            #te_type = line[6] #[3]
            te_type = line[3] if len(line) <= 6 else line[6]
            
            if alu_type in te_type:
                te_inf = data.get(chrom, [])
                te_inf.append((chrom, start, end, strand, te_type))
                data[chrom] = te_inf
    yield data




# 
# ----------------------------------------------------------------------------
# Function: determin the file format 
# ----------------------------------------------------------------------------
# 

def determine_rm_format(file_path, user_format, te_type):
    if user_format == "gff":
        return read_gff_file(file_path, te_type)
    elif user_format == "bed":
        return read_bed_file(file_path, te_type)
    elif user_format == "RM.out":
        return read_repeatmasker_file(file_path, te_type)
    else:
        raise ValueError(f"Unsupported file format: {user_format}")


# 
# -----------------------------------------------------------------------------
# Function: check if if all bases of a TSD are the same
# -----------------------------------------------------------------------------
# 

def all_bases_same(tsd):
    first_char = tsd[0]
    return all(char == first_char for char in tsd)


# 
# -----------------------------------------------------------------------------
# Function: write output of refined detected TSDs into a file
# -----------------------------------------------------------------------------
# 

def write_output_file_tsd(input_file, path, base_name, header):
    
    homo_name = path / (base_name + "_homoTSD.tsv")
    hetero_name = path / (base_name + "_heteroTSD.tsv")
    
    with input_file.open("r") as in_file, \
       hetero_name.open("w") as out_file_het, \
        homo_name.open("w") as out_file_hom:
            
        out_file_het.write(header + "\n")
        out_file_hom.write(header + "\n")
            
    
        for line in in_file:
            if not line.strip().endswith("poly(A)_tail_pos"):
                line = line.strip().split("\t")
                tsd = line[5]
                if not all_bases_same(tsd):
                    out_file_het.write(("\t".join(line)) + "\n")
                else:
                    out_file_hom.write(("\t".join(line)) + "\n")


# 
# -----------------------------------------------------------------------------
# Function: pogress bar
# -----------------------------------------------------------------------------
# 

def progress_bar(progress, total, color=colorama.Fore.YELLOW):
    """ Show progress """
    percent = 100 * (progress / float(total))
    bar = "█" * int(percent)  + "-" * (100 - int(percent))
    msg = "Analyzing TE: "
    print(color + f"\r {msg}|{bar}| {percent:.0f}%", end="\r")
    if progress == total:
        print(colorama.Fore.GREEN + f"\r {msg}|{bar}| {percent:.0f}%", end="\r")


# 
# -----------------------------------------------------------------------------
# Function: spinning bar
# -----------------------------------------------------------------------------
# 

def spinning_bar():
    ''' prints a spinning bar on the terminal while a program is running'''
    
    while True:
        
        for c in "\\|/-":
            time.sleep(.1)
            sys.stdout.write(f"\b{c}")
            sys.stdout.flush()
