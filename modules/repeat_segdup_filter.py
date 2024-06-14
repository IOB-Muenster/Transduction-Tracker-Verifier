#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Module: functions for excluding repeated sequences vverlapping Ssegmental duplications
"""
import subprocess
import sys
import tempfile
from pathlib import Path
import regex
import glob
import os


from modules.utilities_TransductionVerifier_V2 import *


# 
# -----------------------------------------------------------------------------
# Function: Filter repeats overlapping with SegDups
# -----------------------------------------------------------------------------
# 

def filter_SegDups(rm_file, sd_file, output_dir):

    base_name = rm_file.name.split(".")[0]
    
    
    
    #~~~ A dictionary of original RM.out to retrieve data
    original_RM_dict ={}
    segdup_free_keys = []
    
    
    open_func = lambda file_path: gzip.open(file_path, "rt") if file_path.suffix == ".gz" else file_path.open("r")
    
    if sd_file and rm_file:
                
        with open_func(rm_file) as rm_infile, \
            open_func(sd_file) as sd_infile:
            
            # generating a dict for original RM.out
            for line in rm_infile:
                line = line.strip()
                line = regex.split(r"\s+", line)
                if line[0].isdigit():
                    key = f"{line[4]}:{line[5]}-{line[6]}"
                    value = original_RM_dict.get(key, "")
                    value += "\t".join(line)
                    original_RM_dict[key] = value
                    
        
        with open_func(rm_file) as rm_infile, \
            open_func(sd_file) as sd_infile:
                   
                
            with tempfile.NamedTemporaryFile(dir=output_dir, mode="w", delete=False) as rm_f_temp, \
                tempfile.NamedTemporaryFile(dir=output_dir, mode="w", delete=False) as sd_f_temp:
                    rm_f_temp.write(process_ReapetMasker_file_to_bed(rm_infile))
                    sd_f_temp.write(process_SD_file(sd_infile))
        

        
        for i in [rm_f_temp.name, sd_f_temp.name]:
            cmd_sort = f"sort -V -k1,1 -k2,2 -o {i}.tmp.sorted.bed {i}"
            subprocess.run(cmd_sort, shell=True, text=True, universal_newlines=True)
            
        with tempfile.NamedTemporaryFile(dir=output_dir, delete=False, mode="w") as tmp_sd_free_file:
            cmd_intersect = f"bedtools intersect -C -a {rm_f_temp.name}.tmp.sorted.bed  -b {sd_f_temp.name}.tmp.sorted.bed"
            intersect_results = subprocess.run(cmd_intersect, shell=True, text=True, universal_newlines=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            
            if intersect_results.returncode == 0:
                for line in intersect_results.stdout.splitlines():
                    line = line.split("\t")
                    if line[-1] == "0":
                        line = line[:-1]
                        key = f"{line[0]}:{str(int(line[1])+1)}-{line[2]}"
                        segdup_free_keys.append(key)
                        line = "\t".join(line) #convert list to str
                        tmp_sd_free_file.write(line)
                        tmp_sd_free_file.write('\n')
    
    
    segdup_free_RM_out = output_dir / (base_name + "_NoSegDup_RM.out")

    with segdup_free_RM_out.open("w") as out_file:
        for key in segdup_free_keys:
            if key in original_RM_dict:
                out_file.write(original_RM_dict[key] + "\n")
                
    
    
    
    # ~~~ Deleting temporary files
    temp_files_pattern = output_dir / "tmp*"
    temp_files = glob.glob(str(temp_files_pattern)) # str(converts PosixPath to string)
    
    for f in temp_files:
        try:
            os.remove(f)
        except Exception as e:
            print(f"{e}")

        
