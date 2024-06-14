#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#~~~  External modules
import os
import subprocess
import gzip
import regex
import tempfile
import textwrap
import argparse
import glob 
from colorama import Fore, init, Style, Back 
import timeit
import time
from Bio.Seq import Seq
from Bio import SeqIO
from collections import namedtuple
from pathlib import Path
import ast
import pandas as pd

# Internal module(s)
from modules.utilities_TransductionVerifier_V2 import *
from modules.utilities_TransductionTracker_V2 import read_repeatmasker_file
from modules.refine_blat_output import *

#
#------------------------------------------------------------------------------
# Main code
#------------------------------------------------------------------------------
#

def main():
    
    
    # for colorama
    init(autoreset=True)
    start_time = timeit.default_timer()
    
    #~~~ Setting up the parameters
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent("""\
                                    Verifies detected Alu transductions.
                                    """),
        usage=textwrap.dedent("""python TransductionVerifier.py -i <potential_transd.tsv> -r <RepeatMasker.out> -g <reference_genome.fa> -t <Alu family> -s <path/to/pslScore.pl> -c <cpu_number> 
                              """),
        epilog=textwrap.dedent("""
                               Example: python TransductionVerifier.py -i alu_transd_with_pA.tsv -r RM.out -g ref_genome.fa -t AluY -s pslScore.pl -c 4
                               """
        )
        )
        
    parser.add_argument("-i", "--transduction_file", type=str, required=True, metavar="", 
                        help="Path to the potential transduction file (TSV format) detected in by TransductionTracker.py.")
    parser.add_argument("-r", "--repeat_annotation", required=True, type=str, metavar="",
                        help='Path to the segmental duplication-free RepeatMasker output file (.out format). Use "*_NoSegDup_RM.out" if SegDup annotation was provided earlier with "TransductionTracker.V2.py".')
    parser.add_argument("-g", "--reference_genome", required=True, type=str, metavar="",
                        help='Path to the reference genome file.')
    parser.add_argument("-t", "--alu_type", type=str, default="AluY", metavar="",
                               help="Type of Alu element to process. The default is 'AluY'.")
    parser.add_argument("-s", "--pslScore", type=str, required=True, metavar="",
                               help="Path to the pslScore.pl script downloaded from UCSC.")
    parser.add_argument("-c", "--cpus", default=4, type=str, metavar="",
                        help='Number of CPUs for the alignment. The default is 4.')
    
    args = parser.parse_args()
    
    transd_file = args.transduction_file
    RM_file = Path(args.repeat_annotation).resolve()
    ref_genome = args.reference_genome
    alu_type = args.alu_type
    num_cpus = args.cpus
    
    
    #~~~~~~~~~ Filtering out transductions overlapping with repeats
    
    print(f"{Style.BRIGHT}{Fore.YELLOW}Initializing TransductionVerifier.V2 ...")
    time.sleep(5)
    
    
    print(f"{Style.BRIGHT}{Fore.YELLOW}Filltering out transductions overlapping with other TEs ...")
    time.sleep(5)
    
    #~~~ getting the path of the transduction file
    # this path will be used to write the output
    path = os.path.dirname(os.path.abspath(transd_file))
    
    #~~~ reading and processing the transduction file
    header, bed_transduction, raw_input_transduction = process_transduction_file(transd_file)
    
    
    #~~~ reading, processing, and converting the RM file into the BED format.
    if str(RM_file).endswith("out.gz"):
        with gzip.open(RM_file,"rt") as f:
            bed_RM = process_ReapetMasker_file_to_bed(f)
    else:
        with open(RM_file) as f:
            bed_RM = process_ReapetMasker_file_to_bed(f)
            
    
    
    #~~~ create temp files to store the content two input files in the bed format
    with tempfile.NamedTemporaryFile(dir=path, mode="w", delete=False) as transd_temp_file, \
        tempfile.NamedTemporaryFile(dir=path, mode="w", delete=False) as RM_temp_file:
        transd_temp_file.write(bed_transduction)
        RM_temp_file.write(bed_RM)
    
    #~~~ command: sort: transduction file and RM
    for i in [transd_temp_file.name, RM_temp_file.name]:
        cmd = f"sort -V -k1,1 -k2,2 -o {i}.tmp.sorted.bed {i}"
        subprocess.run(cmd, shell=True, text=True, universal_newlines=True)
        
        
     
    #~~~command: bedtools for repeats
    command  = f"bedtools intersect -s -C -a {transd_temp_file.name}.tmp.sorted.bed  -b {RM_temp_file.name}.tmp.sorted.bed"
    result = subprocess.run(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, universal_newlines=True)
    
    #~~~ verify successful execution and store the non-overlapping transductions in a list to extract them from the potential transductions.
    to_grep_RM = []
    if result.returncode == 0:
        #split the stdout into la list of lines and iterate over them
        for line in result.stdout.splitlines():
            line = line.split("\t")
            # Check if there is no overlap with eepeats
            if line[-1] == "0":
                # Making a specific keyword to extract those transductions not overlapping with Repeat from the  
                # .. transduction inpt file
                key = line[0:3]
                formated_key = f'{key[0]}:{(int(key[1]) + 1)}-{key[2]}'   #['chrY', '10884095', '10886133'] ---> chrY:10884096-10886133
                to_grep_RM.append(formated_key)
    else:
        print("Error:", result.stdout)
    
    
    line_match = filter(lambda line: any(locus in line for locus in to_grep_RM), raw_input_transduction.strip().split("\n"))
    sorted_line_match = sorted([line for line in line_match], key=custome_sort_key)

    with open(os.path.join(path, "potential_alu_transd_with_pA_heteroTSD_noTE.tsv"), "w") as f:
        f.write(header + "\n")
        for i in sorted_line_match:
            f.write(i + "\n")
        
    print(f"{Style.BRIGHT}{Fore.YELLOW}Filtering out transductions overlapping with TEs completed.")
    time.sleep(1)
     

    #~~~ Deleting unneeded variables 
    del(bed_transduction)
    del(bed_RM)
    
    #~~~ Deleting temporary files
    temp_files_pattern = os.path.join(path, "tmp*")
    temp_files = glob.glob(temp_files_pattern)
    for f in temp_files:
        try:
            os.remove(f)
        except Exception as e:
            print(f"{e}")
    
            
    
    #~~~ alignment process
    
    print(f"{Style.BRIGHT}{Fore.YELLOW}Initializing the alignment...")
    time.sleep(1)

    if not os.path.exists(path + "/alignment"):
        os.makedirs(path + "/alignment")
        
    alignment_folder = path + "/alignment"
        
    
    # Read the identified_transduction_file, generate headers and sequencing coordinates,
    # and save them in a list of tuples. This list will be used to create a FASTA file 
    # for extraction in the next step.
    print(f"{Style.BRIGHT}{Fore.YELLOW}Preparing the inputs for the alignment...")
    time.sleep(1)
    
    
    trd_no_TE_overlap = os.path.join(path, "potential_alu_transd_with_pA_heteroTSD_noTE.tsv")
    transduction_list = []
    with open(trd_no_TE_overlap) as trd_file:
        for line in trd_file:
            if line.startswith("chr"):
                line_stripped = line.strip()
                line_split = line_stripped.split("\t")
                pos_split = regex.split(r"\W", line_split[3]) # chr1:442696-445883 ==> ["chr1", "442696", "445883"]
                my_header = line_split[2] + "_" + line_split[3] + "_" + line_split[1] # AluY_chr1:442696-445883_-
                point = namedtuple("transd_info", "header, chrom, beg, end, strand")
                transd_info = point(my_header, pos_split[0], int(pos_split[1]), int(pos_split[2]), line_split[1])
                transduction_list.append(transd_info)
    
    
    
    # Reading the ref_genome to extract sequences for queries (transductions)
    with open(ref_genome) as rf_genome, \
        open(alignment_folder +  "/query_alu_transductions.fa", "w") as q_output:
            for record in SeqIO.parse(rf_genome, "fasta"):
                for trd in transduction_list:
                    if trd.chrom == record.id:    
                        if trd.strand == "+":
                            q_output.write(">" + trd.header + "\n")
                            q_output.write(str(record.seq[trd.beg-1:trd.end]) + "\n")
                        else: # minus strand
                            q_output.write(">" + trd.header + "\n")
                            #making reverse complement for transductions  
                            rc_seq = Seq(record.seq[trd.beg:trd.end]).reverse_complement()
                            q_output.write(str(rc_seq) + "\n")
    
    
    # Preparing the target list (list of named tuples)
    target_list = []
    
    for alus in read_repeatmasker_file(RM_file, alu_type):
        for keys, values in alus.items():
            for record in values:
                chromosome, start, finish, strand, te_type = record
                my_header = te_type + "_" + chromosome + ":" + str(start) + "-" + str(finish) + "_" + strand
                point = namedtuple("alu_info", "header, chrom, beg, end, strand")
                alu_info = point(my_header, chromosome, start, finish, strand)
                target_list.append(alu_info)
    
    
    
    # extracting of subject sequences (target)
    with open(ref_genome) as rf_genome, \
        open(alignment_folder + "/target_alus.fa", "w") as t_output:
            for record in SeqIO.parse(rf_genome, "fasta"):
                for target in target_list:
                    if target.chrom == record.id:
                        if target.strand == "+":
                           t_output.write(">" + target.header + "\n")
                           t_output.write(str(record.seq[target.end:target.end+4500]) + "\n")
                        else:  #reverse strand
                            t_output.write(">" + target.header + "\n")
                            rc_seq = Seq(record.seq[(target.beg-1)-4500:target.beg-1]).reverse_complement()
                            t_output.write(str(rc_seq) + "\n")
    
    
    
    # run pblat
    print("Starting the alignment process ...")  
    
    my_target = alignment_folder + "/target_alus.fa"
    my_query = alignment_folder +  "/query_alu_transductions.fa"
    my_output = alignment_folder +  "/rawblatResults.pslx"
    
    
    
    subprocess.run(["pblat",  f"-threads={num_cpus}",  "-stepSize=1", "-oneOff=2", "-minScore=20", "-repMatch=4096", "-out=pslx", "-noHead", my_target, my_query, my_output])
    print(f"{Style.BRIGHT}{Fore.YELLOW}Alignment completed.")
    
    
    

    #~~~
    #~~~ Filtering the alignment output
    #~~~
    print(f"{Style.BRIGHT}{Fore.YELLOW}Filtering out the alignments and generating the final output...")
    pblat_output = alignment_folder +  "/rawblatResults.pslx"
    
    # generating a temp file containing the verified sources
    with open(pblat_output) as in_file, \
        open(os.path.join(alignment_folder, "temp1_filtered.output.pslx"), "w") as out_file:
             
            genome_dict = SeqIO.index(ref_genome, "fasta")  #creates a dictionary out of ref genome
            for line in in_file:
                line_striped = line.strip()
                line_split = line_striped.split("\t")
                Q = line_split[9]
                T = line_split[13]
                if "chrUn" not in Q and "NW" not in Q and "random" not in Q:
                    q = query_instanciation(Q)
                if "chrUn" not in T and "NW" not in T and "random" not in T:
                    t = target_instanciation(T)
                results = filter_blat_record(line_split, q, t, genome_dict)
                for elem in results:
                    line_to_write = "\t".join(_ for _ in elem)
                    out_file.write(line_to_write + "\n")
    
    
    # adding info about the identified terminator within the transduced segment (source element) and 1kb downstream of the transduced seq
    with open(os.path.join(alignment_folder, "temp1_filtered.output.pslx")) as in_file, \
        open(os.path.join(alignment_folder, "temp2_filtered.output.pslx"), "w") as out_file:
            

            for line in in_file:
                line = line.strip()
                line = line.split("\t")
                T = line_split[13]
                #print(line)
                if "chrUn" not in T and "NW" not in T and "random" not in T: 
                    t = target_instanciation(T)
                    result = find_terminator(line, t.strand, line[25], line[26]) # line[25]: transduced seq from source, line[26]: 1kb downstream of thetransduced seq relative to the source
                #check if any terminator was found within the transduced seq; if so, write the coordinate to the final file; if not, write "Notfound"
                    if len(result[0]) == 0 and len(result[1])  > 0:
                        line.append(", ".join(_ for _ in ["Not found", str(result[1])]))
                    elif len(result[0]) > 0 and len(result[1]) == 0:
                        line.append(", ".join(_ for _ in [str(result[0]), "Not found"]))
                    elif len(result[0]) == 0 and len(result[1]) == 0:
                        line.append(", ".join(_ for _ in ["Not found", "Not found"]))
                    else:
                        line.append(", ".join(str(_) for _ in [result[0], result[1]]))
            
                    out_file.write("\t".join(_ for _ in line) + "\n")
    
    
    # check if a terminator is present either within the transduced segment from the source or in the 1 kb downstream sequence.
    # we only want to keep those hits where polyT is NOT found within the transduced segment BUT  present within the 1 kb downstream sequence.
    # additionally, in this step, we calculate the score and identity percentage using the pslScore script from UCSC.
    # the correct script for the operating system from UCSC must have been downloaded.
    
    with open(os.path.join(alignment_folder, "temp2_filtered.output.pslx")) as in_file, \
        open(os.path.join(alignment_folder, "temp3_filtered.output.pslx"), "w") as out_file:
            
            for line in in_file:
                line = line.strip()
                line = line.split("\t")
                polyT_col = regex.split(r", ", line[27])
                if polyT_col[0] == "Not found" and polyT_col[1] != "Not found":
                    out_file.write("\t".join(_ for _ in line) + "\n")
                    
                    
    with open(os.path.join(alignment_folder, "temp3_filtered.output.pslx")) as in_file, \
        open(os.path.join(alignment_folder, "TEMP.pslx"), "w") as out_file:
            #pslScore needs 23 cols from psl file, hence we make a temp file out of our filtered queries
            for line in in_file:
                entire_line = line.strip().split("\t")
                line_for_psl_score = "\t". join(i for i in entire_line[0:23])
                out_file.write(line_for_psl_score + "\n")

    #run pslscore script, keep the results in a variable
    psl_score = args.pslScore
    psl_file = os.path.join(alignment_folder, "TEMP.pslx")
    scores = subprocess.run([psl_score, psl_file], capture_output=True, text=True).stdout.strip("\n")            
    
            
    #add the score and identity percent to the beginning of the filtered blat results
    with open(os.path.join(alignment_folder, "temp3_filtered.output.pslx")) as in_file, \
        open(os.path.join(alignment_folder, "temp4_filtered.output.pslx"), "w") as out_file:
            for line in in_file:
                line_blat = line.strip().split("\t")
                line_score = scores.split("\n")
                for elem in line_score:
                    col = elem.split("\t")
                    query_split = regex.split(r":", col[3])   # AluY_chr1:191553735-191553810_-:0-53 -----> ["AluY_chr1", "191553735-191553810_-", "0-53"]
                    query = query_split[0] + ":" + query_split[1]   # i.e: ["AluY_chr1", "191553735-191553810_-", "0-53"]  --------> AluY_chr1:191553735-191553810_-
                    if col[0] == line_blat[13] and query == line_blat[9]:
                       out_file.write(col[4] + "\t" + col[5] + "\t" + "\t".join(_ for _ in line_blat) + "\n") 
            
    os.remove(os.path.join(alignment_folder, "TEMP.pslx"))
    
    # Keep the best hits and generate the final output
    with open(os.path.join(alignment_folder, "temp4_filtered.output.pslx")) as in_file, \
        open(os.path.join(alignment_folder, "final_filtered.output.pslx"), "w") as out_file:
            
            my_dict = {}
            for line in in_file:
                line = line.strip()
                line = line.split("\t")
                #making a dict keys=query, values=list of hits
                key_dict = line[11]
                hits = my_dict.get(key_dict, [])
                hits.append(line)
                sorted_hits = sorted(hits, key=lambda my_list: my_list[0], reverse=True) #for each query: sort based on the best score (the first column)
                my_dict[key_dict] = sorted_hits
                
            for query, hits in my_dict.items():
                line_to_write = "\t".join(_ for _ in hits[0]) #hits[0] is the best hit 
                out_file.write(line_to_write + "\n")
    print(f"{Style.BRIGHT}{Fore.YELLOW}Filtering out the alignments and generating the final output is done.")
    time.sleep(2)
    
    
    #~~~ Generating final output containing transductions whose source elements have been found
    print(f"{Style.BRIGHT}{Fore.YELLOW}Generating final output containing transductions whose source elements have been found...")
    time.sleep(2)
    
    
    file_base_name = Path(ref_genome).name.split(".")[0]
    
    # converting the potential transduction file to a panda dataframe
    # converting the potential transduction file to a panda dataframe
    transd_pd = pd.read_csv(transd_file, sep="\t", header=None, skiprows=1)
    
    # naming the "transd_pd" in place
    transd_pd.rename(columns={
        0 : "alu_pos",
        1 : "strand",
        2 : "alu_type",
        3 : "transd_pos",
        4 : "transd_len",
        5 : "TSD_seq",
        6 : "TSD_len",
        7 : "left_TSD_pos",
        8 : "right_TSD_pos",
        9 : "polyA_seq",
        10: "polyA_pos"
        }, inplace=True)
    
    # converting the filtered alignment file to a panda dataframe
    source_df = pd.read_csv(final_blat, sep="\t", header=None)

    #~~~ preprocess steps to merge two dataframes
    new_col = []

    for idx in source_df.index:
        value = source_df[11].iloc[idx]
        split_val = value.split("_")
        new_col.append(split_val[1])
        
    source_df.insert(0, "use_to_merge", new_col)

    #merge two dataframes
    transd_source_df = transd_pd.merge(source_df, left_on="transd_pos", right_on="use_to_merge")


    #~~~
    col_name  = [i for i in range(11)] # from 0 to 10 inclusive

    for i in [27, 38, 39, -1]:
        col_name.append(i)


    # Make a copy of the sliced DataFrame to avoid SettingWithCopyWarning
    summary_transd_source_df = transd_source_df.iloc[:, col_name].copy()
        
    #summary_transd_source_df = tmp_summary_transd_source_df.iloc[:, [_ for _ in col_name]]


    # pick up the first ployT within the source  and updating for idx in summary_transd_source_df.index in place
    for idx in summary_transd_source_df.index:
        val = summary_transd_source_df[29].iloc[idx]
        temp_val = val.split("d, ")[1].strip()  #split Not foun"d, "[(...)]
        final_val = ast.literal_eval(temp_val)[0]  #turn the string into python object such as list
        summary_transd_source_df.at[idx, 29] = final_val
        
    #rename the cols
    summary_transd_source_df.rename(columns={
        15 : "Alu_source_info",
        26 : "transduced_pos_from_source",
        27 : "aligned_part",
        29 : "terminator_info (start, end, seq) (coordinates relative to distance after aligned part to the source)"
        }, inplace=True)

    # adjusting the length of transductions
    for idx in summary_transd_source_df.index:
        new_len = int(summary_transd_source_df["transd_len"].iloc[idx]) + 1
        summary_transd_source_df.at[idx, "transd_len"] = new_len
        
        
    #~~~ extracting the sequence of transduction and adding it to our dataframe
    my_ref_genome = SeqIO.index(genome_path, "fasta")

    summary_transd_source_df["transduction_seq"] = None

    for idx in summary_transd_source_df.index:
        strand = summary_transd_source_df["strand"].iloc[idx]
        pos = summary_transd_source_df["transd_pos"].iloc[idx]
        pos = regex.split(r"[:-]", pos)
        chrom, beg, end = pos[0], int(pos[1]), int(pos[2])
        
        if strand == "+":
            seq = my_ref_genome[chrom][beg-1:end].seq
            summary_transd_source_df.at[idx, "transduction_seq"] = seq
        else:
            seq = my_ref_genome[chrom][beg-1:end].seq.reverse_complement()
            seq = seq[::-1]
            summary_transd_source_df.at[idx, "transduction_seq"] = seq
            
    final_df = summary_transd_source_df.drop(columns=["transduced_pos_from_source", "aligned_part", "terminator_info (start, end, seq) (coordinates relative to distance after aligned part to the source)"]).copy()
    
    
    #~~~ generating the final output 
    final_df.to_csv(
        os.path.join(path, file_base_name + "_verified_transductions.tsv"),
        sep="\t",
        index=False
        )
    
    
    
    print(f"{Style.BRIGHT}{Fore.YELLOW}Analysis completed.")
    time.sleep(2)
    
    stop_time = timeit.default_timer()
    print(f"{Back.CYAN}Elapsed time  for the entire processing: {round((stop_time - start_time) / 60, 2)} min")

if __name__ == "__main__":
    main()
            