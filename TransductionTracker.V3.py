#!/usr/bin/env python3
# -*- coding: utf-8 -*-


#~~~ python bulit-in modules
from Bio import SeqIO
import os
import timeit
from pathlib import Path
from colorama import Fore, Style, init, Back
from joblib import Parallel, delayed
import argparse
import textwrap
import time
#from tqdm import tqdm

#~~~ internal modules
from modules.utilities_TransductionTracker_V2 import *
from modules.core_functions_TransductionTracker_V3 import *
from modules.utilities_TransductionVerifier_V2 import * 
from modules.repeat_segdup_filter import *



# 
# -----------------------------------------------------------------------------
# Code
# -----------------------------------------------------------------------------
# 

def main():
    
    #~~~ Initialize colorama
    #~~~ To ensure every print statement automatically resets the color to the default after it's called, so there is no need to manually reset it each time.
    init(autoreset=True) 
    
    
    #~~~ Setting up the argument parser with a description, usage, and epilog for help documentation
    my_parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent("""\
                                    Locates Target Site Duplicates (TSDs) for Alus across a reference genome, utilizing specified repeat annotations.
                                    """),
        usage=textwrap.dedent("""python TransductionTracker.V3.py -r <repeat_annotation_path> -g <reference_genome> -s <segmental_duplications (optional)> -o <output_directory> -u <upstream_flank> -d <downstream_flank> -a <Alu_type> -c <cpu_number>
                              """),
        epilog=textwrap.dedent("""
                               Example: python TransductionTracker.V3.py -r repeats.bed -g genome.fasta -s seg_dup.bed -o output_directory -u 100 -d 4500 -a AluY -c 4
                               """
        )
        )
        
    #~~~ Creating argument groups for better organization
    file_group = my_parser.add_argument_group("Inputs Information")
    output_path = my_parser.add_argument_group("Output Directory")
    analysis_group = my_parser.add_argument_group("Analysis Parameters")
    performance_group = my_parser.add_argument_group("Performance")
    
    
    #~~~ Adding command-line arguments for the script
    file_group.add_argument("-r", "--repeat_annotation", type=str, required=True, metavar="",   
                           help="Path to the repeat annotation file.")
    # file_group.add_argument("-f", "--repeat_annotation_format", type=str, choices=["gff", "bed", "RM.out"], required=True, metavar="",   
    #                        help="Type of repeat annotation format (bed, gff, or RM.out).\n"
    #                        "Recommend using RM.out because the script can automatically extract full-length elements using the length information in RM.out.\n"
    #                        "For other formats i.e gff and bed, ensure full-length elements are extracted.\n")
    file_group.add_argument("-g", "--reference_genome", type=str, required=True, metavar="",
                           help="path to the reference genome in fasta format.")
    file_group.add_argument("-s", "--segmental_duplication", type=str, default=None, metavar="",
                           help="Path to the Segmental Duplication (bed format) file (optional).")
    output_path.add_argument("-o", "--output_directory", type=str, required=True, metavar="",
                             help="Path to the directory where output files will be saved.")
    analysis_group.add_argument("-u", "--upstream_base_count", type=int, default=100, metavar="",
                           help="Number of upstream bases to include in the search for 5' TSDs. The default is 100.")
    analysis_group.add_argument("-d", "--downstream_base_count", type=int, default=4500, metavar="",
                           help="Number of downstream bases to include in the search for 3' TSDs. The default is 4500.")
    analysis_group.add_argument("-a", "--alu_type", type=str, default="AluY", metavar="",
                               help="Type of Alu element to search for. The default is 'AluY'.")
    performance_group.add_argument("-c", "--cpus", type=int, default=4, metavar="",
                           help="Number of CPUs to use for processing. The default is 4.")
    # performance_group.add_argument("-v", "--verbosity_level", type=int, choices=range(1, 11), default=10, metavar="",
    #                        help="Set the verbosity level from 1 (least verbose) to 10 (most verbose). Default is 10.")
                                        
    args = my_parser.parse_args()
    
    #~~~ Initialize the script
    print(f"{Fore.YELLOW}=" * 84)
    print(f"{Style.BRIGHT}{Fore.YELLOW}Initializing TSDetectV2...")
    time.sleep(3)
    

    
    #~~~ Preparing for analysis by reading input files and setting up variables
    TE_type = args.alu_type
    #RM_file = read_repeatmasker_file(args.repeat_annotation, TE_type)
    rm_input_path = Path(args.repeat_annotation)
    file_base_name = rm_input_path.name.split(".")[0]
    reference_genome = args.reference_genome
    up_flank = args.upstream_base_count
    down_flank = args.downstream_base_count
    number_of_cpus = args.cpus
    
    if args.segmental_duplication:
        SD_file = Path(args.segmental_duplication)
    else:
        SD_file = None
        
    #verbosity_level = args.verbosity_level
    
    output_path_dir = Path(args.output_directory).resolve() #resolve convert the provided dir to an absolute path
    
    #~~~ creating a folder to write the results
    # Check if the output directory exists, create it if it does not
    output_folder = output_path_dir / ("output_" + file_base_name)
    output_folder.mkdir(parents=True, exist_ok=True)
    
    
    if SD_file:
        print(f"{Style.BRIGHT}{Fore.YELLOW}Filtering out repeats overlapping with segmental duplicates.")
        time.sleep(5)
        filter_SegDups(rm_input_path, SD_file, output_folder)
        segdup_free_RM_out = output_folder / (file_base_name + "_NoSegDup_RM.out")
        RM_file = read_repeatmasker_file(segdup_free_RM_out, TE_type)
    else:
        print(f"{Back.MAGENTA}{Fore.BLACK}Warning: No segmental duplication annotations were found. Please ensure you filter them out manually.")
        time.sleep(5)
        RM_file = read_repeatmasker_file(rm_input_path, TE_type)
        
        

    #~~~ Notifing the user of the analysis start
    print(f"{Style.BRIGHT}{Fore.YELLOW}Reading repeat annotations and preparing analysis...")
    time.sleep(5)
    
    #~~~ Storing TEs in a list
    #~~~ simple_repeat|low_complexity|satellite are filterd out
    te_list = []
    for te in RM_file:
        for chroms, records in te.items():
            for record in records:
                chrom, start, end, strand, te_type = record
                te_list.append(chrom + '\t' + str(start) + '\t' + str(end) + '\t' +  strand + '\t' + te_type)
    
    print(f"{Style.BRIGHT}{Fore.GREEN}Found {len(te_list)} {TE_type}s to analyze.")
    time.sleep(1)
    print(f"{Back.MAGENTA}{Fore.BLACK}Note: The current version is only optimized for Alus. Using it for other elements is not guaranteed.")
    time.sleep(1)
    print(f"{Style.BRIGHT}{Fore.YELLOW}Starting analysis...")
    time.sleep(1)
    
    
    
    
    #~~~ creating a folder to write the results 
    # current_path = os.getcwd()
    # output_folder = os.path.join(current_path, 'output_' + file_base_name)
    # if not os.path.exists(output_folder):
    #     os.mkdir(output_folder)
    
    #~~~  A vraiable keeping the column names for the output     
    header_primary_results_1 = "\t".join(i for i in ['chr', 'offspring_TE_start', 'offspring_TE_end', 'strand', \
                                'subfamily', 'TSD_seq', 'TSD_len', 'left_TSD_pos', 'right_TSD_pos', \
                                'poly(A)_tail', 'poly(A)_tail_pos']) 
        
        
    #~~~ Setting the timer 
    start_time = timeit.default_timer()
    print(f"{Style.BRIGHT}{Fore.YELLOW}This will take a while!")
    
    
    #~~~ Reading the ref genome
    my_genome = SeqIO.index(reference_genome, "fasta")
    
    # a list that keeps track of those TEs whoes flanking DNA is unkown (in case the ref genome is not T2T)
    up_downstream_gap = [] #contains TEs whose up/downstream sequences ontain  "N"
    
    #~~~ Preparing the delayed funcs for parallelization
    #~~~ Exatrcting the flanking DNA
    delayed_fucns = []
    for locus in te_list:
        locus = locus.split("\t")
        chrom, start, end, strand, te_type = locus[0], int(locus[1]), int(locus[2]), locus[3], locus[4]
        stack = []
        if strand == "+":
            stack.append((str(my_genome[chrom][start-1-up_flank:start-1].seq), str(my_genome[chrom][end-25:end+down_flank].seq)))
            f_upstream_seq = stack[0][0]
            f_downsream_seq = stack[0][1]
            if ("N" in f_upstream_seq) or ("N" in f_downsream_seq):
                up_downstream_gap.append("\t".join([chrom, str(start), str(end), strand, te_type]))
            else:
                try:
                    delayed_fucns.append(delayed(find_putative_TSD_f_V6)(chrom, start, end, strand, te_type, f_upstream_seq, f_downsream_seq, up_flank, down_flank))
                except IndexError:
                    pass
        else:  #analyzing those TEs on the reverse strand
            stack.append((str(my_genome[chrom][start-1-down_flank:start-1+25].seq), str(my_genome[chrom][end:end+up_flank].seq)))
            r_upstream_seq = stack[0][1]
            r_downsream_seq = stack[0][0]
            if ("N" in r_upstream_seq) or ("N" in r_downsream_seq):
                up_downstream_gap.append("\t".join([chrom, str(start), str(end), strand, te_type]))
            else:
                try:
                    delayed_fucns.append(delayed(find_putative_TSD_r_V6)(chrom, start, end, strand, te_type, r_upstream_seq, r_downsream_seq, up_flank, down_flank))
                except IndexError:
                    pass
        stack.pop()
         
    
    
    
    #~~~ Execute parallel jobs using the context manager
    with Parallel(n_jobs=number_of_cpus, verbose=10) as p:
        results = p(delayed_fucns)
    
    print(f"{Style.BRIGHT}{Fore.GREEN}TSD detection step is completed. Generating output files...")
    time.sleep(10)
    
    #~~~ Writing those TEs not analyzed due to a gap
    noTSD_dueToGap = output_folder /  (file_base_name + "_ignored_TEs_due_to_UpDownStreamGap.tsv")
    with noTSD_dueToGap.open("w") as discarded_tes:
        for te in up_downstream_gap:
            discarded_tes.write(te + "\n")
            discarded_tes.flush()
    
             
    #~~~ Writing the final output
    primary_TSDs = output_folder / (file_base_name +  "_results_TSDs.tsv")             
    with primary_TSDs.open('w') as output_TSD_handle_file:
        output_TSD_handle_file.write(header_primary_results_1 + '\n')
        for result in results:
            output_TSD_handle_file.write(result + '\n')
            output_TSD_handle_file.flush()
                
    

    
    
    #~~~ Verify whether the detected TSDs consist of a single repeated base or are heterogeneous
    print(f'{Style.BRIGHT}{Fore.GREEN}Starting a sanity check of detected TSDs...') 
    time.sleep(10)
    write_output_file_tsd(primary_TSDs, output_folder, file_base_name, header_primary_results_1)
    print(f'{Style.BRIGHT}{Fore.GREEN}Sanity check of detected TSDs completed.') 
    time.sleep(10)
    
    
    
    #~~~ identifying transductions
    print(f'{Style.BRIGHT}{Fore.GREEN}Initialzing Transdcution detection...') 
    time.sleep(10)
    
    file_name_detected_hetero_TSDs = output_folder / (file_base_name + "_heteroTSD.tsv")
    alu_tsd_with_pA = output_folder / (file_base_name + "_potential_alu_transd_with_pA.tsv")
    alu_tsd_without_pA = output_folder / (file_base_name + "_potential_alu_transd_without_short_pA.tsv")
    
    with  file_name_detected_hetero_TSDs.open("r") as f_handle, \
          alu_tsd_with_pA.open("w") as output1, \
          alu_tsd_without_pA.open("w") as output2:
            
        # col headers
        col_names = ['Alu_coordinate', 'strand', 'subfamily', 'offspring_trandsuction_coordinate', 
                    'transduction_len','TSD_seq', 'tsd_len', 'left_tsd_pos', 'right_tsd_pos', 
                    'polyA', 'polyA_pos'
                    ]
        header = '\t'.join(i for i in col_names)
            
        # writing headers to each output
        output1.write(header + "\n")
        output2.write(header + "\n")
            
            
            
        for line in f_handle:
            line = line.strip()
            if not line.endswith('_pos'):
                line = line.split("\t")
                try:
                    if len(line) == 11:
                        if not line[9] == 'PolyA_tail_not_found':
                            process_alu_with_polyA(output1, output2, line)
                        else:
                            process_alu_without_polyA(output2, line)
                except IndexError:
                    pass
    
    print(f'{Style.BRIGHT}{Fore.GREEN}Transdcution detection completed.') 
    time.sleep(10)
    
    
    stop_time = timeit.default_timer()
    print(f"{Back.CYAN}Elapsed time  for the entire processing: {round((stop_time - start_time) / 60, 2)} min")
    print(f"{Fore.YELLOW}=" * 84)
    

if __name__ == "__main__":
    main()