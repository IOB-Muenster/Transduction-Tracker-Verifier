#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import re
from Bio.Seq import Seq
from Bio import motifs
import regex



# 
# -----------------------------------------------------------------------------
# Function: generates overlapping  kmers for upstream sequence (forward strand)
# -----------------------------------------------------------------------------
# 

def build_kmer_collection_f_upstream(seq, interval1=5, interval2=45):

    seq = seq.upper()
    
    # a dict in which keys are overlapping k-mers and values are the coordinates of k-mers
    k_mers_f_up = {}
    
    for k in range(interval1, interval2 + 1):
        for i in reversed(range(0, len(seq) - k + 1)):
            kmer = seq[i:i+k]
            kmer_start = i
            kmer_end = i + k
            list_coordinates = k_mers_f_up.get(kmer, [])
            list_coordinates.append((kmer_start, kmer_end))
            k_mers_f_up[kmer] = list_coordinates
    
    return k_mers_f_up

# 
# -----------------------------------------------------------------------------
# Function: generates overlapping kmers for downstream sequence(forward strand)
# -----------------------------------------------------------------------------
# 

def build_kmer_collection_f_downstream(seq, interval1=5, interval2=45):
   
    seq = seq.upper()
    
    # a dict in which keys are overlapping k-mers and values are the coordinates of k-mers
    k_mers_f_down = {}
    
    for k in range(interval1, interval2 + 1):
        for i in range(0, len(seq) - k + 1):
            kmer = seq[i:i+k]
            kmer_start = i
            kmer_end = i + k
            list_coordinates = k_mers_f_down.get(kmer, [])
            list_coordinates.append((kmer_start, kmer_end))
            k_mers_f_down[kmer] = list_coordinates
    
    return k_mers_f_down

# 
# -----------------------------------------------------------------------------
# Function: generates overlapping kmers using upstream sequence (reverse strand)
# -----------------------------------------------------------------------------
# 

def build_kmer_collection_r_upstream(seq, interval1=5, interval2=45):
    
    
    seq = seq.upper()
    
    # a dict in which keys are overlapping k-mers and values are the coordinates of k-mers
    k_mers_r_up = {}
    
    for k in range(interval1, interval2 + 1):
        for i in range(0, len(seq) - k + 1):
            kmer = seq[i:i+k]
            kmer_start = i
            kmer_end = i + k
            list_coordinates = k_mers_r_up.get(kmer, [])
            list_coordinates.append((kmer_start, kmer_end))
            k_mers_r_up[kmer] = list_coordinates
    
    return k_mers_r_up

# 
# -----------------------------------------------------------------------------
# Function: generates overlapping kmers using downstream sequence (reverse strand)
# -----------------------------------------------------------------------------
# 

def build_kmer_collection_r_downstream(seq, interval1=5, interval2=45):
    
    seq = seq.upper()
    
    # a dict in which keys are overlapping k-mers and values are the coordinates of k-mers
    k_mers_r_down = {}
    
    for k in range(interval1, interval2 + 1):
        for i in reversed(range(0, len(seq) - k + 1)):
            kmer = seq[i:i+k]
            kmer_start = i
            kmer_end = i + k
            list_coordinates = k_mers_r_down.get(kmer, [])
            list_coordinates.append((kmer_start, kmer_end))
            k_mers_r_down[kmer] = list_coordinates
    
    return k_mers_r_down

# 
# -----------------------------------------------------------------------------
# Function: calculating hamming distance between two kmers
# -----------------------------------------------------------------------------
# 

def ham_distance(kmer_1, kmer_2):
    # calculating hamming distance for two sequences except for the first and last bases of each sequence
    set_1 = set( (index, value) for index, value in enumerate(kmer_1[1:-1].upper()))
    set_2 = set( (index, value) for index, value in enumerate(kmer_2[1:-1].upper()))
    
    return len(set_1.difference(set_2))

# 
# -----------------------------------------------------------------------------
# Function: Forward strand: polyA tail finder
# -----------------------------------------------------------------------------
# 

def poly_A_tracker_V2_f_strand(seq, pos):
    pA = []
    counter = 0

    # iterate backward from the position before the 3' TSD
    for i in range(pos - 1, -1, -1):
        if seq[i].upper() == "A":
            pA.append(seq[i])
            counter = 0 # Reset counter when an 'A' is found
        else:
            counter += 1
            if counter > 2:
                break # Stop if more than 2 non-'A' nucleotides are found
            pA.append(seq[i]) # Include the non-'A' nucleotide

    # polyA must start with min two As
    polyA = "".join(reversed(pA))
    pattern = r"A{2}.+$"
    for m in regex.finditer(pattern, polyA, regex.IGNORECASE):
        beg = pos - len(m.group())
        end = pos
        return m.group(), (beg, end), len(m.group())
            
            
            
# 
# -----------------------------------------------------------------------------
# Function: Reverse: finding polyA tail
# -----------------------------------------------------------------------------
# 

def poly_A_tracker_V2_r_strand(seq):
    pA = []
    counter = 0
    
    # iterate backward from the position before the 3' TSD
    for i in range(len(seq)):
        if seq[i].upper() == "T":
            pA.append(seq[i])
            counter = 0 # Reset counter when an 'A' is found
        else:
            counter += 1
            if counter > 2:
                break # Stop if more than 2 non-'A' nucleotides are found
            pA.append(seq[i]) # Include the non-'A' nucleotide

    # polyA must start with min two As
    polyA = "".join(pA)
    pattern = r".+T{2}"
    for m in regex.finditer(pattern, polyA, regex.IGNORECASE):
        #coordinate = (m.span()[0] + pos, m.span()[1]+ pos)
        beg = 0
        end = len(m.group())
        return m.group(), (beg, end), len(m.group())    

# 
# -----------------------------------------------------------------------------
# Function: Forward strand: processes each pair of up and downstream sequences 
# of a TE in the RM to identify the putative TSDs, polyAs and transductions
# -----------------------------------------------------------------------------
# 

def find_putative_TSD_f_V6(chrom, start, end, strand, te_type, upstream_seq, downstream_seq, num_upstream_bases, num_downstream_bases):
    
    # generating dicts for up and downstream sequence
    dict_kmers_upstream_f = build_kmer_collection_f_upstream(upstream_seq)
    dict_kmers_downstream_f = build_kmer_collection_f_downstream(downstream_seq)
    
    # finding upstream-downstream matched kmers (1 mismatch is allowed)
    matched_kmers_f = [
        (len(kmer_up), kmer_up, dict_kmers_upstream_f[kmer_up][0], kmer_down, dict_kmers_downstream_f[kmer_down][0]) 
        for kmer_up in dict_kmers_upstream_f.keys()
        for kmer_down in dict_kmers_downstream_f.keys()
        if len(kmer_up) == len(kmer_down)
        if (kmer_up[0] == kmer_down[0] and kmer_up[-1] == kmer_down[-1])
        if ham_distance(kmer_up, kmer_down) <= 1
        ]
     
    # picking the closest kmer pairs to TE and identifying the TSD                    
    sorted_by_pos_f = sorted(matched_kmers_f, key=lambda pos: pos[2][1], reverse=True)
    temp_list_f = [i for i in sorted_by_pos_f if i[2][1] == sorted_by_pos_f[0][2][1]]
    for i in temp_list_f:
        if len(i) != 5:
            print("Forward", i, len(i))
        
    tsd = None   # to avoid UnboundLocalError: for tsd defined in the try block
    try:
        tsd = sorted(temp_list_f, key=lambda i: i[0], reverse=True)[0]
    except IndexError:
        pass
    
    # unpacking the resultant tuple consisting of TSD information
    finally:
        if tsd: 
            seq_object = [Seq(tsd[1]), Seq(tsd[3])]
            motif = motifs.create(seq_object)
            tsd_len, tsd_seq = tsd[0], str(motif.degenerate_consensus)  
            left_tsd_start, left_tsd_end = tsd[2][0], tsd[2][1]
            right_tsd_start, right_tsd_end = tsd[4][0], tsd[4][1]
    
            polyA_search_end_f = right_tsd_start
            
            polyA_f_temp = poly_A_tracker_V2_f_strand(downstream_seq[:polyA_search_end_f], polyA_search_end_f)     
            if polyA_f_temp: 
                polyA_f = polyA_f_temp[0]
                polyA_f_beg = ((end-25)+1+right_tsd_start) - polyA_f_temp[2]  #polyA_f_temp[2] is length polyA
                polyA_f_end = (end-25)+right_tsd_start
            else: 
                polyA_f = 'PolyA_tail_not_found'
            
            
            
            if not polyA_f.startswith('PolyA_tail'):
                line_columns_to_write = (
                    chrom, str(start), str(end), strand, te_type, tsd_seq, str(tsd_len), 
                    str((start-num_upstream_bases)+left_tsd_start) + '-' + str((start-1-num_upstream_bases)+left_tsd_start+tsd_len),
                    str((end-25)+1+right_tsd_start) + '-' + str((end-25)+right_tsd_start+tsd_len), polyA_f,
                    str(polyA_f_beg) + "-" + str(polyA_f_end)
                    )
            else:
                line_columns_to_write = (
                    chrom, str(start), str(end), strand, te_type, tsd_seq, str(tsd_len), 
                    str((start-num_upstream_bases)+left_tsd_start) + '-' + str((start-1-num_upstream_bases)+left_tsd_start+tsd_len),
                    str((end+1)+right_tsd_start) + '-' + str((end-25)+right_tsd_start+tsd_len), polyA_f,
                    '-'
                    )
            
            line_to_write = '\t'.join(i for i in line_columns_to_write)
                
            return line_to_write
        else:
            line_to_write = chrom + '\t' +  str(start) + '\t' + str(end) + '\t' + strand + '\t' + te_type + '\t' + 'TSD_not_found'
            return line_to_write
    
# 
# ------------------------------------------------------------------------------
# Function: Reverse strand: processes each pair of up and downstream sequences 
# of a TE in the RM to identify the putative TSDs, polyAs and transductions
# ------------------------------------------------------------------------------
# 

def find_putative_TSD_r_V6(chrom, start, end, strand, te_type, upstream_seq, downstream_seq, num_upstream_bases, num_downstream_bases):
    
    # generating dicts for up and downstream sequence
    dict_kmers_upstream_r = build_kmer_collection_r_upstream(upstream_seq)
    dict_kmers_downstream_r = build_kmer_collection_r_downstream(downstream_seq)
    
    # finding upstream-downstream matched kmers (1 mismatch is allowed)
    matched_kmers_r = [
        (len(kmer_down), kmer_down, dict_kmers_downstream_r[kmer_down][0], kmer_up, dict_kmers_upstream_r[kmer_up][0]) # 0 is the closest Kmer to the TE
        for kmer_up in dict_kmers_upstream_r.keys()
        for kmer_down in dict_kmers_downstream_r.keys()
        if len(kmer_up) == len(kmer_down)
        if (kmer_up[0] == kmer_down[0] and kmer_up[-1] == kmer_down[-1])
        if ham_distance(kmer_up, kmer_down) <= 1
        ]
    
    # picking the closest kmer pairs to TE and identifying the TSD                    
    sorted_by_pos_r = sorted(matched_kmers_r, key=lambda pos: pos[4][0])
    temp_list_r = [i for i in sorted_by_pos_r if i[4][0] == sorted_by_pos_r[0][4][0]]
    for i in temp_list_r:
        if len(i) != 5:
            print("reverse", i, len(i))
    tsd = None # to avoid UnboundLocalError: for tsd defined in the try block
    try:
        tsd = sorted(temp_list_r, key=lambda i: i[0], reverse=True)[0]
    except IndexError:
        pass
    
    # unpacking the resultant tuple consisting of TSD information
    finally:
        if tsd: 
            seq_object = [Seq(tsd[1]), Seq(tsd[3])]
            motif = motifs.create(seq_object)
            tsd_len, tsd_seq = tsd[0], str(motif.degenerate_consensus)  
            left_tsd_start, left_tsd_end = tsd[2][0], tsd[2][1]
            right_tsd_start, right_tsd_end = tsd[4][0], tsd[4][1]
    
            polyA_search_end_r = left_tsd_end
            
            polyA_r_temp = poly_A_tracker_V2_r_strand(downstream_seq[polyA_search_end_r:])
        
            if polyA_r_temp: 
                polyA_r = str(Seq(polyA_r_temp[0]).reverse_complement())
                polyA_r_beg = ((start+25)-(num_downstream_bases+25))+left_tsd_start+tsd_len
                polyA_r_end = (((start+25)-1-(num_downstream_bases+25))+left_tsd_start+tsd_len) + polyA_r_temp[2]
            else:
                polyA_r = 'PolyA_tail_not_found'
            
     
            if not polyA_r.startswith('PolyA_tail'):
                line_columns_to_write = (
                    chrom, str(start), str(end), strand, te_type,
                    str(Seq(tsd_seq).reverse_complement()), str(tsd_len),
                    str(((start+25)-(num_downstream_bases+25))+left_tsd_start) + '-' +
                    str(((start+25)-1-(num_downstream_bases+25))+left_tsd_start+tsd_len),
                    str(end+1+right_tsd_start) + '-' + str(end+right_tsd_start+tsd_len),polyA_r, 
                    str(polyA_r_beg) + "-" + str(polyA_r_end)
                    )
            else:
                line_columns_to_write = (
                    chrom, str(start), str(end), strand, te_type,
                    str(Seq(tsd_seq).reverse_complement()), str(tsd_len),
                    str(((start+25)-(num_downstream_bases+25))+left_tsd_start) + '-' + 
                    str(((start+25)-1-(num_downstream_bases+25))+left_tsd_start+tsd_len),
                    str(end+1+right_tsd_start) + '-' + str(end+right_tsd_start+tsd_len),polyA_r, 
                    '-'
                    )
            line_to_write = '\t'.join(i for i in line_columns_to_write)
            return line_to_write
        else:
            line_to_write = chrom + '\t' +  str(start) + '\t' + str(end) + '\t' + strand + '\t' + te_type + '\t' + 'TSD_not_found'
            return line_to_write


# 
# -----------------------------------------------------------------------------
# Function: detecting transductions caused by Alu elements with detected pA
# -----------------------------------------------------------------------------
# 

def process_alu_with_polyA(output1, output2, line, dist=70):
    
    strand = line[3]
    if strand == '+':
        #calculates the length between the end of Alu and the start of polyA
        temp_len = int(line[10].split('-')[0])-int(line[2])                    
        if temp_len > dist:
            transD_coordinate = line[0] + ":" + str(int(line[2]) + 1) + "-" + str(int(line[2]) + temp_len -1)
            transD_len = (int(line[2]) + temp_len -1) - (int(line[2]) + 1)
            rows = [line[0]+ ":" + line[1] + "-" + line[2], line[3], line[4], 
                    transD_coordinate, transD_len ,line[5], line[6], line[7], 
                    line[8], line[9], line[10]
                    ]
            
            if len(line[9]) >= 5:
                output1.write("\t".join(str(i) for i in rows) + '\n')     
            else:    #if the length pA is less than 5 (too short), write it into no pA elements
                output2.write("\t".join(str(i) for i in rows) + '\n')
    else: #for minus strand
        temp_len = int(line[1]) - int(line[10].split('-')[1])
        if temp_len > dist:
            transD_coordinate = line[0] + ":" + str((int(line[1])-temp_len) + 1) + "-" + str(int(line[1]) - 1)
            transD_len = (int(line[1]) - 1) - (int(line[1])-temp_len + 1)
            rows = [line[0]+ ":" + line[1] + "-" + line[2], line[3], line[4], 
                    transD_coordinate, transD_len, line[5], line[6], line[7], 
                    line[8], line[9], line[10]
                    ]
            if len(line[9]) >= 5:
                output1.write("\t".join(str(i) for i in rows) + '\n')
            else:    #if the length pA is less than 5 (too short), write it into no pA elements
                output2.write("\t".join(str(i) for i in rows) + '\n')

# 
# -----------------------------------------------------------------------------
# Function: detecting transductions caused by Alu elements lacking a detectable
# poly-A tail
# -----------------------------------------------------------------------------
# 

def process_alu_without_polyA(output2, line, dist=30):
    strand = line[3]
    if strand == '+':
        #calculates the length between the end of Alu and the start of polyA
        temp_len = int(line[8].split('-')[0])-int(line[2])                    
        if temp_len > dist:
            transD_coordinate = line[0] + ":" + str(int(line[2]) + 1) + "-" + str(int(line[2]) + temp_len -1)
            transD_len = (int(line[2]) + temp_len -1) - (int(line[2]) + 1)
            rows = [line[0]+ ":" + line[1] + "-" + line[2], line[3], line[4], 
                    transD_coordinate, transD_len ,line[5], line[6], line[7], 
                    line[8], line[9], line[10]
                    ]
            output2.write("\t".join(str(i) for i in rows) + '\n')
    else:
        temp_len = int(line[1]) - int(line[7].split('-')[1])
        if temp_len > dist:
            transD_coordinate = line[0] + ":" + str((int(line[1])-temp_len) + 1) + "-" + str(int(line[1]) - 1)
            transD_len = (int(line[1]) - 1) - (int(line[1])-temp_len + 1)
            rows = [line[0]+ ":" + line[1] + "-" + line[2], line[3], line[4], 
                    transD_coordinate, transD_len, line[5], line[6], line[7], 
                    line[8], line[9], line[10]
                    ]
            output2.write("\t".join(str(i) for i in rows) + '\n')
