#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from Bio import SeqIO
from Bio.Seq import Seq
import regex


# 
# -----------------------------------------------------------------------------
# Super class: Retroelement
# -----------------------------------------------------------------------------
# 

class Retroelement(object):
    def __init__(self, chrom, beg, end, strand, family_type):
        self.chrom = chrom
        self.beg = beg
        self.end = end
        self.strand = strand
        self.family_type = family_type


# 
# -----------------------------------------------------------------------------
# Subclass: Query
# -----------------------------------------------------------------------------
#        

class Query(Retroelement):
    #inherits instance variables from superclass
    def query_hit_coordinate(self, line_hit): #line_hit: each line of blat output that is stripped and splitted
        """ Using the blat output, find the coordinate of hit (alignment part on the query)"""
        if self.strand == "+":
            hit_beg = self.beg + int(line_hit[11])
            hit_end = self.beg + int(line_hit[12])
            hit_pos = self.chrom + ":" + str(hit_beg) + "-" + str(hit_end)
            return hit_pos
        
        if self.strand == "-":
            hit_beg = self.end - int(line_hit[12])
            hit_end = self.end - int(line_hit[11])
            hit_pos = self.chrom + ":" + str(hit_beg) + "-" + str(hit_end)
            return hit_pos
        
        
# 
# -----------------------------------------------------------------------------
# Subclass: Target
# -----------------------------------------------------------------------------
# 

class Target(Retroelement):
    #inherits instance variables from superclass
    def target_hit_coordinate(self, line_hit): #line_hit: each line of blat output that is stripped and splitted
        """ Using the blat output, find the coordinate of hit (alignment part on the target)"""
    
        if self.strand == "+":
            hit_beg = (self.end + 1) + int(line_hit[15])
            hit_end = self.end + int(line_hit[16])
            hit_pos = self.chrom + ":" + str(hit_beg) + "-" + str(hit_end)
            return hit_pos
        
        if self.strand == "-":
            hit_beg = self.beg - int(line_hit[16])
            hit_end = (self.beg - 1) - int(line_hit[15])
            hit_pos = self.chrom + ":" + str(hit_beg) + "-" + str(hit_end)
            return hit_pos
        
        
        
    def extract_sequence(self, line_hit, ref_genome):
        """
        Extracts the sequence of the transduced segment and the sequence of 1KB downstream of each transduced part.
        ref_genome: s a dictionary.
        """
        #             trandD       1kb
        #Alu-source >>>>>>>>>>> -----------
        
        if self.strand == "+":
            hit_beg = (self.end + 1) + int(line_hit[15])
            hit_end = self.end + int(line_hit[16])
            transD_seq = ref_genome[self.chrom][hit_beg:hit_end].seq
            downstream_seq = ref_genome[self.chrom][hit_end:hit_end+1000].seq
            return transD_seq, downstream_seq
        
        if self.strand == "-":
            hit_beg = self.beg - int(line_hit[16])
            hit_end = (self.beg - 1) - int(line_hit[15])
            transD_seq = Seq(str(ref_genome[self.chrom][hit_beg:hit_end].seq)).reverse_complement()
            downstream_seq = Seq(str(ref_genome[self.chrom][hit_beg-1000:hit_beg].seq)).reverse_complement()
            return transD_seq, downstream_seq

# 
# -----------------------------------------------------------------------------
# Function: query (object) instantiation
# -----------------------------------------------------------------------------
# 

def query_instanciation(query):
    query = query.split("_") #AluY_chr1:1429303-1430783_+ or AluY_d4_chr1:1429303-1430783_+
    query_pos = [regex.split(r"[:-]", i) for i in query if i.startswith("chr")] #splitting the element in the following list starting with chr ["AluY", "chr1:1429303-1430783", "+" ] or ["AluY", "d4", "chr1:1429303-1430783", "+"]
    q = None
    
    try:
        q = Query(query_pos[0][0], int(query_pos[0][1]), int(query_pos[0][2]), query[-1], query[0])
    except IndexError:
        pass
    
    return q

# 
# -----------------------------------------------------------------------------
# Function: target (object) instantiation
# -----------------------------------------------------------------------------
# 

def target_instanciation(target):
    target = target.split("_")  #AluY_chr1:1429303-1430783_+ or AluY_d4_chr1:1429303-1430783_+
    target_pos = [regex.split(r"[:-]", i) for i in target if i.startswith("chr")] #splitting the element in the following list starting with chr ["AluY", "chr1:1429303-1430783", "+" ] or ["AluY", "d4", "chr1:1429303-1430783", "+"]
    t = None
    
    try:
        t = Target(target_pos[0][0], int(target_pos[0][1]), int(target_pos[0][2]), target[-1], target[0])
    except IndexError:
        pass
    
    return t

# 
# -----------------------------------------------------------------------------
# Function: processing and filtering each hit line
# -----------------------------------------------------------------------------
# 

def filter_blat_record(hit_line, query, target, genome_dictionary):
    filter_results = []
    
    try:
        if not (
            (query.end + 1 == target.beg and 
            hit_line[0] == hit_line[10] and #number_match == query_size
            query.chrom == target.chrom and 
            query.strand == target.strand) 
            or                     
            # filter out self-hits on the plus strands
            (query.beg - 1 == target.end and
            hit_line[0] == hit_line[10] and #number_match == query_size
            query.chrom == target.chrom and
            query.strand == target.strand)
           ):
           # both hit and target are on the same strand:
           # Note: Since the reverse complements of both the queries and the targets have already been generated for the minus strand, "+" here means that the
           # sequences are aligned in the same direction.
           if hit_line[8] == "+":
                # check if the alignment starts within 20 nucleotides from the end of the offspring (hit_line[11]) 
                # and within 30 nucleotides on the target (hit_line[15]).
                # the 20-nucleotide accounts for potential diverse polyA matches between source and offspring.
                if (int(hit_line[11]) <= 20 and int(hit_line[15]) <= 30):  
                    # check if more than 60% of the sequence has been aligned
                    if int(hit_line[0]) >= (0.6 * int(hit_line[10])):
                        # check if the gap length in the query is shorter than 8 nt and the gap length in the target is shorter than 20 nt
                        if int(hit_line[5]) <= 8 and int(hit_line[7]) <= 20: 
                        #finally, check if the TE types of the query and the target are the same
                            if query.family_type[0:-1] == target.family_type[0:-1]:
                                line_to_append = "\t".join(_ for _ in hit_line)
                                seqs = "\t".join(str(_) for _ in target.extract_sequence(hit_line, genome_dictionary))
                                filter_results.append((line_to_append, query.query_hit_coordinate(hit_line), target.target_hit_coordinate(hit_line), seqs))
    
    except AttributeError as err:
        print(err.args[0])
        
    return filter_results


# 
# -----------------------------------------------------------------------------
# Function: Finding the first polyT within and downstream of the transduced  
# segment from the  source element.
# Logic: If there is no polyT within the transduced sequence BUT there is one
# downstream of the transduced segment, then we can be confident that the 
# source is legitimate.
# -----------------------------------------------------------------------------
# 

def find_terminator(filter_lines, strand, transD_seq, downstream_1kb_seq):
    """ finds the first polyT downsream of each transduced segment from a
    verified source element, along with its coordinates and distance to the transduced sequence's end."""

    patt_1 = r"T{4,}"


    if strand == "+":
        terminator_within_transD = set()
        terminator_within_downstreamSeq = set()
        
        for match_transd in regex.finditer(patt_1, transD_seq, regex.IGNORECASE):
            terminator_within_transD.add((match_transd.start(), match_transd.end(), match_transd.group()))
    
        for match_down in regex.finditer(patt_1, downstream_1kb_seq, regex.IGNORECASE):
            terminator_within_downstreamSeq.add((match_down.start(), match_down.end(), match_down.group()))
        
        # sort the resulting set based on the closest start of match to the Alu source
        sorted_terminator_within_transD = sorted(terminator_within_transD, key=lambda my_set: my_set[0])
        sorted_terminator_within_downstreamSeq = sorted(terminator_within_downstreamSeq, key=lambda my_set: my_set[0])
        
        return sorted_terminator_within_transD, sorted_terminator_within_downstreamSeq
            
    if strand == "-":
        terminator_within_transD = set()
        terminator_within_downstreamSeq = set()
        
        for match_transd in regex.finditer(patt_1, transD_seq, regex.IGNORECASE):
            terminator_within_transD.add((match_transd.start(), match_transd.end(), match_transd.group()))
    
        for match_down in regex.finditer(patt_1, downstream_1kb_seq, regex.IGNORECASE):
            terminator_within_downstreamSeq.add((match_down.start(), match_down.end(), match_down.group()))

        # sort the resulting set based on the closest start of match to the Alu source
        sorted_terminator_within_transD = sorted(terminator_within_transD, key=lambda my_set: my_set[0])
        sorted_terminator_within_downstreamSeq = sorted(terminator_within_downstreamSeq, key=lambda my_set: my_set[0])

        return sorted_terminator_within_transD, sorted_terminator_within_downstreamSeq