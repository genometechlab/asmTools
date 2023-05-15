#!/usr/bin/python3

import argparse
import os
import time
from typing import Optional, DefaultDict, TextIO
import numpy as np

read_dict_type = dict[Optional[list[Optional[tuple[int, int]]]]]

def parse_cigar(cig: str) -> list[tuple[int,str]]:
    '''
    Return a list of cigar elements as tuples
    '''
    
    #print(f"Parse_Cigar : {cigar}")
    cigar_components: list[Optional[tuple[int,str]]] = []
    current_num: str = ""
    ch: str
    for ch in cig:
        
        #split, add, reset on identified alpha char (I, D, M, H, S)
        if ch.isalpha():
            cigar_components.append((int(current_num), ch))
            current_num = ""
        else:
            current_num += ch
    return cigar_components

def parse_cigar_for_length(cigar: str) -> int:
    #Determine the aligned length with a cigar string
    count:str  = ""
    code: str = ""
    total_dist: int = 0
    c: str
    for c in cigar:
        if c.isalpha():
            if c == 'M' or c == 'D':
                total_dist += int(count)
            count = ""
        else:
            count = count + c
    return total_dist

def is_alignment_contigous(line: str,
                           max_ins: int,
                           max_del: int,
                           max_soft: Optional[int] = None,
                           max_hard: Optional[int] = None
                          ) -> tuple[bool, str]:
    #Check if an alignment is contigous by parsing cigar string elements
    #print(f"Is alignment_contigous\n{line}")
    for element in parse_cigar(line):
        if element[1] == 'I' and element[0] > max_ins:
            return (False, 'I')
        if element[1] == 'D' and element[0] > max_del:
            return (False, 'D')
        if element[1] == 'S' and max_soft and element[0] > max_soft:
            return (False, 'S')
        if element[1] == 'H' and max_hard and element[0] > max_hard:
            return (False, 'H') 
    return (True, 'A')

def coverage_depth(alignment_list_dict: read_dict_type,
                   min_coverage_depth: int
                  ) -> read_dict_type:
    # Sort list of alignments by start position
    print("Calculating coverage depth")
    key : str
    covered_dict: read_dict_type = {}
    for key in list(alignment_list_dict.keys()):
        covered_dict[key] = list()
        print(f"{key=}")
        alignment_list_dict[key] = sorted(alignment_list_dict[key], 
                                          key=lambda x: x[0])
        pair: tuple[int, int]
        max_pos: int = 0
        for pair in alignment_list_dict[key]:
            if pair[1] > max_pos:
                max_pos = pair[1]
        coverage_list: np.array = np.zeros(max_pos+1)
        for pair in alignment_list_dict[key]:
            #trying numpy native addition instead of iteration
            coverage_list[pair[0]:pair[1] + 1] = \
            coverage_list[pair[0]:pair[1] + 1] + 1
        
        i: int
        start: int = -1
        for i in range(len(coverage_list)):
            if i % 10000000 == 0:
                print(i)
            if start == -1 and coverage_list[i] >= min_coverage_depth:
                start = i
            elif start != -1 and coverage_list[i] < min_coverage_depth:
                covered_dict[key].append((start, i-1))         
                start = -1
        if start != -1:
            covered_dict[key].append((start, i-1))
    return covered_dict
            
def write_to_bed(alignment_dict: read_dict_type,
                 output: str
                )-> None:
    
    out_fh: TextIO = open(output, 'w')
    chrom: str
    pair: tuple[int,int]
    for chrom in alignment_dict:
        for pair in alignment_dict[chrom]:
            out_fh.write(f"{chrom}\t{pair[0]}\t{pair[1]}\n")
    
def main(output: str, 
         input_list: list[str], 
         supplementary_output: Optional[str] = None, 
         max_ins: Optional[int] = None, 
         max_del: Optional[int] = None,
         max_soft: Optional[int] = None,
         max_hard: Optional[int] = None,
         min_coverage_depth: Optional[int] = None
        ):

    
    contigous_alignments: dict[Optional[list[Optional[tuple[str, int, int]]]]]\
    = {}
    
    file_path: str
    file_handle: TextIO
    for file_path in input_list:
        print(f"Reading in {file_path}")
        file_handle = open(file_path, 'r')
        
        start: int
        stop: int
        line: str
        
        i: int = 0
        contiguity_count: dict[str, int] = {'A':0, 
                                            'H':0, 
                                            'S':0, 
                                            'I':0, 
                                            'D':0}
        for line in file_handle:
            
            if i % 1000000 == 0:
                print(f"{i} reads processed")
            i += 1
            
            if line[0] == '@':
                continue

            contiguity: tuple[bool, str]
            contiguity = is_alignment_contigous(line.split()[5],
                                                max_ins,
                                                max_del,
                                                max_soft,
                                                max_hard
                                               )
            contiguity_count[contiguity[1]] += 1
            if contiguity[0]:
                chrom: str = line.split()[2]
                
                #Exclude unaligned reads
                if chrom == "*":
                    continue
                    
                start = int(line.split()[3])
                stop = (int(parse_cigar_for_length(line.split()[5])) + 
                        int(line.split()[3]))
                if chrom not in contigous_alignments:
                    contigous_alignments[chrom] = list()
                contigous_alignments[chrom].append((start,stop))
    key: str
    for key in contiguity_count:
        print(f"{contiguity_count[key]} reads were {key}")
    coverage_filtered_dict: read_dict_type = \
    coverage_depth(contigous_alignments,
                   min_coverage_depth)
    write_to_bed(coverage_filtered_dict, output)
    
                                       

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output",
        "-o",
        type=str,
        required=True,
        help="Path to desitnation output bed file"
    )
    
    parser.add_argument(
        "--supplementary_output",
        "-so",
        type=str,
        help="Optional output for regions deemed not covered",
        default=None
    )
    
    parser.add_argument(
        "--input",
        "-i",
        type=str,
        nargs='+',
        required=True,
        help="Any number of input sam paths, must be for same reference"
    )
    
    parser.add_argument(
        "--min_coverage",
        "-m",
        type=int,
        help="The minimum number of contigous reads to count as observed."
    )
    
    parser.add_argument(
        "--max_insertion",
        type=int,
        help="The maximum size of a single insertion event in a given SAM " +
        "alignment."
    )
    
    parser.add_argument(
        "--max_deletion",
        type=int,
        help="The maximum size of a single deletion event in a given SAM " +
        "alignment.",
        default=None
    )
    
    parser.add_argument(
        "--max_softclip",
        type=int,
        help="Maximum allowable size for softclip",
        default=None
    )
    
    parser.add_argument(
        "--max_hardclip",
        type=int,
        help="Maximum allowable size for hardclip",
        default=None
    )
    
    
    Flags, unparsed = parser.parse_known_args()
    
    t1:time
    t2:time
    
    #Check input types
    assert type(Flags.output) == str
    assert (not Flags.supplementary_output or 
            type(Flags.supplementary_output) == str)
    assert (type(Flags.input) == str or
            type(Flags.input) == list or
            not Flags.input)
    assert (not Flags.max_insertion or
            type(Flags.max_insertion) == int)
    assert (not Flags.max_deletion or
            type(Flags.max_deletion) == int)
    assert (not Flags.max_softclip or
            type(Flags.max_softclip) == int)
    assert (not Flags.max_hardclip or
            type(Flags.max_hardclip) == int)
    
    t1 = time.time()
    main(Flags.output, 
         Flags.input, 
         Flags.supplementary_output, 
         Flags.max_insertion, 
         Flags.max_deletion,
         Flags.max_softclip,
         Flags.max_hardclip,
         Flags.min_coverage
        )
    t2 = time.time()
    print(f"{round((t2-t1) / 60, 4)} minutes elapsed")