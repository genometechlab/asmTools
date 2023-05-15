import os
import sys
import pandas as pd
from collections import defaultdict
from typing import DefaultDict, TextIO, Optional
from PIL import Image, ImageDraw, ImageFont
from math import floor
import argparse


'''
This script will output a .png image showing the regions of reference that are 
covered by the assembly (or collection of reads), as well as summarizing the 
missing bases to stdout. 
Up to 4 assemblies can be simultaneously plotted. A .paf or .sam can be used. 
I chose not to use the pysam .bam parsing to prevent dependency on pysam at 
the cost of being able to use arguably the most useful alignment file.

Command Line Arguements / Usage:

    python asmToRefViz -o [outfile] \
                       -i [infile1] ... [infile4] \
                       
                       
            exclusive  -ref [path to ref]
            exclusive  -ref.fai [path to ref.fai]
                       
                       
            optional   -color 128,128,128,256 ... 0,64,256,256 \
            optional   -size 3600 9000
    
'''

def create_fai(ref_path: str) -> bool:
    try:
        os.system(f"samtools faidx {ref_path}")
        return True
    
    except Exception as e:
        print(e)
        return False

def create_chrom_list(ref_path_fai: str,
                      contig_size_cutoff: int
                     ) -> list[Optional[tuple[str,int]]]:
    chrom_list: list[Optional[tuple[str, int]]] = []
    with open(ref_path_fai) as f:
        for line in f:
            if int(line.split()[1]) >= contig_size_cutoff:
                chrom_list.append((line.split()[0], int(line.split()[1])))
    return sorted(chrom_list, key=lambda x: x[1], reverse=True)
        

def parse_cigar_for_length(cigar: str) -> int:
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

def collapse_list(l: list[list[int]]) -> list[list[int]]:
    '''
    Takes a list of coordinate pairs [[int(start1), int(stop1)],
                                      [int(start2), int(stop2)]]
    and collapses overlapping regions into single contigous pairs.
    Returns a list of coordinate pairs.
    '''

    #Final output list
    total_coverage: list[Optional[list[int]]] = []
    
    #High water mark for tracking overlapping coverage
    current_best: list[int] = l[0]

    #Iterate through all pairs of coordinates in input list
    pair: list[int]
    for pair in l[1:]:
        
        #If the current pair does not extend the current best alignment
        #And is starting before the end of the current best alignment
        if pair[0] <= current_best[1] and pair[1] <= current_best[1]:
            continue
            
        #If the current pair extends the best pair, update best pair
        if pair[0] <= current_best[1] and current_best[1] <= pair[1]:
            current_best = [current_best[0], pair[1]]
        
        #If first coordinate of pair is greater than the ending coordinate
        #of the current best, add the current best to the collapsed alignment
        #output and restart with the current pair
        else:
            total_coverage.append(current_best)
            current_best = pair
    #Add the final updated alignment pair to the output
    total_coverage.append(current_best)
    
    return total_coverage

def make_chrom_dict(file_path: str) -> DefaultDict[str, 
                                                   list[Optional[list[int]]]]:
    '''
    Read in a paf file, create a list of alignment coordinates for each 
    chromosome.
    Return a dictionary with chromosomes as keys and a list of alignment
    coordinate pairings
    '''
    
    #Defaultdict for final output, keys are chromosomes, lists contain 
    #alignment pairs
    chrom_dict: DefaultDict[str, list[Optional[list[int]]]] = defaultdict(list)
    
    f: TextIO
    
    with open(file_path, 'r') as f:
        line: str
        line = f.readline()
        
        #check if file is a .paf file
        if len(line.split()) > 4 and (line.split()[4] == '+' or 
                                      line.split()[4] == '-'):
            if 'tp:A:S' in line.split()[12]:
                pass
            else:
                chrom_dict[line.split()[5]].append([int(line.split()[7]), 
                                                    int(line.split()[8])])
                
            for line in f:
                #Ignore secondary alignments, we're assuming that there is a 1:1
                #ground truth alignment from the assembly to chm13
                if 'tp:A:S' in line.split()[12]:
                    continue
            
                #Append each line to the chromdict coordinates list
                chrom_dict[line.split()[5]].append([int(line.split()[7]), 
                                                    int(line.split()[8])])
        
        #Check if file is a simple .bed file
        elif line[0] != '@' and len(line.split('\t')) == 3:
            chrom_dict[line.split('\t')[0]].append([int(line.split('\t')[1]),
                                                    int(line.split('\t')[2])])
            for line in f:
                chrom_dict[
                    line.split('\t')[0]].append([int(line.split('\t')[1]),
                                                 int(line.split('\t')[2])])                
        
        #Else for sam file
        else:
            #bitwise check if secondary alignment
            if line[0] == '@' or int(line.split()[1]) & (1 << 8):
                pass
            else:
                align_len = (int(line.split()[3]) + 
                             parse_cigar_for_length(line.split()[5]) - 1)
                chrom_dict[line.split()[2]].append([int(line.split()[3]) - 1,
                                                    align_len])
            for line in f:
                if line[0] == '@' or int(line.split()[1]) & (1 << 8):
                    continue
                align_len = (int(line.split()[3]) + 
                             parse_cigar_for_length(line.split()[5]) - 1)
                chrom_dict[line.split()[2]].append([int(line.split()[3]) - 1,
                                                    align_len])
    
    #For each chromosome sort the list by the first coordinate and collapse
    #the coverage of overlapping segments
    chrom: str
    for chrom in list(chrom_dict.keys()):
        chrom_dict[chrom] = sorted(chrom_dict[chrom], 
                                   key=lambda x: x[0])
        chrom_dict[chrom] = collapse_list(chrom_dict[chrom])
    
    return chrom_dict

def output_total_missing_bases(asm_name: str , 
                               chrom_dict:\
                               DefaultDict[str,list[Optional[list[int]]]],
                               chrom_list: list[Optional[tuple[str, int]]]
                              )-> None:
    total_missing_bases: int = 0
    chrom: str
    
    print(f"Alignment Summary for {asm_name}")
    #Using max length of chm13 cacluate missing bases by subtracting
    #each aligned pairing from collapsed list
    for chrom in chrom_list:
        total_length: int = chrom[1]
        pair: list[int]
        for pair in chrom_dict[chrom[0]]:
            total_length -= (pair[1] - pair[0])
        total_missing_bases += total_length
        
        #Formatting string
        length_string: str = f"{chrom[1] - total_length}/{chrom[1]}"
        
        #Output missing bases for each chrom to stdout
        print(f"\tFor {chrom[0]:>5} {length_string:>19} bases" +
              f" were aligned for a total of " +
              f"{round((chrom[1] - total_length)/chrom[1] * 100, 2)}%")   
    print(f"\tTotal missing bases for {asm_name} : {total_missing_bases}")

def multiple_assembly(input_dicts: \
                      list[DefaultDict[str,list[Optional[list[int]]]]],
                      dict_order_size: list[tuple[str,int]],
                      size: tuple[[int, int]] = None,
                      chrom_size_adjustment: int =  250000000) -> \
                      tuple[list[Optional[tuple[list[float], int]]],
                            list[Optional[list[float]]],
                            list[Optional[list[float]]],
                            list[Optional[list[float]]]]:

    #set default size
    if size is None:
        size = [3600, 9000]
    xsize: int
    ysize: int
    
    xsize = int(size[1])
    ysize = int(size[0])
    
    vertical_margin: float = 0.001
    
    #Determine the total number of zones needed
    num_zones: int = 0
    
    plot_zones: set[str] = set([x[0] for x in dict_order_size])
    input_dict : DefaultDict[str,list[Optional[list[int]]]]
    for input_dict in input_dicts:
        for key in input_dict:
            if key in plot_zones:
                num_zones += 1
    dict_count: int = len(input_dicts)
    
    #Determine how much vertical room each chromosome will get on the canvas
    zone_height: float = ysize * (1-vertical_margin) / num_zones
    zone_buffer: float = (ysize / num_zones) * 0.2
    
    #Calculate size per nucleotide
    xbuffer: float = 0.01
    per_nucleotide: float = ((1 - 2 * xbuffer) * xsize) / chrom_size_adjustment
    xstart = xbuffer * xsize
    
    #Construct Frames and Text
    frames: list[Optional[tuple[list[float], int]]] = []
    text: list[Optional[list[float]]] = []
    start_end: list[Optional[list[float]]] = []
    legend: list[Optional[list[float]]] = [None] * len(input_dicts)
    

    i: int
    j: int
    key_pair: tuple[str, int]
    pair: list[int]
    
    #Iterate through chromosomes in order
    for i, key_pair in enumerate(dict_order_size):
        
        #Iterate through each asm dictionary first determining text box,
        #then start and stop, then iterating through key pairs
        for j, input_dict in enumerate(input_dicts):
            key: str = key_pair[0]
            x1: float = 0.0
            x2: float = xbuffer * xsize * 0.9
            y2: float = zone_height * (i * dict_count + j) + (vertical_margin * 
                                                              ysize + 
                                                              zone_buffer)
            y1: float = y2 + zone_height - zone_buffer
            text.append([x1, y1, x2, y2])
            
            x1 = xbuffer * xsize * 0.9
            x2 = xbuffer * xsize
            start_end.append([x1, y1, x2, y2])
            
            x1 = xstart + key_pair[1] * per_nucleotide
            x2 = x1 + xbuffer * xsize * 0.1
            start_end.append([x1, y1, x2, y2])
            
            #Iterate through each pair in the assocatied dictionary
            for pair in input_dict[key]:
                x1 = xstart + per_nucleotide * pair[0]
                x2 = xstart + per_nucleotide * pair[1]
                frames.append(([x1, y1, x2, y2], j))
            legend[j] = (float(xsize) - 0.2 * float(xsize), 
                            frames[-1][0][1], 
                            float(xsize) - 0.11 * float(xsize), 
                            frames[-1][0][3])
    
    return (frames, text, start_end, legend)

def draw_image(f: list[Optional[tuple[list[float], int]]], 
               t: list[Optional[list[float]]], 
               start_end: list[Optional[list[float]]],
               legend: list[Optional[list[float]]],
               outpath: str,
               asm_list: list[str],
               chrom_list: list[Optional[tuple[str, int]]],
               color: Optional[list[list[int]]] = None,
               size: Optional[list[int]] = None,
               fontsize: Optional[int] = None
              ) -> None:

    '''
    Draw image takes the list of text boxes, start and stop coordinates,
    and alignment boxes and using PIL draws rectangles representing these 
    features. Gaps in the rectangles are areas where no alignment was found.
    '''
    if fontsize is None:
        fontsize = 15
    font = ImageFont.truetype("Gidole-Regular.ttf", 
                              fontsize)
    if size is None:
        size = [3600, 9000]
    xsize: int
    ysize: int
    
    xsize = int(size[1])
    ysize = int(size[0])
    
    #Check if the number of colors
    if color is not None and len(color) != len(asm_list):
        print(f"Warning: Assembly Count ({len(asm_list)}) and Color Count " + 
              f"({len(color)}) should match for best results.")
    
    #Default color dictionary
    color_dict: dict[Optional[str, tuple[int, int, int, int]]]
    if color is None:
        color_dict = {0:(116, 34, 128, 256),
                      1:(178, 125, 28, 256),
                      2:(12, 122, 17, 256),
                      3:(50, 97, 169, 256)
                  }
    else:
        color_dict = {}
        i: int
        c: tuple[int, int, int, int]
        for i, c in enumerate(color):
            color_dict[i] = tuple([int(x) for x in c.split(',')])

    canvas = Image.new('RGBA', 
                       (xsize, ysize), 
                       (256, 256, 256, 256))
    draw = ImageDraw.Draw(canvas)
    
    for i, text_box in enumerate(t):
        txt = f"{chrom_list[floor(i/len(asm_list))][0]}"
        draw.rectangle(text_box, 
                       fill = (128, 128, 128, 64))
        

        draw.text((text_box[0], 
                   text_box[3]), 
                  txt, 
                  font=font,
                  fill=(0, 0, 0, 256)
                 )
    
    for start_stop in start_end:
        draw.rectangle(start_stop, 
                       fill = (184, 79, 72, 128)
                      )
    
    for i, leg_box in enumerate(legend):
        txt = f"{asm_list[i]}"
        draw.rectangle(leg_box,
                       fill = color_dict[i]
                      )
        draw.text((leg_box[0] + .1 * float(xsize),
                   leg_box[3]),
                   txt,
                   font=font,
                   fill=(0,0,0,256)
                  )
                       
    
    frame: tuple[list[float], int]
    for frame in f:
        asm_color: tuple[int, int, int, int] = color_dict[frame[1]]
        draw.rectangle(frame[0], 
                       fill = asm_color)
    canvas.save(outpath, 
                "PNG")
         
def main(outpath: str, 
         infile_list: list[str],
         color: Optional[list[list[str]]],
         size: Optional[list[str]],
         contig_size_cutoff: Optional[int],
         fontsize: Optional[int],
         ref_fai: str
        ) -> None:
    '''
    Driver function, taking the path for the output file, and
    a list of input files. No return value.
    '''
    
    if not contig_size_cutoff:
        contig_size_cutoff = 0
    
    chrom_list: list[tuple[str, int]] = create_chrom_list(ref_fai,
                                                          contig_size_cutoff
                                                         )
    chrom_dict_list: list[DefaultDict[str,list[Optional[list[int]]]]] = \
                               [make_chrom_dict(x) for x in infile_list]
    name: str
    chrom_dict: DefaultDict[str,list[Optional[list[int]]]]
    for name, chrom_dict in zip(infile_list, chrom_dict_list):
        name = name.split('/')[-1]
        output_total_missing_bases(name, chrom_dict, chrom_list)

    f: list[Optional[tuple[list[float], int]]] 
    t: list[Optional[list[float]]] 
    s: list[Optional[list[float]]] 
    f,t,s, l = multiple_assembly(chrom_dict_list,
                                 chrom_list,
                                 size=size,
                                 chrom_size_adjustment=(1.01 * max(
                                     [x[1] for x in chrom_list])))
    draw_image(f,
               t,
               s,
               l,
               outpath,
               [x.split('/')[-1] for x in infile_list],
               chrom_list,
               color,
               size,
               fontsize
              )
    

    
if __name__ == "__main__":
    
    # Parse arguements
    
    parser = argparse.ArgumentParser()
    
    parser.add_argument(
        "--outpath",
        "-o",
        type=str,
        required=True,
        help = "Path for saving .png file"
    )
    
    parser.add_argument(
        "--inpath",
        "-i",
        type=str,
        required=True,
        nargs='+',
        help = "File paths for asm comparisson"
    )
                              
    arg_group1 = parser.add_mutually_exclusive_group(
        required=True
    )
                              
    arg_group1.add_argument(
        "--ref",
        type=str,
        help="Path to reference file, this will create an .fai file"
    )
    
    arg_group1.add_argument(
        "--ref_fai",
        type=str,
        help="Path to .fai file for reference"
    )
    
    parser.add_argument(
        "--contig_cutoff",
        type=int,
        default=None,
        help="Filter out all contigs smaller than this value"
    )
        
    parser.add_argument(
        "--color",
        "-c",
        type = str,
        nargs='+',
        help = 'RGBA color values seperated by , and a space between' +
        'colors ex. 156,156,128,256 212,10,15,256\n The number of colors '+
        'should match the number of assemblies'
    )
    
    parser.add_argument(
        "--size",
        "-s",
        type=str,
        nargs=2,
        help = 'xsize  ysize'
    )
    
    parser.add_argument(
        "--fontsize",
        "-f",
        type=int,
        help = "Font size, the larger the number of contigs the smaller you " +
        "should go."
    )
        
    FLAGS, unparsed = parser.parse_known_args()
    
    if FLAGS.ref:
        if create_fai(FLAGS.ref):
            ref_fai = f"{FLAGS.ref}.fai"
        else:
            print("Could not produce index file for reference")
            sys.exit()
    else:
        ref_fai = FLAGS.ref_fai
                            
                              
    main(
        FLAGS.outpath, 
        FLAGS.inpath,
        FLAGS.color,
        FLAGS.size,
        FLAGS.contig_cutoff,
        FLAGS.fontsize,
        ref_fai
        )