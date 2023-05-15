# Tools for Assembly Analysis
## Contents
- asmToRefViz.py
- Contigous_alignment.py
### asmToRefViz.py
This script produces a visualization of alignments from a sam, bed, or paf file to the associated reference.<br>
Input: alignment file (.sam, .bed, .paf), reference file (.fa or .fai)<br>
Output: .png file showing alignments<br>

usage: asmToRefViz.py [-h] --outpath OUTPATH --inpath INPATH [INPATH ...] (--ref REF | --ref_fai REF_FAI) [--contig_cutoff CONTIG_CUTOFF] [--color COLOR [COLOR ...]] [--size SIZE SIZE] [--fontsize FONTSIZE]<br>
<br>
optional arguments:<br>
  -h, --help<br>            &emsp;&emsp;&emsp;show this help message and exit<br>
  --outpath OUTPATH, -o OUTPATH<br>
                        &emsp;&emsp;&emsp;Path for saving .png file<br>
  --inpath INPATH [INPATH ...], -i INPATH [INPATH ...]<br>
                        &emsp;&emsp;&emsp;File paths for asm comparisson<br>
  --ref REF <br>             &emsp;&emsp;&emsp;Path to reference file, this will create an .fai file<br>
  --ref_fai REF_FAI <br>     &emsp;&emsp;&emsp;Path to .fai file for reference<br>
  --contig_cutoff CONTIG_CUTOFF<br>
                        &emsp;&emsp;&emsp;Filter out all contigs smaller than this value<br>
  --color COLOR [COLOR ...], -c COLOR [COLOR ...]<br>
                        &emsp;&emsp;&emsp;RGBA color values seperated by , and a space between colors ex. 156,156,128,256 212,10,15,256 The number of colors should <br>&emsp;&emsp;&emsp;match the number of assemblies<br>
  --size SIZE SIZE, -s SIZE SIZE<br>
                        &emsp;&emsp;&emsp;xsize ysize<br>
  --fontsize FONTSIZE, -f FONTSIZE<br>
                        &emsp;&emsp;&emsp;Font size, the larger the number of contigs the smaller you should go.<br>


### Contigous_alignments.py
This script filters reads based on a variety of parameterized contiguity measures<br>
Input: alignment file (.sam)<br>
Output: filtered alignment file (.sam)<br>

usage: Contigous_alignment.py [-h] --output OUTPUT [--supplementary_output SUPPLEMENTARY_OUTPUT] --input INPUT [INPUT ...] [--min_coverage MIN_COVERAGE] [--max_insertion MAX_INSERTION] [--max_deletion MAX_DELETION]
                              [--max_softclip MAX_SOFTCLIP] [--max_hardclip MAX_HARDCLIP]<br>
<br>
optional arguments:<br>
  -h, --help<br>            &emsp;&emsp;&emsp;show this help message and exit<br>
  --output OUTPUT, -o OUTPUT<br>
                        &emsp;&emsp;&emsp;Path to desitnation output bed file<br>
  --supplementary_output SUPPLEMENTARY_OUTPUT, -so SUPPLEMENTARY_OUTPUT<br>
                        &emsp;&emsp;&emsp;Optional output for regions deemed not covered<br>
  --input INPUT [INPUT ...], -i INPUT [INPUT ...]<br>
                        &emsp;&emsp;&emsp;Any number of input sam paths, must be for same reference<br>
  --min_coverage MIN_COVERAGE, -m MIN_COVERAGE<br>
                        &emsp;&emsp;&emsp;The minimum number of contigous reads to count as observed.<br>
  --max_insertion MAX_INSERTION<br>
                        &emsp;&emsp;&emsp;The maximum size of a single insertion event in a given SAM alignment.<br>
  --max_deletion MAX_DELETION<br>
                        &emsp;&emsp;&emsp;The maximum size of a single deletion event in a given SAM alignment.<br>
  --max_softclip MAX_SOFTCLIP<br>
                        &emsp;&emsp;&emsp;Maximum allowable size for softclip<br>
  --max_hardclip MAX_HARDCLIP<br>
                        &emsp;&emsp;&emsp;Maximum allowable size for hardclip<br>
