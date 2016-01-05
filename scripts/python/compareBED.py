#!/usr/local/bin/python
"""
***********************************************
- PROGRAM: compareBED.py
- CONTACT: Bryan lajoie (bryan.lajoie@umassmed.edu)
***********************************************
"""

  
from __future__ import print_function
from __future__ import division

import numpy as np
import scipy as sp
import pdb
import h5py
import sys
import argparse
import logging
import time
import gzip
import re
import os
import math
import uuid
import socket
from itertools import izip
from collections import defaultdict
from datetime import datetime

import matplotlib.pyplot as plt

__version__ = "1.00"

def main():
    
    parser=argparse.ArgumentParser(description='Compare BED files',formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('-1', '--input_1', dest='input_1', type=str, required=True, help='input bed file 1')
    parser.add_argument('-2', '--input_2', dest='input_2', type=str, required=True, help='input bed file 2')
    parser.add_argument('--cm', '--compareMode', dest='compare_mode', type=str, required=True, help='bed compare function [subtract/log2ratio/add/multiple/min/max]')
    parser.add_argument('--ya','--yaxisrange', dest='yaxisrange', type=float, nargs='+', required=False, default=[-1,-1], help='y-axis for aggregrate plot')
    parser.add_argument('-v', '--verbose', dest='verbose', action='count', help='Increase verbosity (specify multiple times for more)')
    parser.add_argument('--version', action='version', version='%(prog)s 1.0')
    
    args=parser.parse_args()

    input_1=args.input_1
    input_2=args.input_2
    compare_mode=args.compare_mode
    yaxisrange=args.yaxisrange
    verbose=args.verbose
    
    log_level = logging.WARNING
    if verbose == 1:
        log_level = logging.INFO
    elif verbose >= 2:
        log_level = logging.DEBUG
    logging.basicConfig(level=log_level)
    
    verboseprint = print if verbose else lambda *a, **k: None
    
    if not os.path.isfile(input_1):
        sys.exit('invalid input file! (non-existant)')
    if not os.path.isfile(input_2):
        sys.exit('invalid input file! (non-existant)')
    
    print(input_1)
    print(input_2)
    
    bin_midpoints=list()
    bin_scores=list()
    
    with input_wrapper(input_1) as fh1, input_wrapper(input_2) as fh2: 
        for i,(line_1,line_2) in enumerate(izip(fh1,fh2)):
            line_1=line_1.rstrip("\n")
            if line_1.startswith("#") or line_1.startswith("track"):
                continue
            line_2=line_2.rstrip("\n")
            if line_2.startswith("#") or line_2.startswith("track"):
                continue
                
            line_1_list=line_1.split("\t")
            line_1_chr=line_1_list[0]
            line_1_start=int(line_1_list[1])
            line_1_end=int(line_1_list[2])
            if(len(line_1_list) == 4):
                line_1_score=float(line_1_list[3])
                line_1_name="NA"
            elif(len(line_1_list) > 4):
                line_1_name=line_1_list[3]
                line_1_score=float(line_1_list[4])
            
            line_2_list=line_2.split("\t")
            line_2_chr=line_2_list[0]
            line_2_start=int(line_2_list[1])
            line_2_end=int(line_2_list[2])
            if(len(line_2_list) == 4):
                line_2_name="NA"
                line_2_score=float(line_2_list[3])
            elif(len(line_2_list) > 4):
                line_2_name=line_2_list[3]
                line_2_score=float(line_2_list[4])
            
            if(line_1_chr != line_2_chr):
                sys.exit('error: bed file not aligned/equal\n\t'+'line_num:'+str(i)+"\t"+line_1_chr+' != '+line_2_chr+'\n')
            if(line_1_start != line_2_start):
                sys.exit('error: bed file not aligned/equal\n\t'+'line_num:'+str(i)+"\t"+line_1_start+' != '+line_2_start+'\n')
            if(line_1_end != line_2_end):
                sys.exit('error: bed file not aligned/equal\n\t'+'line_num:'+str(i)+"\t"+line_1_end+' != '+line_2_end+'\n')
            if(line_1_name != line_2_name):
                sys.exit('error: bed file not aligned/equal\n\t'+'line_num:'+str(i)+"\t"+line_1_name+' != '+line_2_name+'\n')
             
            chr=line_1_chr=line_2_chr
            start=line_1_start=line_2_start
            end=line_1_end=line_2_end
            name=line_1_name=line_2_name
            
            compare_score="NA"
            if compare_mode == "subtract":
                compare_score=line_1_score-line_2_score
            if compare_mode == "min":
                compare_score=min(line_1_score,line_2_score)
            if compare_mode == "max":
                compare_score=max(line_1_score,line_2_score)
            if compare_mode == "add":
                compare_score=line_1_score+line_2_score
            if compare_mode == "divide":
                compare_score="NA"
                if line_2_score != 0:
                    compare_score=line_1_score/line_2_score
            if compare_mode == "log2ratio":
                if line_1_score != 0 and line_2_score != 0:
                    compare_score=math.log((line_1_score/line_2_score),2)
            
            print(chr,start,end,name,line_1_score,line_2_score,compare_score)
            
            if compare_score != "NA":
                bin_midpoints.append((start+end)/2)
                bin_scores.append(compare_score)
            
    if((len(yaxisrange) != 2) or (yaxisrange[0] == yaxisrange[1])):
        yaxisrange=[min(bin_scores),max(bin_scores)]
        
    plt.figure(figsize=(20,10))
    plt.plot(bin_midpoints,bin_scores)
    plt.xlabel('genome distance from element-anchor')
    plt.ylabel('aggregrate signal')
    plt.title("bed comparison")
    plt.grid(True)
    axisbounds=[min(bin_midpoints),max(bin_midpoints)]+yaxisrange
    plt.axis(axisbounds)
    plt.savefig('test.aggregrate.vector.png')
            
    print("")

def input_wrapper(infile):
    if infile.endswith('.gz'):
        fh=gzip.open(infile,'r')
    else:
        fh=open(infile,'r')
        
    return fh
    
def output_wrapper(outfile,append=False,suppress_comments=False):
    
    if outfile.endswith('.gz'):
        if append:
            fh=gzip.open(outfile,'a')
        else:
            fh=gzip.open(outfile,'w')   
    else:
        if append:
            fh=open(outfile,'a')
        else:
            fh=open(outfile,'w')
    
    # disable comment(s)if (UCSC format file)
    if outfile.endswith('.bed'):
        suppress_comments = True
    if outfile.endswith('.bed.gz'):
        suppress_comments = True
    if outfile.endswith('.bedGraph'):
        suppress_comments = True
    if outfile.endswith('.bedGraph.gz'):
        suppress_comments = True
    if outfile.endswith('.wig'):
        suppress_comments = True
    if outfile.endswith('.wig.gz'):
        suppress_comments = True
    if outfile.endswith('.sam'):
        suppress_comments = True
    if outfile.endswith('.sam.gz'):
        suppress_comments = True
    if outfile.endswith('.bam'):
        suppress_comments = True
    if outfile.endswith('.fastq'):
        suppress_comments = True
    if outfile.endswith('.fastq.gz'):
        suppress_comments = True

    if not suppress_comments:
        print("## ",os.path.basename(__file__),sep="",file=fh)
        print("## ",sep="",file=fh)
        print("## Dekker Lab",sep="",file=fh)
        print("## Contact: Bryan R. Lajoie",sep="",file=fh)
        print("## https://github.com/blajoie",sep="",file=fh)
        print("## ",sep="",file=fh)
        print("## Version:\t",__version__,sep="",file=fh)
        print("## Date:\t",get_date(),sep="",file=fh)
        print("## Host:\t",get_compute_resource(),sep="",file=fh)
    
    return(fh)

def get_date():
    time=datetime.now()
    date=time.strftime('%I:%M:%S %p, %m/%d/%Y')
    
    return date

def get_compute_resource():
    return(socket.gethostname())
    
if __name__=="__main__":
      main()
