#!/usr/local/bin/python
"""
***********************************************
- PROGRAM: boundary2tad.py
- CONTACT: Bryan lajoie (bryan.lajoie@umassmed.edu)
***********************************************
"""

from __future__ import print_function

# Built in modules
import argparse
import os.path
import sys
import gzip
import re
import itertools
import time
import gzip
from operator import itemgetter

from scipy.stats.stats import nanmean

import numpy as np
import scipy as sp
import matplotlib
matplotlib.use('pdf')
import matplotlib.pyplot as plt

# For eigenvectors and eigenvalues
from scipy import linalg as la
from math import cos,log,sin,sqrt 
from scipy import weave 

def main():
    print("")
    
    # Get input options
    args = check_options()

    # Store the variables
    boundaryFile = args.boundaryFile
    insulationFile = args.insulationFile
    
    if not os.path.isfile(boundaryFile):
        sys.exit('invalid input file! (non-existant)')
    if not os.path.isfile(insulationFile):
        sys.exit('invalid input file! (non-existant)')
    
    jobName=os.path.basename(boundaryFile)
    jobName=re.sub(".gz", "", jobName)    
    jobName=re.sub(".matrix", "", jobName)
    
    print("")
    
    nBoundaries,chr2index,index2chr=validate_boundary_file(boundaryFile)
    
    #print("loading boundaries ... ",end="")
    if boundaryFile.endswith('.gz'):
        b_fh=gzip.open(boundaryFile,'r')
    else:
        b_fh=open(boundaryFile,'r')
    
    field=b_fh.readline().rstrip("\n").split("\t")
    index2field=dict(enumerate(field))
    field2index=dict((value, key) for key, value in index2field.iteritems())
        
    # create empty numpy structured array
    boundary_ref=np.empty(nBoundaries, 
        dtype={'names':['index', 'header', 'boundaryHeader', 'start', 'end', 'boundaryStrength','chr','available'],
        'formats':['int64','a500','a500','int64','int64','float64','int64','bool']})
    
    # load boundaries
    
    for i,line in enumerate(b_fh):
        l=line.rstrip("\n").split("\t")
        
        header=l[field2index["header"]]
        header_name,header_assembly,header_chr,header_start,header_end=splitHeader(header)
        
        chr_id=chr2index[header_chr]
        boundary_ref[i]=(i,l[field2index["header"]],l[field2index["boundaryHeader"]],l[field2index["start"]],l[field2index["end"]],l[field2index["boundaryStrength"]],chr_id,1)
        
    b_fh.close()
    
    nchrs=len(chr2index)
    print("nChrs=",nchrs)
    chr_range=np.zeros((nchrs,2),dtype='int32')
    for chr_idx in xrange(nchrs):
        chr=index2chr[chr_idx]
        chr_boundaries=np.nonzero(boundary_ref["chr"]==chr_idx)
        chr_range[chr_idx]=np.min(chr_boundaries),np.max(chr_boundaries)
        
    # shift all field2index keys (due to adding in i=index)
    for key in field2index:
        value=field2index[key]
        field2index[key]=value+1
    field2index["index"]=0
           
    print("")
    
    tads=[]
    
    print("assmebing tads") 
    for chr_idx in xrange(nchrs):
        chr=index2chr[chr_idx]
        chr_start,chr_end=chr_range[chr_idx]
        #print(chr_idx,chr,chr_range[chr_idx])
        
        chr_ref=boundary_ref[chr_start:chr_end]
        chr_boundaries=chr_ref[["boundaryStrength","index"]]
        chr_boundaries=chr_boundaries.tolist()
        chr_boundaries=sorted(chr_boundaries,key=lambda chr_boundaries:chr_boundaries[0], reverse=True)
        
        boundary_ref,chr_boundaries,chr_range,chr_idx,tads,field2index=assemble_tads(boundary_ref,chr_boundaries,chr_range,chr_idx,tads,field2index)
        
        tad_bedFile=jobName+"__tads.bed"
        tad_fh=open(tad_bedFile,"w")
        print("track name='tad_bed' description='tad_bed' visibility=dense",end="\n",file=tad_fh)

        for i in xrange(len(tads)):
            tad=tads[i]
            
            chr_1=splitHeader(tad[0])[2]
            chr_2=splitHeader(tad[0])[2]
            
            if chr_1 != chr_2:
                sys.exit('error - inter TAD detected')
            chr=chr_1=chr_2
            chr=deGroupChr(chr)
            name="___".join(tad)
            tad_start=splitHeader(tad[0])[3]
            tad_end=splitHeader(tad[1])[4]
            
            print(chr,tad_start,tad_end,name,1000,file=tad_fh)
        
        tad_fh.close()
    print("")

def deGroupChr(chr_id):
    return(chr_id.split('-')[0])
    
def assemble_tads(boundary_ref,chr_boundaries,chr_range,chr_idx,tads,field2index):
    
    chr_start,chr_end=chr_range[chr_idx]
    
    if len(chr_boundaries) == 0:
        return boundary_ref,chr_boundaries,chr_range,chr_idx,tads,field2index
    else:
        tmp_strength,tmp_index=chr_boundaries.pop()
        boundary_ref[tmp_index]["available"]=0
        
        tads=create_tad(boundary_ref,tmp_index,chr_range[chr_idx],tads)
        
        return assemble_tads(boundary_ref,chr_boundaries,chr_range,chr_idx,tads,field2index)

def create_tad(boundary_ref,anchor_idx,chr_bound,tads):
    
    strength_variance=0.25
    
    left_bound=None
    left_ref=np.nonzero( 
        (boundary_ref["available"]==True) & 
        ((boundary_ref["index"]>=chr_bound[0]) & (boundary_ref["index"]<=chr_bound[1])) & 
        (boundary_ref["index"]<anchor_idx) )
    left_size=np.shape(left_ref)[1]
    if left_size > 0:
        left_bound=np.max(left_ref)
        
    right_bound=None
    right_ref=np.nonzero( 
        (boundary_ref["available"]==True) & 
        ((boundary_ref["index"]>=chr_bound[0]) & (boundary_ref["index"]<=chr_bound[1])) & 
        (boundary_ref["index"]>anchor_idx) )
    right_size=np.shape(right_ref)[1]
    if right_size > 0:
        right_bound=np.min(right_ref)
    
    anchor_header=boundary_ref[anchor_idx]["header"]
    anchor_strength=boundary_ref[anchor_idx]["boundaryStrength"]
    #print("ANCHOR",anchor_idx,anchor_header,anchor_strength)
    #print(left_bound,anchor_idx,right_bound)
    
    left_idx=None
    if left_bound != None:
        for i in xrange(anchor_idx-1,left_bound-1,-1):
            tmp_strength=boundary_ref[i]["boundaryStrength"]
            #print("\t",anchor_idx,"-",i,"\t",anchor_strength,tmp_strength)
            if(abs(anchor_idx-left_bound) > 1):
                strength_diff=abs(anchor_strength-tmp_strength)
                #print("\t\t",strength_diff)
                if(strength_diff < strength_variance):
                    break
            else:
                left_idx=i
                
        if left_idx != None:
            left_tad=[boundary_ref[left_idx]["header"],boundary_ref[anchor_idx]["header"]]
            tads.append(left_tad)
            #print(left_tad)
            #print("")
        
    if right_bound != None:
        right_tad=[boundary_ref[anchor_idx]["header"],boundary_ref[right_bound]["header"]]
        tads.append(right_tad)
        #print(right_tad)
    
    return(tads)
    
def splitHeader(header):
    m=re.search(r'(\S+)\|(\S+)\|(\S+):(\d+)-(\d+)',header)
    if m==None:
        sys.exit('error: incorrect header format ['+str(header)+']!')  

    header_name,header_assembly,header_chr,header_start,header_end=m.groups()
    
    return(header_name,header_assembly,header_chr,header_start,header_end)
        
def validate_boundary_file(boundaryFile):

    if boundaryFile.endswith('.gz'):
        b_fh=gzip.open(boundaryFile,'r')
    else:
        b_fh=open(boundaryFile,'r')
    
    field=b_fh.readline().rstrip("\n").split("\t")
        
    index2field=dict(enumerate(field))
    field2index=dict((value, key) for key, value in index2field.iteritems())
    
    current_chr=0
    chr2index=dict()
    index2chr=dict()
    for i, line in enumerate(b_fh):
        l=line.rstrip("\n").split("\t")
        header=l[field2index["header"]]
        
        header_name,header_assembly,header_chr,header_start,header_end=splitHeader(header)
        
        if(header_chr not in chr2index):
            chr2index[header_chr]=current_chr
            index2chr[current_chr]=header_chr
            current_chr+=1
            
        chr_id=chr2index[header_chr]
        if(chr_id < (current_chr-1)):
            sys.exit('improperly sorted boundary file!')
    
    b_fh.close()
    
    nLines=i+1
    
    return(nLines,chr2index,index2chr)
    
def check_options():
    ''' Checks the options to the program '''

    # Create parser object
    parser = argparse.ArgumentParser(add_help=True,
        formatter_class=argparse.RawDescriptionHelpFormatter)

    # Add arguments 
    parser.add_argument('-bf' , metavar='--boundaryFile'  , help="boundary input file", dest="boundaryFile", type=str, default="")
    parser.add_argument('-if' , metavar='--insulationFile'  , help="insulaiton input file", dest="insulationFile", type=str, default="")
    
    # Parse command line with parse_args and store it in an object
    args = parser.parse_args()

    return args

if __name__=="__main__":
      main()
