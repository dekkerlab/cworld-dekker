#!/usr/local/bin/python
"""
***********************************************
- PROGRAM: boundary2tad.py
- CONTACT: Bryan lajoie (bryan.lajoie@umassmed.edu)
***********************************************
"""

from __future__ import print_function
from __future__ import division

import sys
import argparse
import subprocess
import shlex
import logging
import itertools
import time
import gzip
import re
import os
import math
import uuid
import socket
from collections import defaultdict
from collections import Counter
from datetime import datetime
from operator import itemgetter

from scipy.stats.stats import nanmean

import numpy as np
import scipy as sp

# For eigenvectors and eigenvalues
from scipy import linalg as la
from math import cos,log,sin,sqrt 
from scipy import weave 

verboseprint=lambda *a, **k: None
__version__ = "1.0"
debug = None

def main():
    
    parser=argparse.ArgumentParser(description='Extract data from hdf5 file.',formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('-bf' , metavar='--boundaryFile', help="boundary input file", required=True, dest="boundaryFile", type=str, default="")
    parser.add_argument('-if' , metavar='--insulationFile', help="insulation input file", required=True, dest="insulationFile", type=str, default="")
    parser.add_argument('-bn' , metavar='--boundaryNoise', help="boundary noise estimated", default=0.25, dest="boundary_noise", type=float)
    parser.add_argument('-v', '--verbose', dest='verbose',  action='count', help='Increase verbosity (specify multiple times for more)')
    parser.add_argument('--version', action='version', version='%(prog)s '+__version__)
    
    args=parser.parse_args()

    boundaryFile = args.boundaryFile
    insulationFile = args.insulationFile
    boundary_noise = args.boundary_noise
    verbose=args.verbose
    
    log_level = logging.WARNING
    if verbose == 1:
        log_level = logging.INFO
    elif verbose >= 2:
        log_level = logging.DEBUG
    logging.basicConfig(level=log_level)
    
    global verboseprint
    verboseprint = print if verbose else lambda *a, **k: None
    
    if not os.path.isfile(boundaryFile):
        sys.exit('invalid input file! (non-existant)')
    if not os.path.isfile(insulationFile):
        sys.exit('invalid input file! (non-existant)')
    
    jobName=os.path.basename(boundaryFile)
    jobName=re.sub(".gz", "", jobName)    
    jobName=re.sub(".matrix", "", jobName)
    
    verboseprint("")
    
    # load boundary file
    nBoundaries,chr2index,index2chr,boundary_index2field,boundary_field2index,boundary_ref = load_boundary_file(boundaryFile)

    # load insulation file
    header_rows,header_cols,header2idx,insulation,nan_rowcols = load_insulation_file(insulationFile)
    
    # nchrs
    nchrs=len(chr2index)
    
    # calculate all chr ranges
    chr_range=np.zeros((nchrs,2),dtype='int32')
    for chr_idx in xrange(nchrs):
        chr=index2chr[chr_idx]
        chr_boundaries=np.nonzero(boundary_ref["chr"]==chr_idx)
        chr_range[chr_idx]=np.min(chr_boundaries),np.max(chr_boundaries)
        
    # shift all boundary_field2index keys (due to adding in i=index)
    for key in boundary_field2index:
        value=boundary_field2index[key]
        boundary_field2index[key]=value+1
    boundary_field2index["index"]=0
           
    verboseprint("")
    
    tads=[]
    
    # now assemble tads from boundaries
    #verboseprint("assmebing tads ... ",end="") 
    for chr_idx in xrange(nchrs):
        chr=index2chr[chr_idx]
        chr_start,chr_end=chr_range[chr_idx]
        
        chr_ref=boundary_ref[chr_start:chr_end]
        chr_boundaries=chr_ref[["boundaryStrength","index"]]
        chr_boundaries=chr_boundaries.tolist()
        chr_boundaries=sorted(chr_boundaries,key=lambda chr_boundaries:chr_boundaries[0], reverse=True)
        
        boundary_ref,chr_boundaries,chr_range,chr_idx,tads,boundary_field2index=assemble_tads(boundary_noise,boundary_ref,chr_boundaries,chr_range,chr_idx,tads,boundary_field2index,header2idx,insulation)
        
        tad_bedFile=jobName+"__nested-tads.bed"
        tad_fh=output_wrapper(tad_bedFile,suppress_comments=True)
        print("track name='",jobName,"__nested-tads' description='",jobName,"__nested-tads' visibility=squish",sep="",end="\n",file=tad_fh)

        for i in xrange(len(tads)):
            tad,tad_headers,tad_strength=tads[i]
            
            chr_1=splitHeader(tad[0])[2]
            chr_2=splitHeader(tad[0])[2]
            
            if chr_1 != chr_2:
                sys.exit('error - inter-chr-TAD detected, cis only please!')
                
            chr=chr_1=chr_2
            chr=deGroupChr(chr)
            name="___".join(tad)
            tad_start=splitHeader(tad[0])[3]
            tad_end=splitHeader(tad[1])[4]
            
            print(chr,tad_start,tad_end,name,tad_strength,sep="\t",file=tad_fh)
        
        tad_fh.close()
    verboseprint("done")
    
    verboseprint("")
    
    # build matrix
    rows=len(header_rows)
    cols=len(header_cols)
    verboseprint("building matrix ","[",rows,"x",cols,"]",sep="")
    matrix=np.zeros((rows,cols),dtype="float32")
    
    matrix=fillTadMatrix(header_rows,header_cols,matrix,tads,header2idx,nan_rowcols)
    
    tadMatrixFile=jobName+'__nested-tads.matrix.gz'
    verboseprint("writing tad matrix ...",end="")
    writeMatrix(header_rows,header_cols,matrix,tadMatrixFile)
    verboseprint("done")
    
def fillTadMatrix(header_rows,header_cols,matrix,tads,header2idx,nan_rowcols):
    
    for i,t in enumerate(tads):
        tad_start_header,tad_end_header=t[1]
        tad_strength=t[2]
        
        tad_start_idx=header2idx[tad_start_header]
        tad_end_idx=header2idx[tad_end_header]
        
        #print(i,tad_start_header,tad_start_idx,tad_end_header,tad_end_idx,tad_strength,sep="\t")
        
        for y in xrange(tad_start_idx,tad_end_idx):
            for x in xrange(tad_start_idx,tad_end_idx):
                matrix[y,x]+=tad_strength
        
        # fill all nans
        for i in nan_rowcols:
            idx=header2idx[i]
            matrix[idx,:]=np.nan
            matrix[:,idx]=np.nan
            
    return matrix
    
def deGroupChr(chr_id):
    return(chr_id.split('-')[0])
    
def assemble_tads(boundary_noise,boundary_ref,chr_boundaries,chr_range,chr_idx,tads,field2index,header2idx,insulation):

    chr_start,chr_end=chr_range[chr_idx]
    
    while(len(chr_boundaries) != 0):
        tmp_strength,tmp_index=chr_boundaries.pop()
       
        #print("")
        #print("STARTING TAD #",tmp_index,boundary_ref[tmp_index]["header"]," ... ")
        
        boundary_ref[tmp_index]["available"]=False
        tads=create_tad(boundary_noise,boundary_ref,tmp_index,chr_range[chr_idx],tads,header2idx,insulation)
    
    return boundary_ref,chr_boundaries,chr_range,chr_idx,tads,field2index
    
def create_tad(boundary_noise,boundary_ref,anchor_idx,chr_bound,tads,header2idx,insulation):
    
    left_bound=None
    left_bound_header=None
    left_bound_boundary_header=None
    left_ref=np.nonzero(
        (boundary_ref["available"]==True) & 
        ((boundary_ref["index"]>=chr_bound[0]) & (boundary_ref["index"]<=chr_bound[1])) & 
        (boundary_ref["index"]<anchor_idx) )[0]
    left_size=len(left_ref)
    
    if left_size > 0:
        left_bound=np.max(left_ref)
        left_bound_boundary_header=boundary_ref[left_bound]["header"]
        left_bound_header=boundary_ref[left_bound]["boundaryHeader"]
        
    right_bound=None
    right_bound_header=None
    right_bound_boundary_header=None
    right_ref=np.nonzero( 
        (boundary_ref["available"]==True) & 
        ((boundary_ref["index"]>=chr_bound[0]) & (boundary_ref["index"]<=chr_bound[1])) & 
        (boundary_ref["index"]>anchor_idx) )[0]
    right_size=len(right_ref)
    
    if right_size > 0:
        right_bound=np.min(right_ref)
        right_bound_boundary_header=boundary_ref[right_bound]["header"]
        right_bound_header=boundary_ref[right_bound]["boundaryHeader"]
    
    anchor_boundary_header=boundary_ref[anchor_idx]["header"]
    anchor_header=boundary_ref[anchor_idx]["boundaryHeader"]
    anchor_strength=boundary_ref[anchor_idx]["boundaryStrength"]
    
    verboseprint("\t","ANCHOR",anchor_idx,anchor_strength)
    verboseprint("\t",left_bound_boundary_header,"::",anchor_boundary_header,"::",right_bound_boundary_header)
    verboseprint("\t",left_bound_header,"::",anchor_header,"::",right_bound_header)
    verboseprint("\tleft",left_ref)
    verboseprint("\tright",right_ref)
    
    # search potential boundaries to the right
    
    left_idx=None
    if left_bound != None:
        for i in left_ref[::-1]:
            #print("\t\tleft searching",i,"...")
            tmp_strength=boundary_ref[i]["boundaryStrength"]
            #print("\t\t\t","left",anchor_idx,left_bound,i,anchor_strength,"vs",tmp_strength,"(",anchor_strength-tmp_strength,") [",boundary_noise,"]")
 
            if(tmp_strength > (anchor_strength-boundary_noise)) or (abs(anchor_strength-tmp_strength) > boundary_noise):
                #print("\t\t\t\t","found a potential TAD!",anchor_strength,tmp_strength)
                left_idx=i
                break
            # else keep looking
                
        if left_idx != None:
            #print("\t\t\t\t","good left idx",left_idx)
            tad_start_header_idx=header2idx[boundary_ref[left_idx]["boundaryHeader"]]
            tad_end_header_idx=header2idx[boundary_ref[anchor_idx]["boundaryHeader"]]
            tad_insulation=np.array(insulation[tad_start_header_idx:tad_end_header_idx+1])
            na_pc=np.float((np.count_nonzero(np.isnan(tad_insulation)))/tad_insulation.shape[0])
            tad_strength=np.nanmax(tad_insulation)-np.nanmin(tad_insulation)
            con_nan=num_consecutive_nan(tad_insulation)
            
            #print("\t\t\t\t",na_pc,con_nan)
            if na_pc < 0.25 and con_nan < 10:
                left_tad=[boundary_ref[left_idx]["header"],boundary_ref[anchor_idx]["header"]]
                left_tad_headers=[boundary_ref[left_idx]["boundaryHeader"],boundary_ref[anchor_idx]["boundaryHeader"]]
                tads.append((left_tad,left_tad_headers,tad_strength))
                verboseprint("\tleft_tad",left_tad,tad_strength)
    
    # search potential boundaries to the right
    
    right_idx=None    
    if right_bound != None:
        for i in right_ref:
           # print("\t\tright searching",i,"...")
            tmp_strength=boundary_ref[i]["boundaryStrength"]
            #print("\t\t\t","right",anchor_idx,right_bound,i,anchor_strength,"vs",tmp_strength,"(",anchor_strength-tmp_strength,") [",boundary_noise,"]")
            
            if(tmp_strength > (anchor_strength-boundary_noise)) or (abs(anchor_strength-tmp_strength) > boundary_noise):
                #print("\t\t\t\t","found a potential TAD!",anchor_strength,tmp_strength)
                right_idx=i
                break
                
        if right_idx != None:
            #print("\t\t\t\t","good right idx",right_idx)
            tad_start_header_idx=header2idx[boundary_ref[anchor_idx]["boundaryHeader"]]
            tad_end_header_idx=header2idx[boundary_ref[right_idx]["boundaryHeader"]]
            tad_insulation=np.array(insulation[tad_start_header_idx:tad_end_header_idx+1])
            na_pc=np.float((np.count_nonzero(np.isnan(tad_insulation)))/tad_insulation.shape[0])
            tad_strength=np.nanmax(tad_insulation)-np.nanmin(tad_insulation)
            con_nan=num_consecutive_nan(tad_insulation)
            
            #print("\t\t\t\t",na_pc,con_nan)
            if na_pc < 0.25 and con_nan < 10:
                right_tad=[boundary_ref[anchor_idx]["header"],boundary_ref[right_idx]["header"]]
                right_tad_headers=[boundary_ref[anchor_idx]["boundaryHeader"],boundary_ref[right_idx]["boundaryHeader"]]
                tads.append((right_tad,right_tad_headers,tad_strength))
                verboseprint("\tright_tad",right_tad,tad_strength)
                    
    return(tads)
    
def num_consecutive_nan(arr):
    max_con_nan = 0
    con_nan = [len(list(v)) for i, v in itertools.groupby(np.isnan(arr)) if i]
    if(len(con_nan) > 0):
        max_con_nan = max(con_nan)
    
    return max_con_nan
     
def splitHeader(header):
    m=re.search(r'(\S+)\|(\S+)\|(\S+):(\d+)-(\d+)',header)
    if m==None:
        sys.exit('error: incorrect header format ['+str(header)+']!')  

    header_name,header_assembly,header_chr,header_start,header_end=m.groups()
    
    return(header_name,header_assembly,header_chr,header_start,header_end)
        
def load_boundary_file(boundaryFile):

    init=1
    current_chr=0
    chr2index=dict()
    index2chr=dict()
    
    index2field=dict()
    field2index=dict()
    num_boundaries=0
    
    b_fh=input_wrapper(boundaryFile)
    for line in b_fh:
        l=line.rstrip("\n").split("\t")
        
        if line.startswith("#"):
            continue
        
        if init == 1:
            index2field=dict(enumerate(l))
            field2index=dict((value, key) for key, value in index2field.iteritems())
            init=0
            continue
      
        header=l[field2index["header"]]
        header_name,header_assembly,header_chr,header_start,header_end=splitHeader(header)
        
        if(header_chr not in chr2index):
            chr2index[header_chr]=current_chr
            index2chr[current_chr]=header_chr
            current_chr+=1
            
        chr_id=chr2index[header_chr]
        if(chr_id < (current_chr-1)):
            sys.exit('improperly sorted boundary file!')
        
        num_boundaries += 1
    
    b_fh.close()
    
    boundary_ref=np.empty(num_boundaries, 
        dtype={'names':['index', 'header', 'boundaryHeader', 'start', 'end', 'boundaryStrength','chr','available'],
        'formats':['int64','a500','a500','int64','int64','float64','int64','bool']})
   
    init=1
    i=0
    b_fh=input_wrapper(boundaryFile)
    for line in b_fh:
        l=line.rstrip("\n").split("\t")
        
        if line.startswith("#"):
            continue
        
        if init == 1:
            index2field=dict(enumerate(l))
            field2index=dict((value, key) for key, value in index2field.iteritems())
            init=0
            continue
      
        header=l[field2index["header"]]
        header_name,header_assembly,header_chr,header_start,header_end=splitHeader(header)
        
        if(header_chr not in chr2index):
            chr2index[header_chr]=current_chr
            index2chr[current_chr]=header_chr
            current_chr+=1
            
        chr_id=chr2index[header_chr]
        if(chr_id < (current_chr-1)):
            sys.exit('improperly sorted boundary file!')
        

        boundary_ref[i]=(i,l[field2index["header"]],l[field2index["boundaryHeader"]],l[field2index["start"]],l[field2index["end"]],l[field2index["boundaryInsulation"]],chr_id,1)
        
        i += 1
    
    b_fh.close()
    
    # if using insulation score
    boundary_ref["boundaryStrength"]=boundary_ref["boundaryStrength"]*-1
    boundary_ref["boundaryStrength"]=boundary_ref["boundaryStrength"]+abs(min(boundary_ref["boundaryStrength"]))
        
    return(num_boundaries,chr2index,index2chr,index2field,field2index,boundary_ref)

      
def load_insulation_file(insulationFile):

    i_fh=input_wrapper(insulationFile)
    
    headers=[]
    
    header2idx=dict()
    insulation=[]
    
    init=1
    index2field=dict()
    field2index=dict()
    nan_rowcols=[]
    
    i=0
    for line in i_fh:
        l=line.rstrip("\n").split("\t")
        
        if line.startswith("#"):
            continue
        
        if init == 1:
            index2field=dict(enumerate(l))
            field2index=dict((value, key) for key, value in index2field.iteritems())
            init=0
            continue
      
        header=l[field2index["header"]]
        insulation_score=(l[field2index["insulationScore"]])
        if insulation_score == 'NA':
            insulation_score=np.nan
            nan_rowcols.append(header)
        insulation_score=float(insulation_score)
        
        header2idx[header]=i
        insulation.append(insulation_score)
        headers.append(header)
        
        i += 1
        
        
    return(headers,headers,header2idx,insulation,nan_rowcols)
        
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

def writeMatrix(header_rows,header_cols,matrix,matrixFile,precision=4):
    """
    write a np matrix with row/col headers - my5C file format - txt formatted gzipped file
    """
    
    nrows=len(header_rows)
    ncols=len(header_cols)
    
    # interaction matrix output
    out_fh=gzip.open(matrixFile,"wb")
    
    # write matrix col headers
    header=[str(i) for i in header_cols]
    print(str(nrows)+"x"+str(ncols)+"\t"+"\t".join(header),file=out_fh)

    format_func=("{:0."+str(precision)+"f}").format
    
    k=0
    
    for i in xrange(nrows):
        print(header_rows[i]+"\t"+"\t".join(map(format_func,matrix[i,:])),file=out_fh)
    
    out_fh.close()

if __name__=="__main__":
      main()
