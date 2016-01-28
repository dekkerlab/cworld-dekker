#!/usr/local/bin/python

"""
vector to matrix
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
from datetime import datetime

import numpy as np
import scipy as sp
import scipy.stats
import itertools

from collections import *
from math import cos,log,sin,sqrt 
from sklearn.decomposition import PCA
from sklearn import decomposition

verboseprint=lambda *a, **k: None
__version__ = "1.0"
debug = None

def main():
    
    parser=argparse.ArgumentParser(description='convert 1 or 2 vectors into a matrix (TXT - matrix.gz)',formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('--vxy', '--vector_xy', dest='vector_xy', type=str, default=None, help='vector file - bedGraph or 1 col')
    parser.add_argument('--vx', '--vector_x', dest='vector_x', type=str, default=None, help='x axis vector file - bedGraph or 1 col')
    parser.add_argument('--vy', '--vector_y', dest='vector_y', type=str, default=None, help='y axis vector file - bedGraph or 1 col')
    parser.add_argument('-a', '--assembly', dest='assembly', type=str, default="NA", help='genome assembly')
    parser.add_argument('--cm', '--contrustion_mode', dest='construction_mode', type=str, default="mean", choices=['mean','multiply','add','subtract'], help='matrix constructor mode')
    
    parser.add_argument('-v', '--verbose', dest='verbose', action='count', help='Increase verbosity (specify multiple times for more)')
    parser.add_argument('--version', action='version', version='%(prog)s 1.0')
    
    args=parser.parse_args()

    vector_xy=args.vector_xy
    vector_x=args.vector_x
    vector_y=args.vector_y
    assembly=args.assembly
    construction_mode=args.construction_mode
    verbose=args.verbose
    
    log_level = logging.WARNING
    if verbose == 1:
        log_level = logging.INFO
    elif verbose >= 2:
        log_level = logging.DEBUG
    logging.basicConfig(level=log_level)
    
    global verboseprint
    verboseprint = print if verbose else lambda *a, **k: None
    
    if(vector_xy == None) and (vector_x != None) and (vector_y != None):
        if not os.path.isfile(vector_x):
            sys.exit('invalid vector_x file! (non-existant)')
        if not os.path.isfile(vector_y):
            sys.exit('invalid vector_y file! (non-existant)')
    elif(vector_x == None) and (vector_y == None) and (vector_xy != None):
        if not os.path.isfile(vector_xy):
            sys.exit('invalid vector_xy file! (non-existant)')
        vector_x=vector_xy
        vector_y=vector_xy
    else:
        sys.exit('incorrect usage')
        
    scriptPath=os.path.realpath(__file__)
    scriptPath="/".join(scriptPath.split("/")[0:-2])
    
    name_y=os.path.basename(vector_y)
    name_y=re.sub(".gz", "", name_y)    
    name_y=re.sub(".matrix", "", name_y)
    
    name_x=os.path.basename(vector_x)
    name_x=re.sub(".gz", "", name_x)    
    name_x=re.sub(".matrix", "", name_x)
    
    if(name_y == name_x):
        name=name_y+'__'+construction_mode
    else:
        name=name_y+'__'+name_x+'__'+construction_mode
        
    verboseprint("")
    
    verboseprint(name)

    verboseprint("")
    
    header_rows,vy=load_vector(vector_y,assembly)
    header_cols,vx=load_vector(vector_x,assembly)
    
    rows=vy.shape[0]
    cols=vx.shape[0]
    
    matrix=np.zeros((rows,cols),dtype="float32")
    
    verboseprint("building matrix ... ",end="")
    for vy_idx,vy_v in enumerate(vy):    
        for vx_idx,vx_v in enumerate(vx):
            if(construction_mode == 'mean'):
                matrix[vy_idx,vx_idx]=(vy_v+vx_v)/2
            elif(construction_mode == 'multiply'):
                matrix[vy_idx,vx_idx]=(vy_v*vx_v)
                # speed up later via matrix math
            elif(construction_mode == 'add'):
                matrix[vy_idx,vx_idx]=(vy_v+vx_v)
            elif(construction_mode == 'subtract'):
                matrix[vy_idx,vx_idx]=(vy_v-vx_v)
    verboseprint("done")
    
    verboseprint("")
    
    matrixFile=name+'.matrix.gz'
    verboseprint("writing matrix ... ",end="")
    writeMatrix(header_rows,header_cols,matrix,matrixFile)
    verboseprint("done")
    
def load_vector(v,assembly):
    
    with input_wrapper(v)as fh:
        
        headers=[]
        vector=[]
        i=0
        for line in fh:
            l=line.rstrip("\n").split("\t")
        
            if line.startswith("#") or line.startswith("track"):
                continue
            
            if(len(l) == 4):
                score=l[3]
                header='vector_'+str(i)+'|'+str(assembly)+'|'+str(l[0])+':'+str(l[1])+'-'+str(l[2])
            elif(len(l) == 1):
                score=l[1]
                header='vector_'+str(i)
                
            if score == 'NA':
                score=np.nan
            else:
                score=np.float(score)
                
            headers.append(header)
            vector.append(score)
            i+=1
    
    return(headers,np.array(vector))
    
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
    out_fh=output_wrapper(matrixFile)
    
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