#!/usr/local/bin/python
"""
***********************************************
- PROGRAM: matrix2tab.py
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
import shutil
import gzip
import re
import os
import math
import uuid
import socket
import itertools
import collections
from datetime import datetime

import re
import os

from scipy.stats.stats import nanmean
import numpy as np
import scipy as sp

# For eigenvectors and eigenvalues
from scipy import linalg as la
from math import cos,log,sin,sqrt 
from scipy import weave 

verboseprint=lambda *a, **k: None
__version__ = "1.0"

def main():
    print("")
    
    parser=argparse.ArgumentParser   (description='convert a cworld/my5C fomatted tsv matrix file to a 3 col tab matrix [symmetrical only]',formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('-i',dest='inputMatrix',type=str,required=True,help='interaction matrix hdf5 file')
    parser.add_argument('-v', '--verbose', dest='verbose',  action='count', help='Increase verbosity (specify multiple times for more)')
    
    args=parser.parse_args()
    
    inputMatrix=args.inputMatrix
    verbose=args.verbose

    log_level = logging.WARNING
    if verbose == 1:
        log_level = logging.INFO
    elif verbose >= 2:
        log_level = logging.DEBUG
    logging.basicConfig(level=log_level)
    
    global verboseprint
    verboseprint = print if verbose else lambda *a, **k: None
    
    if not os.path.isfile(inputMatrix):
        sys.exit('invalid input file! (non-existant)')
    
    scriptPath=os.path.realpath(__file__)
    scriptPath="/".join(scriptPath.split("/")[0:-2])
   
    verboseprint("inputMatrix",inputMatrix)
    inputMatrix_name=os.path.basename(inputMatrix)
    inputMatrix_name=re.sub(".gz", "", inputMatrix_name)    
    inputMatrix_name=re.sub(".matrix", "", inputMatrix_name)
    verboseprint("inputMatrix_name",inputMatrix_name)
    
    verboseprint("")
    
    verboseprint("loading matrix ... ",end="")
    if inputMatrix.endswith('.gz'):
        infh=gzip.open(inputMatrix,'r')
    else:
        infh=open(inputMatrix,'r')
    
    matrix,header_rows,header_cols = load_matrix((l for l in infh if not l.startswith('#')), hrows=1, hcols=1) # since this returns data, header_rows and header_cols
    infh.close()
    verboseprint("done")
    
    verboseprint("")
    
    # get number of rows/col (assuming symmetrical)
    nrows=len(header_rows)
    ncols=len(header_cols)
    nmatrix_rows=matrix.shape[0]
    nmatrix_cols=matrix.shape[1]
    verboseprint("rows:",nrows,nmatrix_rows)
    verboseprint("cols:",ncols,nmatrix_cols)
    
    if(nrows != ncols):
        sys.exit('non-symmetrical matrix!')
    if(nmatrix_rows != nmatrix_cols):
        sys.exit('non-symmetrical matrix!')
    if(nrows != nmatrix_rows):
        sys.exit('non-symmetrical matrix!')    
    if(ncols != nmatrix_cols):
        sys.exit('non-symmetrical matrix!')
    
    assembly=getHeaderAssembly(header_rows)
    equalSpacingFlag,equalSizingFlag,header_spacing,header_sizing=getHeaderSpacing(header_rows)
    num_headers=len(header_rows)
    
    verboseprint("")
    
    matrixFile=inputMatrix_name+".tab.gz"
    writeTab(header_rows,matrix,matrixFile)

    

def getHeaderAssembly(header_rows):
    
    assembly=None
    for i,header in enumerate(header_rows):
        m=re.search(r'(\S+)\|(\S+)\|(\S+):(\d+)-(\d+)',header)
        if m==None:
            sys.exit('error: incorrect input format!')  

        bin_id,genome,chr_id,bin_start,bin_end=m.groups()
        if assembly==None:
            assembly=genome
        else:
            if assembly!=genome:
                sys.exit('assembly/genome is not constant!')
            assembly=genome
    
    return(assembly)
    
def getHeaderSpacing(header_rows):
    
    globalHeaderSpacingArr=[]
    globalHeaderSizingArr=[]
    equalSpacingFlag=1
    equalSizingFlag=1
    
    bin_size=None
    bin_step=None
    globalHeaderSpacing=globalHeaderSizing=-1
    for i in xrange(len(header_rows)-1):
        
        header=header_rows[i]
        headerObject=getHeaderObject(header)
        headerRegion=headerObject["region"]
        headerStart=headerObject["start"]
        headerEnd=headerObject["end"]
        headerSize=headerObject["size"]
        
        nextHeader=header_rows[i+1]
        nextHeaderObject=getHeaderObject(nextHeader)
        nextHeaderRegion=nextHeaderObject["region"]
        nextHeaderStart=nextHeaderObject["start"]
        nextHeaderEnd=nextHeaderObject["end"]
        nextHeaderSize=nextHeaderObject["size"]
        
        if((nextHeaderRegion != headerRegion) or (headerEnd == nextHeaderEnd) or (headerStart == nextHeaderStart)):
            continue
        
        if((globalHeaderSpacing != ((nextHeaderStart-headerStart))) and (globalHeaderSpacing != -1)):
            equalSpacingFlag=0 
        if((globalHeaderSizing != (headerSize)) and (globalHeaderSizing != -1)):
            equalSizingFlag=0 
        
        globalHeaderSpacing=(nextHeaderStart-headerStart)
        globalHeaderSizing=headerSize
                
        globalHeaderSpacingArr.append(globalHeaderSpacing)
        globalHeaderSizingArr.append(globalHeaderSizing)
    
    meanGlobalHeaderSpacing=int(np.mean(globalHeaderSpacingArr))
    meanGlobalHeaderSizing=int(np.mean(globalHeaderSizingArr))
    
    return(equalSpacingFlag,equalSizingFlag,meanGlobalHeaderSpacing,meanGlobalHeaderSizing)
        
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

def deGroupChr(chr_id):
    return(chr_id.split('-')[0])
    
def getHeaderObject(header,enforceValidHeaders=0):
    
    subName=assembly=coords=None
    tmp=header.split('|')
    if enforceValidHeaders and len(tmp) != 3:
        sys.exit('error: incorrect input format!')  
    subName,assembly,coords=tmp
    
    chromosome=pos=None
    tmp=coords.split(':')
    if enforceValidHeaders and len(tmp) != 2:
        sys.exit('error: incorrect input format!')  
    chromosome,pos=tmp
    
    region=chromosome
    primerType=fragmentNumber=None
    if("__" in subName):
        tmp=subName.split("__")
        region=tmp
    else:
        tmp=subName.split("_")
        if(len(tmp) == 5 and tmp[0] == "5C"):
            region=str(tmp[1])+"_"+str(tmp[2])
            primerType=int(tmp[3])
            fragmentNumber=int(tmp[4])
    
    start=end=None
    tmp=pos.split("-")
    if enforceValidHeaders and len(tmp) != 2:
        sys.exit('error: incorrect input format!')  
    start=int(tmp[0])
    end=int(tmp[1])
    
    size=((end-start)+1) # add to for 1-based positioning
    midpoint=((end+start)/2)
        
    headerObject=dict()
    
    headerObject["subName"]=subName
    headerObject["fragmentNumber"]=fragmentNumber
    headerObject["primerType"]=primerType
    headerObject["assembly"]=assembly
    headerObject["chromosome"]=chromosome
    headerObject["coords"]=coords
    headerObject["region"]=region
    headerObject["start"]=start
    headerObject["end"]=end
    headerObject["midpoint"]=midpoint
    headerObject["size"]=size
        
    return(headerObject)

def headers2tabs(headers):
    
    header_tabs=[]
    
    for i in headers:
        header_object=getHeaderObject(i)
        header_tabs.append([header_object["chromosome"],header_object["start"],header_object["end"]])

    return np.array(header_tabs)
    
def writeTab(header_rows,matrix,matrixFile,precision=4):
    """
    write a np matrix into tab 3-col matrix format 
    """
    
    nrows=len(header_rows)
    
    header_tabs=headers2tabs(header_rows)
    
    # interaction matrix output
    out_fh=gzip.open(matrixFile,"wb")
    
    format_func=("{:0."+str(precision)+"f}").format
    
    k=0
    
    for i in xrange(nrows):
        print("\t".join(header_tabs[i])+"\t"+"\t".join(map(format_func,matrix[i,:])),file=out_fh)
    
    out_fh.close()
    
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

    
def load_matrix(fh,hrows=0,hcols=0,np_dtype='float32',row_block_size=1000,numpy_mode=True,max_rows=None,verbose=False,return_all=False,pad=None):
    """
    load a np.array or a list of lists from a text file handle (but works with any iterator) or filename, more memory efficient than numpy.genfromtxt(), headers are returned as lists of strings
    """
    fh_from_filename=False
    
    
    if type(fh)==str:
        if (fh=='-'):
            fh=sys.stdin
        else:
            fh=open(fh,'r')
            fh_from_filename=True

    original_fh=fh
        
    # init

    firstline=fh.next()
    
    fh=itertools.chain([firstline],fh)
    
    cols=len(firstline.rstrip("\n").split("\t"))
    rows=row_block_size

    if (max_rows!=None and max_rows<row_block_size):
        rows=max_rows

    if(hcols):
        cols-=hcols

   
    if numpy_mode:
        data=np.zeros((rows,cols),dtype=np_dtype)
    else:
        data=[]

    header_rows=[[] for i in range(hrows)]

    for i in range(hrows):
        header_rows[i]=fh.next().rstrip("\n").split("\t")[hcols:]
 
    header_cols=[[] for i in range(hcols)]
    
    # fill one line at a time

    prev_cols=-1

    r=0

    if (max_rows==None or r<max_rows):
  
        for i in fh:
            line=i.rstrip("\n").split("\t")

            cols=len(line)-hcols

           # if(cols==0):
           #     sys.exit('no valid columns in input line '+str(r))

            if(prev_cols>-1 and cols!=prev_cols):
                if(pad and cols<prev_cols):
                    line=line+['']*(prev_cols-cols)
                    cols=len(line)-hcols
                else:
                    sys.exit('inconsistent number of columns in input line '+str(r))

            prev_cols=cols

            if numpy_mode:
                not_allowed = ['','NA']
                try: # if np_dtype does not except ''or 'NA' as a value
                    np.dtype(np_dtype).type(not_allowed)
                except ValueError:
                    try:
                        np.dtype(np_dtype).type('nan')
                        line=[('nan' if i in not_allowed else i) for i in line] # '' or 'NA' are replaced with 'nan'
                    except ValueError:
                        pass
        
                
            for j in range(hcols):
                header_cols[j].append(line[j])

            if numpy_mode:
                data[r,:]=line[hcols:]

                # enlarge data if needed
                if(r==(data.shape[0]-1)):
                    data=np.resize(data,(data.shape[0]+row_block_size,cols))
                    rows=data.shape[0]

            else:
                data.append(line[hcols:]) 

            r+=1

            if (max_rows!=None and r>=max_rows):
                break

    rows=r

    if numpy_mode:
        data=np.resize(data,(rows,cols))

    if (fh_from_filename):
        original_fh.close()

    if (hcols==1):
        header_cols=header_cols[0]
        
    if (hrows==1):
        header_rows=header_rows[0]

    if(verbose):
        sys.stderr.write("loaded matrix with dimensions ("+str(len(data))+","+str(cols)+")\n")
    
    if (return_all or (hrows and hcols)):
        return data,header_rows,header_cols
    if(hrows):
        return data,header_rows
    if(hcols):
        return data,header_cols

    
    return data

if __name__=="__main__":
      main()
