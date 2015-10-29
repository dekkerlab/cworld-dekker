#!/usr/local/bin/python
"""
***********************************************
- PROGRAM: getEigenVectors.py
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
from collections import defaultdict
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
    
    parser=argparse.ArgumentParser   (description='covert a matrix into a correlation matrix',formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('-i',dest='inputMatrix',type=str,required=True,help='interaction matrix hdf5 file',)
    parser.add_argument('-v', '--verbose', dest='verbose',  action='count', help='Increase verbosity (specify multiple times for more)')
    parser.add_argument('--fillnan',dest='fill_nan',action='store_true',help='fill all NAN with corr value')
    
    args=parser.parse_args()
    
    inputMatrix=args.inputMatrix
    verbose=args.verbose
    fill_nan=args.fill_nan

    log_level = logging.WARNING
    if verbose == 1:
        log_level = logging.INFO
    elif verbose >= 2:
        log_level = logging.DEBUG
    logging.basicConfig(level=log_level)
    
    global verboseprint
    verboseprint = print if verbose else lambda *a, **k: None
    
    verboseprint("\n",end="")

    if not os.path.isfile(inputMatrix):
        sys.exit('invalid input file! (non-existant)')
    
    verboseprint("inputMatrix",inputMatrix)
    inputMatrixName=os.path.basename(inputMatrix)
    inputMatrixName=re.sub(".gz", "", inputMatrixName)    
    inputMatrixName=re.sub(".matrix", "", inputMatrixName)
    verboseprint("inputMatrixName",inputMatrixName)
    
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

    verboseprint("matrix assembly:",assembly)
   
    verboseprint("")
    
    # find nan rows
    verboseprint("finding nan rows ... ",end="\n")
    nan_rows=np.sum(np.isnan(matrix),axis=0)==matrix.shape[0]
    nan_cols=np.sum(np.isnan(matrix),axis=1)==matrix.shape[1]
    nan_rowcols=nan_rows | nan_cols
    
    # convert all nan to 0
    verboseprint("converting all nans to 0 ... ",end="")
    matrix = np.nan_to_num(matrix)
    verboseprint("done")
    
    # calculate corrcoef matrix
    verboseprint("calculating coorcoef ... ",end="")
    corrMatrix = np.corrcoef(matrix)
    verboseprint("done")
   
    if not fill_nan:
        verboseprint("replacing nan ... ",end="")
        corrMatrix[nan_rows,:]=np.nan
        corrMatrix[:,nan_cols]=np.nan
        verboseprint("done")
    
    corrMatrixFile=inputMatrixName+'.correlation.matrix.gz'
    verboseprint("writing corrcoef matrix ...",end="")
    writeMatrix(header_rows,header_cols,corrMatrix,corrMatrixFile)
    verboseprint("done")
    
    verboseprint("")


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
