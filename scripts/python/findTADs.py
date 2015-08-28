#!/usr/local/bin/python
"""
***********************************************
- PROGRAM: findTADs.py
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

from scipy.stats.stats import nanmean
from scipy.spatial.distance import squareform
from scipy.cluster.hierarchy import linkage,dendrogram
from scipy.cluster.vq import kmeans, vq, whiten
import scipy.spatial.distance as dist

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
    inputMatrix = args.inputMatrix
    nclust = args.nclust
        
    if not os.path.isfile(inputMatrix):
        sys.exit('invalid input file! (non-existant)')
    
    print("inputMatrix",inputMatrix)
    inputMatrixName=os.path.basename(inputMatrix)
    inputMatrixName=re.sub(".gz", "", inputMatrixName)    
    inputMatrixName=re.sub(".matrix", "", inputMatrixName)
    print("inputMatrixName",inputMatrixName)
    
    print("")
    
    print("loading matrix ... ",end="")
    if inputMatrix.endswith('.gz'):
        infh=gzip.open(inputMatrix,'r')
    else:
        infh=open(inputMatrix,'r')
    
    matrix,header_rows,header_cols = load_matrix((l for l in infh if not l.startswith('#')), hrows=1, hcols=1) # since this returns data, header_rows and header_cols
    infh.close()
    print("done")
    
    nrows=len(header_rows)
    ncols=len(header_cols)
    
    print("")
    
    # enfore symmetrical matrices only
    enforceSymmetrical(matrix)
    
    # get assembly from all headers
    assembly=getMatrixAssembly(header_rows+header_cols)
    
    # find nan rows
    print("finding nan rows ... ",end="\n")
    nan_rows=np.sum(np.isnan(matrix),axis=0)==matrix.shape[0]
    nan_cols=np.sum(np.isnan(matrix),axis=1)==matrix.shape[1]
    nan_rowcols=nan_rows | nan_cols
    
    nans=(np.isnan(matrix))
    matrix[nans]=np.nan
    
    print("converting all 2D nan to 0 ... ",end="\n")
    matrix = np.nan_to_num(matrix)
    
    #z = linkage(matrix)
    #codebook, _ = kmeans(matrix,nclust,iter=1)  # three clusters
    #cluster_indices,cluster_distances = vq(matrix, codebook)
    
    
    for i in xrange(nrows):
        tmp_header=header_rows[i]
        header_name,header_assembly,header_chr,header_start,header_end=splitHeader(tmp_header)
        
        tmp_cluster=cluster_indices[i]
        tmp_cluster_dist=cluster_distances[i]
        print(header_chr+"\t"+header_start+"\t"+header_end+"\t"+str(tmp_cluster))
        
    #d = dendrogram(z)
    #matplotlib.pyplot.gcf()
    #matplotlib.pyplot.savefig("tmp.pdf")

    
    print("")


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

def enforceSymmetrical(matrix):
    
    nmatrix_rows=matrix.shape[0]
    nmatrix_cols=matrix.shape[1]
    
    if(nmatrix_rows != nmatrix_cols):
        sys.exit('non-symmetrical matrix!')
    
    all((matrix[y,:].T==matrix[y,:]).all() for y in range(matrix.shape[0]))
    all((matrix[:,x].T==matrix[:,x]).all() for x in range(matrix.shape[1]))

def splitHeader(header):
    m=re.search(r'(\S+)\|(\S+)\|(\S+):(\d+)-(\d+)',header)
    if m==None:
        sys.exit('error: incorrect input format!')  

    header_name,header_assembly,header_chr,header_start,header_end=m.groups()
    
    return(header_name,header_assembly,header_chr,header_start,header_end)
        
def getMatrixAssembly(headers):
   
    assembly=None
    for i,header in enumerate(headers):
        m=re.search(r'(\S+)\|(\S+)\|(\S+):(\d+)-(\d+)',header)
        if m==None:
            sys.exit('error: incorrect input format!')  

        header_name,header_assembly,header_chr,header_start,header_end=m.groups()
        if assembly==None:
            assembly=header_assembly
        else:
            if assembly!=header_assembly:
                sys.exit('assembly/header_assembly is not constant!')
            assembly=header_assembly

    print("matrix assembly:",assembly)
    return(assembly)
   
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
    
    
    nrows=len(header_rows)
    ncols=len(header_cols)
    ndata_rows=data.shape[0]
    ndata_cols=data.shape[1]
    
    if(nrows != ndata_rows):
        sys.exit('header/data mismatch!')    
    if(ncols != ndata_cols):
        sys.exit('header/data mismatch!')    
        
    if (return_all or (hrows and hcols)):
        return data,header_rows,header_cols
    if(hrows):
        return data,header_rows
    if(hcols):
        return data,header_cols

    
    return data

def check_options():
    ''' Checks the options to the program '''

    # Create parser object
    parser = argparse.ArgumentParser(add_help=True,
        formatter_class=argparse.RawDescriptionHelpFormatter)

    # Add arguments 
    parser.add_argument('-i' , metavar='--inputMatrix'  , help="*Input matrix file", dest="inputMatrix", type=str, default="")
    parser.add_argument('-n' , metavar='--nclust'  , help="number of clusters", dest="nclust", type=int, default=1)
    
    # Parse command line with parse_args and store it in an object
    args = parser.parse_args()

    return args

if __name__=="__main__":
      main()
