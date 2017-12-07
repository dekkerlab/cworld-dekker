#!/usr/local/bin/python

"""
PCA on supplied matrix.  Extract PC1, PC2, PC3.  Works best on distance normalized matrix.
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

# HAS BEEN COMMENTED LONG BEFORE 2017
# For eigenvectors and eigenvalues
#from scipy.stats.stats import nanmean
#from scipy import linalg as la
#from scipy import weave 

verboseprint=lambda *a, **k: None
__version__ = "1.0"
debug = None

def main():
    
    parser=argparse.ArgumentParser(description='Extract c-data from HDF5 file into TXT (matrix.gz)',formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('-i', '--input', dest='inputMatrix', type=str, required=True, help='interaction matrix hdf5 file')
    parser.add_argument('-r', '--refseq', dest='refSeqFile', type=str, required=True, help='refseq file to calculate gene density per bin/PC')
    parser.add_argument('-v', '--verbose', dest='verbose', action='count', help='Increase verbosity (specify multiple times for more)')
    parser.add_argument('--version', action='version', version='%(prog)s 1.0')
    
    args=parser.parse_args()

    inputMatrix=args.inputMatrix
    refSeqFile=args.refSeqFile
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
    
    inputMatrix_name=os.path.basename(inputMatrix)
    inputMatrix_name=re.sub(".gz", "", inputMatrix_name)    
    inputMatrix_name=re.sub(".matrix", "", inputMatrix_name)
    
    verboseprint("",file=sys.stderr)
    
    verboseprint("loading matrix ... ",end="",file=sys.stderr)
    infh=input_wrapper(inputMatrix)
    
    matrix,header_rows,header_cols = load_matrix((l for l in infh if not l.startswith('#')), hrows=1, hcols=1) # since this returns data, header_rows and header_cols
    infh.close()
    verboseprint("done",file=sys.stderr)
    
    verboseprint("",file=sys.stderr)
    if(len(header_rows) != len(header_cols)):
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
            
    try:
        with input_wrapper(refSeqFile) as rsfh:
            pass
    except IOError as e:
        sys.exit("invalid refSeq file! ("+refSeqFile+")")
    
    # get number of rows/col (assuming symmetrical)
    nrows=matrix.shape[0]
    ncols=matrix.shape[1]
    
    # find nan rows
    verboseprint("finding nan rows ... ",end="",file=sys.stderr)
    nan_rowcols = np.sum(np.isnan(matrix),0)==matrix.shape[0]
    valid_rowcols=np.invert(nan_rowcols)
    
    # remove nan rows
    # numpy negation with "~" more common, "-" deprecated
    matrix=matrix[~nan_rowcols,:][:,~nan_rowcols]
    verboseprint("done",file=sys.stderr)
    
    # convert all nan to 0
    verboseprint("converting all 2D nan to 0 ... ",end="",file=sys.stderr)
    matrix = np.nan_to_num(matrix)
    verboseprint("done",file=sys.stderr)
    
    # calculate corrcoef matrix
    verboseprint("calculating coorcoef ... ",end="",file=sys.stderr)
    corrMatrix = np.corrcoef(matrix)
    verboseprint("done",file=sys.stderr)
    
    verboseprint("")
    
    # do eigenvector analysis
    verboseprint("running PCA ... ",end="",file=sys.stderr)
    pca_score,pca_v = calculate_eigen(corrMatrix, 3)
    verboseprint("done",file=sys.stderr)
    
    verboseprint("\teigen1","{:.12%}".format(pca_score[0]))
    verboseprint("\teigen2","{:.12%}".format(pca_score[1]))
    verboseprint("\teigen3","{:.12%}".format(pca_score[2]))
    
    verboseprint("")
    
    eig1=pca_v[0]
    eig2=pca_v[1]
    eig3=pca_v[2]
    
    # pre-populate with nan
    egv1=np.nan*np.ones(nrows)
    egv2=np.nan*np.ones(nrows)
    egv3=np.nan*np.ones(nrows)
    
    egv1[~nan_rowcols]=eig1
    egv2[~nan_rowcols]=eig2
    egv3[~nan_rowcols]=eig3

    geneDensity=np.nan*np.ones(nrows)
    
    compartmentFile=inputMatrix_name+".compartments"    
    writeCompartmentFile(egv1,egv2,egv3,pca_score[0:3],geneDensity,header_rows,compartmentFile)
    
    verboseprint("")
    
    verboseprint("intersecing compartments with ref seq ... ",end="",file=sys.stderr)
    compartmentRefSeqFile=compartmentFile+".refSeq.txt"
    os.system("bedtools intersect -a "+compartmentFile+" -b "+refSeqFile+" -c > "+compartmentRefSeqFile)
    verboseprint("done",file=sys.stderr)
    
    eigenMultiplier,geneDensity = detectActiveCompartment(compartmentRefSeqFile)
    #os.system("rm "+compartmentRefSeqFile)

    verboseprint("\tflipping vectors by",eigenMultiplier," ... ",end="",file=sys.stderr)
    egv1 *= eigenMultiplier
    egv2 *= eigenMultiplier
    egv3 *= eigenMultiplier
    verboseprint("done",file=sys.stderr)
    
    verboseprint("")
    
    writeCompartmentFile(egv1,egv2,egv3,pca_score[0:3],geneDensity,header_rows,compartmentFile)
    
    eig1BedGraphFile=inputMatrix_name+".eigen1.bedGraph"    
    writeBedGraphFile(egv1,pca_score[0],header_rows,inputMatrix_name,eig1BedGraphFile)
    
    verboseprint("drawing eigen plot (",inputMatrix_name,") ... ",end="",file=sys.stderr)
    eigenPlot = scriptPath+"/R/plotEigen.R"
    os.system("Rscript "+eigenPlot+" `pwd` "+compartmentFile+" "+inputMatrix_name+" > /dev/null")
    verboseprint("done",file=sys.stderr)
    
    verboseprint("drawing evr plot (",inputMatrix_name,") ... ",end="",file=sys.stderr)
    evrFile=inputMatrix_name+".evr.txt"    
    writePCAevr(pca_score,evrFile)
    evrPlot = scriptPath+"/R/plotEVR.R"
    os.system("Rscript "+evrPlot+" `pwd` "+evrFile+" "+inputMatrix_name+" > /dev/null")
    verboseprint("done",file=sys.stderr)
    
    collapsed_corrMatrixFile=inputMatrix_name+'.collapsed.correlation.matrix.gz'
    verboseprint("writing collapsed_corrcoef matrix ...",end="",file=sys.stderr)
    writeMatrix(header_rows[np.where(valid_rowcols)],header_cols[np.where(valid_rowcols)],corrMatrix,collapsed_corrMatrixFile)
    verboseprint("done",file=sys.stderr)
    
    valid_rowcols=np.c_[valid_rowcols].T
    expanded_corrMatrix=np.zeros([nrows,ncols])
    expanded_corrMatrix.fill(np.nan)
    expanded_corrMatrix[np.where(valid_rowcols&valid_rowcols.T)]=corrMatrix.flatten()
    
    corrMatrixFile=inputMatrix_name+'.correlation.matrix.gz'
    verboseprint("writing corrcoef matrix ...",end="",file=sys.stderr)
    writeMatrix(header_rows,header_cols,expanded_corrMatrix,corrMatrixFile)
    verboseprint("done",file=sys.stderr)
    
    verboseprint("",file=sys.stderr)


def writePCAevr(pca_score,outfile):
    out_fh=output_wrapper(outfile)
    print("eigenvector","\t","evr",sep="",file=out_fh)
    
    for i,evr in enumerate(pca_score):
        print(i+1,evr,sep="\t",file=out_fh)
    
    out_fh.close()
    
        
def writeBedGraphFile(egv1,evr,header_rows,name,outfile):
    "write the compartment file"
    
    verboseprint("writing bed graph file (",outfile,") ... ",end="",file=sys.stderr)
    
    # "nanmax" expects iterable as 1st argument, second is axis
    # wrapping abs(min(egv1)) and max(egv1) as a 2-element list:
    yBound=np.nanmax([abs(np.nanmin(egv1)),np.nanmax(egv1)])
    yBound *= 1.25
    yBound=round(yBound,5)
    
    out_fh=output_wrapper(outfile,suppress_comments=True)
    print("track type=bedGraph name='"+name+"-evr:"+str(evr)+"%' description='"+name+"-evr:"+str(evr)+"%' maxHeightPixels=128:64:32 visibility=full autoScale=off viewLimits="+str(-yBound)+":"+str(yBound)+" color=0,255,0 altColor=255,0,0",end="\n",file=out_fh)

    for i,header in enumerate(header_rows):
        m=re.search(r'(\S+)\|(\S+)\|(\S+):(\d+)-(\d+)',header)
        if m==None:
            sys.exit('error: incorrect input format!')  

        bin_id,genome,chr_id,bin_start,bin_end=m.groups()
        
        chr_id=chr_id.split("-")[0]
        
        eigen1=round(egv1[i],5)
        
        print(str(chr_id)+"\t"+str(bin_start)+"\t"+str(bin_end)+"\t"+str(eigen1),end="\n",file=out_fh)
  
    out_fh.close()
    verboseprint("done",file=sys.stderr)
    
def writeCompartmentFile(egv1,egv2,egv3,pca_score,geneDensity,header_rows,outfile):
    "write the compartment file"
        
    nan_geneDensity=np.sum(np.isnan(geneDensity))
    
    out_fh=output_wrapper(outfile,suppress_comments=True)
    verboseprint("writing eigenvector file (",outfile,") ... ",end="",file=sys.stderr)
    
    if len(geneDensity)==nan_geneDensity:
        print("#chr\tstart\tend\tname\tindex\teigen1\teigen1evr\teigen2\teigen2evr\teigen3\teigen3evr",end="\n",file=out_fh)
    else:
        print("#chr\tstart\tend\tname\tindex\teigen1\teigen1evr\teigen2\teigen2evr\teigen3\teigen3evr\tgeneDensity",end="\n",file=out_fh)
        
    for i,header in enumerate(header_rows):
        m=re.search(r'(\S+)\|(\S+)\|(\S+):(\d+)-(\d+)',header)
        if m==None:
            sys.exit('error: incorrect input format!')  

        bin_id,genome,chr_id,bin_start,bin_end=m.groups()
        
        eigen1=egv1[i]
        eigen2=egv2[i]
        eigen3=egv3[i]
        
        if len(geneDensity)==nan_geneDensity:
            print(str(chr_id)+"\t"+str(bin_start)+"\t"+str(bin_end)+"\t"+str(header)+"\t"+str(i)+"\t"+str(eigen1)+"\t"+str(pca_score[0])+"\t"+str(eigen2)+"\t"+str(pca_score[1])+"\t"+str(eigen3)+"\t"+str(pca_score[2]),end="\n",file=out_fh)
        else:
            nGenes=geneDensity[i]
            print(str(chr_id)+"\t"+str(bin_start)+"\t"+str(bin_end)+"\t"+str(header)+"\t"+str(i)+"\t"+str(eigen1)+"\t"+str(pca_score[0])+"\t"+str(eigen2)+"\t"+str(pca_score[1])+"\t"+str(eigen3)+"\t"+str(pca_score[2])+"\t"+str(nGenes),end="\n",file=out_fh)
    
    out_fh.close()
    
    verboseprint("done",file=sys.stderr)
    
def detectActiveCompartment(file):
    "detect the active compartment - overlap +/- with gene density"
        
    eigenMultiplier=1
    
    infh=input_wrapper(file)
    #chr    start   end     name    eigen1  eigen2  eigen3 geneCount    
    
    geneDensity=[]
    
    posSum=0
    posCount=0
    negSum=0
    negCount=0
    
    for i,x in enumerate(infh):
        a=x.rstrip("\n").split("\t")
        
        nGenes=float(a[11])
        geneDensity.append(nGenes)
        
        if a[5] == 'nan':
            continue
        
        eigen1=float(a[5])
        
        # skip eigen == 0
        if eigen1 > 0:
            posSum += (nGenes*abs(eigen1))
            posCount += 1
            
        if eigen1 < 0:
            negSum += (nGenes*abs(eigen1))
            negCount += 1
            
    posAvg=1
    if posCount > 0:
        posAvg=(posSum/posCount)
    negAvg=1
    if negCount > 0:
        negAvg=(negSum/negCount)
    
    if negSum > posSum:
        eigenMultiplier=-1
        
    verboseprint("\tposSum",posSum,"posCount",posCount,"posAvg",posAvg,file=sys.stderr)
    verboseprint("\tnegSum",negSum,"negCount",negCount,"negAvg",negAvg,file=sys.stderr)
    verboseprint("\teigenMultiplier",eigenMultiplier,file=sys.stderr)
    
    return eigenMultiplier,geneDensity
    
def calculate_eigen(A, numPCs = 3):
    """performs eigen vector analysis, and returns 3 best principal components
    result[0] is the first PC, etc"""    
    
    #A = np.array(A,float)
    #M = (A-np.mean(A.T,axis=1)).T 
    #covM = np.dot(M,M.T)
    #[latent,coeff] =  scipy.sparse.linalg.eigsh(covM,numPCs)
    #return (np.transpose(coeff[:,::-1]),latent[::-1])
    
    #egv_data = calculate_eigen(corrMatrix, 3)
    
    #eig1=egv_data[0][0]
    #eig2=egv_data[0][1]
    #eig3=egv_data[0][2]
 
    # pre-populate with nan
    #egv1=np.nan*np.ones(nrows)
    #egv2=np.nan*np.ones(nrows)
    #egv3=np.nan*np.ones(nrows)
    
    #egv1[~nan_rowcols]=eig1
    #egv2[~nan_rowcols]=eig2
    #egv3[~nan_rowcols]=eig3
    
    ncomp=min(100,A.shape[0])
    pca = decomposition.PCA(n_components=ncomp)
    pca.fit(A)
    PCA(copy=True, n_components=3, whiten=False)
    pca_score=pca.explained_variance_ratio_
    pca_v = pca.components_[0:3]
        
    return(pca_score,pca_v)
    

def load_matrix(fh,hrows=0,hcols=0,np_dtype='float32',row_block_size=1000,numpy_mode=True,max_rows=None,verbose=False,return_all=False,pad=None):
    """
    From Noam Kaplan (noamlib)
    load a np.array or a list of lists from a text file handle (but works with any iterator) or filename, more memory efficient than numpy.genfromtxt(), headers are returned as lists of strings
    """
    fh_from_filename=False
    
    if type(fh)==str:
        if (fh=='-'):
            fh=sys.stdin
        else:
            fh=input_wrapper(fh)
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
        return data,np.asarray(header_rows),np.asarray(header_cols)
    if(hrows):
        return data,np.asarray(header_rows)
    if(hcols):
        return data,np.asarray(header_cols)
    
    return data

  
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
