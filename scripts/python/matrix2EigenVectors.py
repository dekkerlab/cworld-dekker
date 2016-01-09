#!/usr/local/bin/python

"""
PCA on supplied matrix.  Extract PC1, PC2, PC3.  Works best on distance normalized matrix.
"""

from __future__ import print_function
from __future__ import division

# Built in modules
import argparse
import os.path
import sys
import gzip
import re
import logging

import numpy as np
import scipy as sp
import scipy.stats
import itertools

from collections import *
from math import cos,log,sin,sqrt 
from sklearn.decomposition import PCA
from sklearn import decomposition

# For eigenvectors and eigenvalues
#from scipy.stats.stats import nanmean
#from scipy import linalg as la
#from scipy import weave 

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
    if inputMatrix.endswith('.gz'):
        infh=gzip.open(inputMatrix,'r')
    else:
        infh=open(inputMatrix,'r')
    
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
        with open(refSeqFile) as rsfh:
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
    matrix=matrix[-nan_rowcols,:][:,-nan_rowcols]
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
    writeCompartmentFile(egv1,egv2,egv3,pca_score,geneDensity,header_rows,compartmentFile,verbose)
    
    verboseprint("")
    
    verboseprint("intersecing compartments with ref seq ... ",end="",file=sys.stderr)
    compartmentRefSeqFile=compartmentFile+".refSeq.txt"
    os.system("bedtools intersect -a "+compartmentFile+" -b "+refSeqFile+" -c > "+compartmentRefSeqFile)
    verboseprint("done",file=sys.stderr)
    
    eigenMultiplier,geneDensity = detectActiveCompartment(compartmentRefSeqFile,verbose)
    os.system("rm "+compartmentRefSeqFile)

    verboseprint("\tflipping vectors by",eigenMultiplier," ... ",end="",file=sys.stderr)
    egv1 *= eigenMultiplier
    egv2 *= eigenMultiplier
    egv3 *= eigenMultiplier
    verboseprint("done",file=sys.stderr)
    
    verboseprint("")
    
    writeCompartmentFile(egv1,egv2,egv3,pca_score,geneDensity,header_rows,compartmentFile,verbose)
    
    eig1BedGraphFile=inputMatrix_name+".eigen1.bedGraph"    
    writeBedGraphFile(egv1,pca_score[0],header_rows,inputMatrix_name,eig1BedGraphFile,verbose)
    
    verboseprint("drawing plot (",inputMatrix_name,") ... ",end="",file=sys.stderr)
    drawPlot = scriptPath+"/R/plotEigen.R"
    os.system("Rscript "+drawPlot+" `pwd` "+compartmentFile+" "+inputMatrix_name+" 0.1 > /dev/null")
    verboseprint("done",file=sys.stderr)
    
    #valid_rowcols=np.where(nan_rowcols == False)[0]
    
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

def writeBedGraphFile(egv1,evr,header_rows,name,outfile,verbose=0):
    "write the compartment file"
    
    verboseprint = print if verbose else lambda *a, **k: None
    
    verboseprint("writing bed graph file (",outfile,") ... ",end="",file=sys.stderr)
    
    yBound=np.nanmax(abs(np.nanmin(egv1)),np.nanmax(egv1))
    yBound *= 1.25
    yBound=round(yBound,5)
    
    out_fh=open(outfile,"w")
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
    
def writeCompartmentFile(egv1,egv2,egv3,pca_score,geneDensity,header_rows,outfile,verbose=0):
    "write the compartment file"
    
    verboseprint = print if verbose else lambda *a, **k: None
    
    nan_geneDensity=np.sum(np.isnan(geneDensity))
    
    out_fh=open(outfile,"w")
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
    
def detectActiveCompartment(file,verbose=0):
    "detect the active compartment - overlap +/- with gene density"
    
    verboseprint = print if verbose else lambda *a, **k: None
    
    eigenMultiplier=1
    
    infh=open(file,'r')
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
    
    pca = decomposition.PCA(n_components=3)
    pca.fit(A)
    PCA(copy=True, n_components=3, whiten=False)
    pca_score=pca.explained_variance_ratio_
    pca_v = pca.components_
    
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
        return data,np.asarray(header_rows),np.asarray(header_cols)
    if(hrows):
        return data,np.asarray(header_rows)
    if(hcols):
        return data,np.asarray(header_cols)
    
    return data

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