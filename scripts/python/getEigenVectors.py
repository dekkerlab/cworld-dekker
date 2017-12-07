#!/usr/local/bin/python
"""
***********************************************
- PROGRAM: getEigenVectors.py
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

# deprecated from scipy and unusde in the script:
# from scipy.stats.stats import nanmean
import numpy as np
import scipy as sp
from  collections import *

# user defined modules
from noamlib  import *

# For eigenvectors and eigenvalues
from scipy import linalg as la
from math import cos,log,sin,sqrt 
# deprecated from scipy to be replaced with Cython:
from scipy import weave 

def main():
    print("")
    
    # Get input options
    args = check_options()

    # Store the variables
    inputMatrix = args.inputMatrix
        
    if not os.path.isfile(inputMatrix):
        sys.exit('invalid input file! (non-existant)')
        
    print("inputMatrix",inputMatrix)
    inputMatrix_name=os.path.basename(inputMatrix)
    inputMatrix_name=re.sub(".gz", "", inputMatrix_name)    
    inputMatrix_name=re.sub(".matrix", "", inputMatrix_name)
    print("inputMatrix_name",inputMatrix_name)
    
    print("loading matrix ... ",end="")
    if inputMatrix.endswith('.gz'):
        infh=gzip.open(inputMatrix,'r')
    else:
        infh=open(inputMatrix,'r')
   
    matrix,header_rows,header_cols = load_matrix((l for l in infh if not l.startswith('#')), hrows=1, hcols=1) # since this returns data, header_rows and header_cols
    infh.close()
    print("done")
    
    print("")
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

    print("matrix assembly:",assembly)
    refSeqDir="/cShare/tools/geneDensity/"
    refSeqFile=refSeqDir+assembly+".refseq.txt"
    
    try:
        with open(refSeqFile) as rsfh:
            pass
    except IOError as e:
        sys.exit("invalud refSeq file! ("+refSeqFile+")")

    print("refSeqFile",refSeqFile)
    
    print("")
    # get number of rows/col (assuming symmetrical)
    nrows=matrix.shape[0]
    
    # find nan rows
    print("finding nan rows ... ",end="")
    nan_rowcols=np.sum(np.isnan(matrix),0)==matrix.shape[0]
    # remove nan rows
    matrix=matrix[-nan_rowcols,:][:,-nan_rowcols]
    print("done")
    
    # convert all nan to 0
    print("converting all 2D nan to 0 ... ",end="")
    matrix = np.nan_to_num(matrix)
    print("done")
    
    # calculate obs/exp matrix
    print("calculating obs/exp ... ",end="")
    matrix = observedOverExpected(matrix)
    print("done")
    
    # calculate corrcoef matrix
    print("calculating coorcoef ... ",end="")
    matrix = np.corrcoef(matrix)
    print("done")
    
    # do eigenvector analysis
    print("running PCA ... ",end="")
    egv_data = calculate_eigen(matrix, 3)
    print("done")
    
    eig1=egv_data[0][0]
    eig2=egv_data[0][1]
    eig3=egv_data[0][2]
 
    # pre-populate with nan
    egv1=np.nan*np.ones(nrows)
    egv2=np.nan*np.ones(nrows)
    egv3=np.nan*np.ones(nrows)
    
    egv1[~nan_rowcols]=eig1
    egv2[~nan_rowcols]=eig2
    egv3[~nan_rowcols]=eig3

    geneDensity=np.nan*np.ones(nrows)
    
    compartmentFile=inputMatrix_name+".compartments"    
    writeCompartmentFile(egv1,egv2,egv3,geneDensity,header_rows,compartmentFile)
    
    print("intersecing compartments with ref seq ... ",end="")
    compartmentRefSeqFile=compartmentFile+".refSeq.txt"
    os.system("bedtools intersect -a "+compartmentFile+" -b "+refSeqFile+" -c > "+compartmentRefSeqFile)
    print("done")
    
    eigenMultiplier,geneDensity = detectActiveCompartment(compartmentRefSeqFile)
    os.system("rm "+compartmentRefSeqFile)

    print("flipping vectors by",eigenMultiplier," ... ",end="")
    egv1 *= eigenMultiplier
    egv2 *= eigenMultiplier
    egv3 *= eigenMultiplier
    print("done")
    
    writeCompartmentFile(egv1,egv2,egv3,geneDensity,header_rows,compartmentFile)
    
    eig1BedGraphFile=inputMatrix_name+".eigen1.bedGraph"    
    writeBedGraphFile(egv1,header_rows,inputMatrix_name,eig1BedGraphFile)
    
    print("drawing plot (",inputMatrix_name,") ... ",end="")
    drawPlot = "/cShare/tools/Rscripts/plotEigen.R"
    os.system("Rscript "+drawPlot+" `pwd` "+compartmentFile+" "+inputMatrix_name+" 0.1 > /dev/null")
    print("done")
    
    print("")

def writeBedGraphFile(egv1,header_rows,name,outfile):
    "write the compartment file"
    
    print("writing bed graph file (",outfile,") ... ",end="")
    
    yBound=np.nanmax(abs(np.nanmin(egv1)),np.nanmax(egv1))
    yBound *= 1.25
    yBound=round(yBound,5)
    
    out_fh=open(outfile,"w")
    print("track type=bedGraph name='"+name+"' description='"+name+"' visibility=full autoScale=off viewLimits="+str(-yBound)+":"+str(yBound)+" color=0,0,0 altColor=100,100,100",end="\n",file=out_fh)

    for i,header in enumerate(header_rows):
        m=re.search(r'(\S+)\|(\S+)\|(\S+):(\d+)-(\d+)',header)
        if m==None:
            sys.exit('error: incorrect input format!')  

        bin_id,genome,chr_id,bin_start,bin_end=m.groups()
        
        chr_id=chr_id.split("-")[0]
        
        eigen1=round(egv1[i],5)
        
        print(str(chr_id)+"\t"+str(bin_start)+"\t"+str(bin_end)+"\t"+str(eigen1),end="\n",file=out_fh)
  
    out_fh.close()
    print("done")

  
def writeCompartmentFile(egv1,egv2,egv3,geneDensity,header_rows,outfile):
    "write the compartment file"
    
    nan_geneDensity=np.sum(np.isnan(geneDensity))
    
    out_fh=open(outfile,"w")
    print("writing eigenvector file (",outfile,") ... ",end="")
    
    if len(geneDensity)==nan_geneDensity:
        print("#chr\tstart\tend\tname\tindex\teigen1\teigen2\teigen3",end="\n",file=out_fh)
    else:
        print("#chr\tstart\tend\tname\tindex\teigen1\teigen2\teigen3\tgeneDensity",end="\n",file=out_fh)
        
    for i,header in enumerate(header_rows):
        m=re.search(r'(\S+)\|(\S+)\|(\S+):(\d+)-(\d+)',header)
        if m==None:
            sys.exit('error: incorrect input format!')  

        bin_id,genome,chr_id,bin_start,bin_end=m.groups()
        
        eigen1=egv1[i]
        eigen2=egv2[i]
        eigen3=egv3[i]
        
        if len(geneDensity)==nan_geneDensity:
            print(str(chr_id)+"\t"+str(bin_start)+"\t"+str(bin_end)+"\t"+str(header)+"\t"+str(i)+"\t"+str(eigen1)+"\t"+str(eigen2)+"\t"+str(eigen3),end="\n",file=out_fh)
        else:
            nGenes=geneDensity[i]
            print(str(chr_id)+"\t"+str(bin_start)+"\t"+str(bin_end)+"\t"+str(header)+"\t"+str(i)+"\t"+str(eigen1)+"\t"+str(eigen2)+"\t"+str(eigen3)+"\t",str(nGenes),end="\n",file=out_fh)
    
    out_fh.close()
    print("done")
    
def detectActiveCompartment(file):
    "detect the active compartment - overlap +/- with gene density"
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
        
        nGenes=float(a[8])
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
    
    if negAvg > posAvg:
        eigenMultiplier=-1
        
    print("\tposSum",posSum,"posCount",posCount,"posAvg",posAvg)
    print("\tnegSum",negSum,"negCount",negCount,"negAvg",negAvg)
    print("\teigenMultiplier",eigenMultiplier)
    
    #print("\tposAvg",posAvg)
    #print("\tnegAvg",negAvg)
    #print("\teigenMultiplier",eigenMultiplier)
    
    return eigenMultiplier,geneDensity
    
def calculate_eigen(A, numPCs = 3):
    """performs eigen vector analysis, and returns 3 best principal components
    result[0] is the first PC, etc"""    
    A = np.array(A,float)
    M = (A-np.mean(A.T,axis=1)).T 
    covM = np.dot(M,M.T)
    [latent,coeff] =  scipy.sparse.linalg.eigsh(covM,numPCs)
    return (np.transpose(coeff[:,::-1]),latent[::-1])

def observedOverExpected(matrix):
    "Calculates observedOverExpected of any contact map"
    data = np.asarray(matrix, dtype = float, order = "C")
    N = data.shape[0]
    bins = logbins(1,N,1.2)
    bins = [(0,1)] + [(bins[i],bins[i+1]) for i in xrange(len(bins)-1)]
    bins = np.array(bins,order = "C")
    M = len(bins)
    code = r"""
    #line 50 "binary_search.py"
    using namespace std;
    for (int bin = 0; bin < M; bin++)
    {
        int start = bins[2 * bin];
        int end = bins[2 * bin + 1];
        
        double ss = 0 ;
        int count   = 0 ;  
        for (int offset = start; offset < end; offset ++)
        {
            for (int j = 0; j < N - offset; j++)
            {
                ss += data[(offset + j) * N + j];
                count += 1;                            
            }
        }
        double meanss = ss / count;
        //printf("%lf\n",meanss); 
        for (int offset = start; offset < end; offset ++)
        {
            for (int j = 0; j < N - offset; j++)
            {
                if (meanss !=  0)
                {
                    data[(offset + j) * N + j] /= meanss;                                                             
                    if (offset > 0) {data[(offset + j)  + j*N] /= meanss;}
                }
            }
        }
    }
    
    """
    support = """
    #include <math.h>  
    """
    weave.inline(code, ['data', 'bins' , 'N' ,'M'], extra_compile_args=['-march=native -malign-double -O3'],support_code =support )
    return data

def logbins(a, b, pace, N_in=0):
    "create log-spaced bins"
    a = int(a)
    b = int(b) 
    beg = log(a)
    end = log(b - 1)
    pace = log(pace)
    N = int((end - beg) / pace)
     
    if N > 0.8 * (b-a): 
        return np.arange(a,b+1)
    
    if N_in != 0: N = N_in  
    pace = (end - beg) / N
    mas = np.arange(beg, end + 0.000000001, pace)
    ret = np.exp(mas)
    ret = np.array([int(i) for i in ret])
    ret[-1] = b 
    for i in xrange(len(ret) - 1):
        if ret[i + 1] <= ret[i]:
            ret[i + 1: - 1] += 1
    return [int(i) for i in ret]


def check_options():
    ''' Checks the options to the program '''

    # Create parser object
    parser = argparse.ArgumentParser(add_help=True,
        formatter_class=argparse.RawDescriptionHelpFormatter)

    # Add arguments 
    parser.add_argument('-i' , metavar='--inputMatrix'  , help="*Input matrix file", dest="inputMatrix", type=str, default="")
    
    # Parse command line with parse_args and store it in an object
    args = parser.parse_args()

    return args

if __name__=="__main__":
      main()
