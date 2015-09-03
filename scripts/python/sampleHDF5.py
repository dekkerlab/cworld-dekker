
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
import shutil
import random
import copy
from collections import defaultdict
from datetime import datetime

verboseprint=lambda *a, **k: None
__version__ = "1.0"


def main():

    parser=argparse.ArgumentParser(description='Extract c-data from HDF5 file into TXT (matrix.gz)',formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('-i', '--input', dest='in_file', type=str, required=True, help='interaction matrix hdf5 file')
    parser.add_argument('-v', '--verbose', dest='verbose',  action='count', help='Increase verbosity (specify multiple times for more)')
    parser.add_argument('--info',dest='info', action='store_true', help='interaction matrix hdf5 file')
    parser.add_argument('-n',dest='num_reads', type=int, required=True,help='how many reads to sample from matrix')
    parser.add_argument('-o', '--output', dest='out_file', type=str, help='interaction matrix output file')
    parser.add_argument('-b','--blocksize', dest='blocksize', type=int, default=None, help='block size of HDF5 file')
    parser.add_argument('-p', dest='precision', type=int, default=4, help='output precision (# of digits)')
    parser.add_argument('--version', action='version', version='%(prog)s '+__version__)
    
    args=parser.parse_args()

    in_file=args.in_file
    verbose=args.verbose
    info=args.info
    num_reads=args.num_reads
    out_file=args.out_file
    blocksize=args.blocksize
    precision=args.precision
    
    log_level = logging.WARNING
    if verbose == 1:
        log_level = logging.INFO
    elif verbose >= 2:
        log_level = logging.DEBUG
    logging.basicConfig(level=log_level)
    
    verbose = info if info else verbose
    global verboseprint
    verboseprint = print if verbose else lambda *a, **k: None
    format_func=("{:."+str(precision)+"f}").format
    
    verboseprint("\n",end="")
    
    in_file_name=os.path.basename(in_file)
    inhdf=h5py.File(in_file,'r')
    
    # attrs
    genome=inhdf.attrs['genome'][:]
    # datasets
    bin_positions=inhdf['bin_positions'][:]
    chr_bin_range=inhdf['chr_bin_range'][:]
    chrs=list(inhdf['chrs'][:])
    # matrix shape
    dim=inhdf['interactions'].shape
    
    # ensure symmetrical
    ensure_symmetrical(dim)
    nrow,ncol=dim
    n=nrow=ncol
    
    verboseprint("num_reads",num_reads)
    
    # calculate optimal block size
    hdf_blocksize=inhdf['interactions'].chunks[0]
    blocksize=get_blocksize(hdf_blocksize,blocksize)
    
    # get _optional_ y axis headers
    headers = np.empty(ncol)
    headers[:] = np.nan
    if "headers" in inhdf.keys():
        headers=inhdf['headers'][:]
    else:
        headers=np.array([str(i)+'|'+genome+'|'+str(chrs[bin_positions[i,0]])+':'+str(bin_positions[i,1])+'-'+str(bin_positions[i,2]) for i in np.arange(n)])
    
    # build chr lookup dict
    chr_dict={}
    for i,c in enumerate(chrs):
        chr_dict[c]=i
   
    B=np.zeros((nrow,ncol))
    
    sampleList=[]
    for i in xrange(0,n,blocksize):
        current_block=inhdf['interactions'][i:i+blocksize,:]
        for y in xrange(current_block.shape[0]):
            for x in xrange(current_block.shape[1]):
                count=current_block[y][x]
                sampleList.append((count,(i+y,x)))
                
    last_fib_i=0
    for fib_n,fib_i in enumerate(fib(num_reads)):
        
        print(fib_n,fib_i)
        
        out_file="sample"+str(fib_n)+".hdf5"
        verboseprint(out_file)
        
        verboseprint("copying hdf file")
        shutil.copy(in_file,out_file)
        
        outhdf=h5py.File(out_file)
            
        verboseprint("expanding matrix into list")
        y_offset=0
            
        verboseprint("building output matrix")
        for i in weighted_sample(sampleList,fib_i-last_fib_i):
            x,y=i
            B[y][x] += 1
            B[x][y] += 1
        
        verboseprint("writing matrix to hdf")    
        for i in xrange(0,n,blocksize):
            outhdf['interactions'][i:i+blocksize,:]=B[i:i+blocksize,:]
          
        outhdf.close()
    
        last_fib_i=fib_i
        verboseprint("")
    
    inhdf.close()
    
def fib(n):
    a,b = 1,1
    
    while(b < n):
        a,b = b,a+b
        yield a
    yield n
    
def weighted_sample(items, n):
    total = float(sum(w for w, v in items))
    i = 0
    w, v = items[0]
    while n:
        x = total * (1 - random.random() ** (1.0 / n))
        total -= x
        while x > w:
            x -= w
            i += 1
            w, v = items[i]
        w -= x
        yield v
        n -= 1
        
def input_wrapper(in_file):
    if in_file.endswith('.gz'):
        fh=gzip.open(in_file,'r')
    else:
        fh=open(in_file,'r')
        
    return fh
    
def output_wrapper(out_file):
    
    if out_file.endswith('.gz'):
        fh=gzip.open(out_file,'wb')
    else:
        fh=open(out_file,'w')
    
    suppress_comments=0
    
    # disable comment(s)if (UCSC format file)
    if out_file.endswith('.bed'):
        suppress_comments = 1
    if out_file.endswith('.bed.gz'):
        suppress_comments = 1
    if out_file.endswith('.bedGraph'):
        suppress_comments = 1
    if out_file.endswith('.bedGraph.gz'):
        suppress_comments = 1
    if out_file.endswith('.wig'):
        suppress_comments = 1
    if out_file.endswith('.wig.gz'):
        suppress_comments = 1

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
    
def write_factors(out_file,x_bin_mask,y_bin_mask,chrs,bin_positions,factors):
    """output x/y axis ICE factors (after chr/zoom subset)
    """
    
    out_fh=output_wrapper(out_file+'.xfactors')
    print("# binIndex\tbinChr\tbinStart\tbinEnd",file=out_fh)
    for i in np.nonzero(x_bin_mask)[0]:
        print(str(i)+"\t"+chrs[bin_positions[i,0]]+"\t"+str(bin_positions[i,1])+"\t"+str(bin_positions[i,2])+"\t"+str(factors[i]),file=out_fh)
    out_fh.close()
    
    out_fh=output_wrapper(out_file+'.yfactors')
    print("# binIndex\tbinChr\tbinStart\tbinEnd",file=out_fh)
    for i in np.nonzero(y_bin_mask)[0]:
        print(str(i)+"\t"+chrs[bin_positions[i,0]]+"\t"+str(bin_positions[i,1])+"\t"+str(bin_positions[i,2])+"\t"+str(factors[i]),file=out_fh)
    out_fh.close()    
           
def write_bins(out_file,x_bin_mask,y_bin_mask,chrs,bin_positions):
    """output x/y axis bins (after chr/zoom subset)
    """
    
    out_fh=output_wrapper(out_file+'.xbins')
    print("# binIndex\tbinChr\tbinStart\tbinEnd",file=out_fh)
    for i in np.nonzero(x_bin_mask)[0]:
        print(str(i)+"\t"+chrs[bin_positions[i,0]]+"\t"+str(bin_positions[i,1])+"\t"+str(bin_positions[i,2]),file=out_fh)
    out_fh.close()
    
    out_fh=output_wrapper(out_file+'.ybins')
    print("# binIndex\tbinChr\tbinStart\tbinEnd",file=out_fh)
    for i in np.nonzero(y_bin_mask)[0]:
        print(str(i)+"\t"+chrs[bin_positions[i,0]]+"\t"+str(bin_positions[i,1])+"\t"+str(bin_positions[i,2]),file=out_fh)
    out_fh.close()    
            
def build_bin_mask(n,chrs,zoom_dict,chr_dict,chr_bin_range,bin_positions,axis=None):
    """build a 1D mask (x or y axis) based upon user chr/coor selection
    """
    
    bin_mask=np.zeros(n,dtype=bool)
   
    # build bin mask based on user chr/zoom selection
    for c in chrs:
        c_ind=chr_dict[c]
        r=chr_bin_range[chr_dict[c]]
        if c in zoom_dict:
            zoom_coord_arr=zoom_dict[c]
            for zoom_coord in zoom_coord_arr:
                tmp_bin_positions=bin_positions[r[0]:r[1]+1]
                for i,b in enumerate(tmp_bin_positions):
                    if b[2] < zoom_coord[1]: continue
                    if b[1] > zoom_coord[2]: break
                    overlap=is_overlap([zoom_coord[1],zoom_coord[2]], [b[1],b[2]])
                    if(overlap > 0):
                        bin_mask[r[0]+i]=True
        else:
            bin_mask[r[0]:r[1]+1]=True
    verboseprint("\t",axis," bin_mask\t",np.sum(bin_mask),sep="")
    
    return(bin_mask)
    
def dump_hdf_info(in_file,in_file_name,nrow,ncol,genome,hdf_blocksize,blocksize,chrs,chr_dict,chr_bin_range,bin_positions):
    """dump hdf info
    """
    
    verboseprint("inputFile",in_file,sep="\t")
    verboseprint("inputFileName",in_file_name,sep="\t")
    verboseprint("matrix shape\t",nrow," x ",ncol,sep="")
    verboseprint("assembly",genome,sep="\t")
    verboseprint("h5 chunk",hdf_blocksize,sep="\t")
    verboseprint("user chunk",blocksize,sep="\t")
    verboseprint("\nchrs",sep="\t")
    
    for i,c in enumerate(chrs):
        cbr=chr_bin_range[chr_dict[c]]
        start,end=bin_positions[cbr[0]][1],bin_positions[cbr[1]][2]
        size=(end-start)+1
        nbins=(cbr[1]-cbr[0])+1
        verboseprint("\t",i,"\t",c,":",start,"-",end,"\t(",size,")\t",cbr,"\t",nbins,sep="")
        
    verboseprint("")
    quit()
    
    
def get_blocksize(hdf_blocksize,blocksize):
    """adjust blocksize to be evenly divisible by hdf_blocksize
    """
    
    if blocksize == None:
        blocksize = hdf_blocksize
    else:
        if blocksize%hdf_blocksize != 0:
            blocksize=int(math.ceil(blocksize/hdf_blocksize)*hdf_blocksize)
            
    verboseprint("hdf_blocksize",hdf_blocksize)
    verboseprint("blocksize",blocksize)
    verboseprint("")
    
    return(blocksize)
    
def ensure_symmetrical(dim):
    """ensure nrow=ncol [symmetrical]
    """
    nrow,ncol=dim
    
    if nrow!=ncol:
        sys.exit('error: non-symmetrical matrix found!')
        
def subset_by_coords(zoom_chrs,zoom_dict,coord):
    """read UCSC coordinates, extract chr, coordinates 
    """
    
    # process zoom coordinates
    if(coord!=None):
        for z in coord:
            coord=split_coord(z)
            
            if coord==None:
                verboseprint("invalid coord",z)
                continue
                
            coord_chr,coord_start,coord_end=coord
            if coord_chr not in zoom_chrs:
                zoom_chrs += [coord_chr]
            zoom_dict[coord_chr].append(coord)
            
    return zoom_chrs,zoom_dict
        
def subset_by_bed(bed_chrs,bed_dict,bed_file,element_exten):
    """read bed file, extract chr, coordinates 
    """
    
    num_elements=0
    for b in bed_file:
        e_fh=input_wrapper(b)
        for i,li in enumerate(e_fh):
            li=li.rstrip("\n")
            if li.startswith("#") or li.startswith("track"):
                continue
            
            lineList=li.split("\t")
            lineList[1]=max(1,(int(lineList[1])-element_exten))
            lineList[2]=(int(lineList[2])+element_exten)
            z=lineList[0]+':'+str(lineList[1])+'-'+str(lineList[2])
            
            bed_coord=split_coord(z)
            
            if bed_coord==None:
                verboseprint("invalid coord",z)
                continue
                
            bed_chr,bed_start,bed_end=bed_coord
            if bed_chr not in bed_chrs:
                bed_chrs += [bed_chr]
            bed_dict[bed_chr].append(bed_coord)
            num_elements += 1
        e_fh.close()
        
    return bed_chrs,bed_dict
        
def is_float(s):
    try:
        float(s)
        return True
    except ValueError:
        return False
        
def getSmallUniqueString():  
    tmp_uniq=str(uuid.uuid4())
    tmp_uniq=tmp_uniq.split('-')[-1]
    return(tmp_uniq)
    
def bin2header(bin,genome,chrs,index=getSmallUniqueString()):
    #name|assembly|chr:start-end
    header=str(index)+'|'+genome+'|'+str(chrs[bin[0]])+':'+str(bin[1])+'-'+str(bin[2])
    return(header)

def deGroupChr(chr_id):
    return(chr_id.split('-')[0])
    
def deGroupHeader(header,extractBy="liteChr",index=getSmallUniqueString()):
    m=re.search(r'(\S+)\|(\S+)\|(\S+):(\d+)-(\d+)',header)
    if m==None:
        sys.exit('error: incorrect input format!')
                
    bin_id,genome,chr_id,bin_start,bin_end=m.groups()
    chr_id=chr_id.split('-')[0]

    header=str(bin_id)+'|'+genome+'|'+str(chr_id)+':'+str(bin_start)+'-'+str(bin_end)
    
    return(header)
    
def split_coord(z):
    """validate and split zoom coordinate.
    coordinates must be UCSC formatted.
    e.g. chr1:500-1000
    chr(colon)start(hyphen)end where start <= end
    """
    z=z.replace(',','')
    zoom_coord=re.search(r'(\S+):(\d+)-(\d+)',z)
    
    if zoom_coord==None:
        return None
        
    zoom_chr,zoom_start,zoom_end=zoom_coord.groups()
    zoom_start=int(zoom_start)
    zoom_end=int(zoom_end)
    
    if(zoom_start > zoom_end):
        return None
        
    return [zoom_chr,zoom_start,zoom_end]
        
def de_dupe_list(input):
    """de-dupe a list, preserving order.
    """
    
    output = []
    for x in input:
        if x not in output:
            output.append(x)
    return output

def byte_to_megabyte(byte):
    """convert bytes into megabytes.
    """
    
    return round(((byte / 1000) / 1000),4) # megabyte
    # return float(((byte / 1024) / 1024),4) # mebibyte

    
def flip_intervals(a,b):
    """flip intervals, to ensure a < b
    """
    
    return(b,a)
    
def is_overlap(a, b):
    """test to for overlap between two intervals.
    """
    
    if(a[0] > a[1]):
        sys.exit('\nerror: incorrectly formated interval! start '+str(a[0])+' > end '+str(a[1])+'!\n\t'+str(a)+' '+str(b)+'\n')
    if(b[0] > b[1]):
        sys.exit('\nerror: incorrectly formated interval! start '+str(b[0])+' > end '+str(b[1])+'!\n\t'+str(a)+' '+str(b)+'\n')
    
    if a[0] < b[0] and a[1] > b[1]:
        return((b[1]-b[0])+1)
    
    if b[0] < a[0] and b[1] > a[1]:   
        return((a[1]-a[0])+1)
        
    if b[0] < a[0]:
        a,b=flip_intervals(a,b)
           
    return max(0, ( min(a[1],b[1]) - max(a[0],b[0]) ) ) 
    
if __name__=="__main__":
    main()

   