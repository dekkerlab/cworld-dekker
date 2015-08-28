
from __future__ import print_function

import numpy as np
import scipy as sp
import cv2
import sys
import argparse
import logging
import gzip
import re
import os
import math
import matplotlib.pyplot as plt

def main():

    parser=argparse.ArgumentParser(description='Extract c-data from HDF5 file into TXT (matrix.gz)',formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('-i', '--input', dest='infile', type=str, required=True, help='input image')
    parser.add_argument('-v', '--verbose', dest='verbose',  action='count', help='Increase verbosity (specify multiple times for more)')
    parser.add_argument('--version', action='version', version='%(prog)s 1.0')

    args=parser.parse_args()

    infile=args.infile
    verbose=args.verbose
    
    log_level = logging.WARNING
    if verbose == 1:
        log_level = logging.INFO
    elif verbose >= 2:
        log_level = logging.DEBUG
    logging.basicConfig(level=log_level)
    
    verboseprint = print if verbose else lambda *a, **k: None
    
    verboseprint("\n",end="")
    
    infile_name=os.path.basename(infile)
    if output_relative:
       infile_name=infile
       
    img = cv2.imread(infile)

    kernel = np.ones((5,5),np.float32)/25
    dst = cv2.filter2D(img,-1,kernel)

    plt.subplot(121),plt.imshow(img),plt.title('Original')
    plt.xticks([]), plt.yticks([])
    plt.subplot(122),plt.imshow(dst),plt.title('Averaging')
    plt.xticks([]), plt.yticks([])
    plt.show()

if __name__=="__main__":
    main()