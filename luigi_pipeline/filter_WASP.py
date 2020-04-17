#!/usr/bin/env python

import sys
import argparse
import getopt
import pysam
import os

def get_args():
    parser=argparse.ArgumentParser(description='filter BAM files')
    parser.add_argument('--BAMfile', '-in', type=str, required=True, help="BAM file")
    parser.add_argument('--outfile', '-out', type=str, required=True, help = 'Name for output file')

    args=parser.parse_args()
    return args

def remove_WASPAlignments(ReadTagsEntry, KeepInfo):
    keepRead = KeepInfo
    
    if 'vW' in ReadTagsEntry and ReadTagsEntry[1] > 1:
        keepRead=False
    
    return(keepRead)

def main():
    args = get_args()
    
    infile = pysam.Samfile(args.BAMfile, "rb")
    out = pysam.Samfile(args.outfile , "wb", template = infile )
    
    keep = True

    for DNAread in infile.fetch():
        intags = DNAread.tags 
        for entry in intags:
            keep = remove_WASPAlignments(entry, keep)

        if keep:
            out.write(DNAread)

        keep = True

    out.close()

if __name__ == '__main__':
    main()
