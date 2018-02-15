#! /usr/bin/env python3

# this script is run in a directory that includes all of the assemblies that were multiplexed and run on the same Illumina lane
# it generates all the pairwise CD-HIT commands, which are designed to remove index crosstalk/barcode bleeding contigs from each sample

# about CD-HIT-EST-2D output files
# "the output are two files: a fasta file of sequences in db2 that are not similar to db1 and a text file that lists similar sequences between db1 & db2"

# load the system module
import sys
import os
import re

# get current working directory
cwd = os.getcwd()
    
# read in all the contents of the directory
directory_listing = os.listdir(cwd)

# loop over all the files in this directory, considering only *.fa files
for file1 in directory_listing:
    filename_re = '([RL]\d+)\_.*\.fa$'
   
    if re.search( filename_re, file1 ):
        print("######")
        for file2 in directory_listing:
            if re.search( filename_re, file2 ) and ( file1 != file2 ):
                outfile_pt1 = re.search( filename_re, file1 ).group(1)
                outfile_pt2 = re.search( filename_re, file2 ).group(1)
                outfile = outfile_pt1 + '-' + outfile_pt2 + '.cdhit'
                # print CD-HIT-EST-2D command (nucleotide sequences)
                print('cd-hit-est-2d', '-c 1 -T 0 -d 300 -aS 1', '-n 10', '-p 1', '-i', file2, '-i2 ', file1, '-o', outfile)

                
