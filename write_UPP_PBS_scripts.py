#! /usr/bin/env python3

import sys
import os
import re
import argparse

# create an ArgumentParser object ('parser') that will hold all the information necessary to parse the command line
parser = argparse.ArgumentParser()

# add arguments
parser.add_argument( "fasta", help="Unaligned FASTA file" )

args = parser.parse_args()

# set walltime and number of cpu's
wall    = 6
num_cpu = 32

# set input directory
indir = '/storage/aja/orthofinder/orthogroups/orthogroup_fasta_occupancy_10'

# set output directory
outdir = '/storage/aja/orthofinder/orthogroups/orthogroup_alignments_occupancy_10'

# parse name of input file in order to set prefix for output file names
a = args.fasta.split(".")
prefix = a[0]

# print top of PBS script
print('#PBS -N', prefix+'.UPP')
print('#PBS -q q06h32c')
print('#PBS -j oe')
print('#PBS -o',  prefix+'.UPP.$PBS_JOBID')
print('#PBS -l nodes=1:ppn=32')
print('#PBS -l walltime=' + str(wall) + ':00:00')
print()

print('cd /scratch/$PBS_JOBID')
print()

print('# load modules')
print('module load gcc/7.2.1 python/3.6.0 java/sunjdk_1.8.0 sepp/4.3.5')
print()

print('# run UPP')
print('run_upp.py', '-m amino', '-x', str(num_cpu), '-s', '/'.join([indir, args.fasta]), '-d', outdir, '-o', prefix)

