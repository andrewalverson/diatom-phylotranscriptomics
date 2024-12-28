#! /usr/bin/env python3

import sys
import os
import re
import argparse
import statistics
from Bio import SeqIO

# create an ArgumentParser object ('parser') that will hold all the information necessary to parse the command line
parser = argparse.ArgumentParser(description = "write a PBS job script to run UPP with the '-M' setting")

# add arguments
parser.add_argument( "fasta", help="unaligned FASTA file" )
parser.add_argument( "-q", "--queue", help="trestles queue, default = aja", default = 'aja' )
parser.add_argument( "-w", "--wall",  help="walltime", default = '6' )

args = parser.parse_args()

### compute median sequence length for this dataset ###

lengths = [] # list of sequence lengths
upper_half_lengths = [] # list of sequences lengths >= median sequence length
median_length = 0 # store median sequence length
upper_median  = 0 # calculate the median length for all sequence, then calculate the median for all sequences >= the global median length

# open and parse FASTA file, storing sequence lengths
fasta_sequences = SeqIO.parse(open(args.fasta),'fasta')
for record in fasta_sequences:
    lengths.append(len(record.seq))
fasta_sequences.close()

median_length = int(statistics.median(lengths))

for i in lengths:
	if i >= median_length:
		upper_half_lengths.append(i)

upper_median = int(statistics.median(upper_half_lengths))

### write UPP PBS script for this dataset ###

# set input directory
indir = '/storage/aja/orthofinder/orthogroup_fasta_occupancy_20_no_stops'

# set output directory
outdir = '/storage/aja/orthofinder/gene_trees_round_1/upp'

# number of cpus
num_cpu = 8

# parse name of input file in order to set prefix for output file names
a = args.fasta.split(".")
prefix = a[0]

# print top of PBS script
print('#PBS -N', prefix+'.UPP')
print('#PBS -q', args.queue)
print('#PBS -j oe')
print('#PBS -o',  prefix+'.UPP.$PBS_JOBID')
print('#PBS -l nodes=1:ppn=' + str(num_cpu))
print('#PBS -l walltime=' + str(args.wall) + ':00:00')
print()

print('TMPDIR=/scratch/$PBS_JOBID')
print('cd /scratch/$PBS_JOBID')
print()

print('# load modules')
print('module load gcc/7.2.1 java/sunjdk_1.8.0 python/3.7.3-anaconda-razor')
print()

print('conda init bash')
print('source activate /home/aja/.conda/envs/sepp')


print()
print('# run UPP')
print('run_upp.py', '-m amino', '-x', str(num_cpu), '-s', '/'.join([indir, args.fasta]), '-d', outdir, '-o', prefix, '-M', upper_median)
