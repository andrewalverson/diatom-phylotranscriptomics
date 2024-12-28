#! /usr/bin/env python3

import sys
import os
import re
import argparse
import statistics
from Bio import SeqIO

# create an ArgumentParser object ('parser') that will hold all the information necessary to parse the command line
parser = argparse.ArgumentParser()

# add arguments
parser.add_argument( "fasta", help="unaligned FASTA file" )

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
# indir = '/storage/aja/orthofinder/orthogroups/orthogroup_fasta_occupancy_10'
# indir = '/storage/aja/orthofinder/orthogroups/orthogroup_alignments_occupancy_10_round2/AA_unaligned'
indir = '/home/alverson/storage/orthogroups/orthogroup_alignments_occupancy_10_round2/AA_unaligned'

# set output directory
# outdir = '/storage/aja/orthofinder/orthogroups/orthogroup_alignments_occupancy_10'
# outdir = '/storage/aja/orthofinder/orthogroups/orthogroup_alignments_occupancy_10_round2/UPP_M_setting'
outdir = '/home/alverson/storage/orthogroups/orthogroup_alignments_occupancy_10_round2/UPP_M_setting'

# number of cpu's
num_cpu = 1

# parse name of input file in order to set prefix for output file names
a = args.fasta.split(".")
prefix = a[0]

# print command

print('# run UPP')
print('run_upp.py', '-m amino', '-x', str(num_cpu), '-s', '/'.join([indir, args.fasta]), '-d', outdir, '-o', prefix, '-M', upper_median, '2>', '.'.join([prefix, 'err']))
print()

print('# remove temporary files')
print('rm -r /tmp/sepp/' + prefix + '.*')
print()
