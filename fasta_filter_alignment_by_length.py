#! /usr/bin/env python3

# load required modules
import argparse
import statistics
import re
from Bio import SeqIO

# create an ArgumentParser object ('parser') that will hold all the information necessary to parse the command line
parser = argparse.ArgumentParser(description="OUTPUT: sequence length, fraction of max sequence length in the file")

# add arguments
parser.add_argument('fasta', help='FASTA file')
parser.add_argument('-c', '--cutoff', help="Filter sequences with >X fraction characters (default = 0.8)", type=float, default=0.8)

args = parser.parse_args()

# open and parse FASTA file, storing sequence lengths
fasta_sequences = SeqIO.parse(open(args.fasta),'fasta')
for record in fasta_sequences:

	alignment_length = len(str(record.seq))

	gap_fraction = str(record.seq).count('-') / alignment_length

	if(gap_fraction > args.cutoff):
		continue
	else:
		print(record.format("fasta"), end = '')
