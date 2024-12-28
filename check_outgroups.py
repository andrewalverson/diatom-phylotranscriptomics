#! /usr/bin/env python3

# load required modules
import sys
import os
import re
import csv
import argparse

# create an ArgumentParser object ('parser') that will hold all the information necessary to parse the command line
parser = argparse.ArgumentParser(description="Summarize outgroup presence/absence in FASTA file")

# add arguments
parser.add_argument( "fasta", help="FASTA file" )
parser.add_argument( "out", help="File containing list (one species per row) of outgroup species" )

args = parser.parse_args()

# list containing names of outgroup taxa
outgroups = []

# get list of outgroup taxa
with open(args.out, 'r') as outgroup_file:
	for line in outgroup_file:
		line = line.rstrip()
		outgroups.append(line)

# total number of outgroups
num_outgroups = len(outgroups)
print("num_outgroups:", num_outgroups)

# print header line
s = "\t"
print(s.join(['header1', 'header2', 'header3' ]))

# open and parse FASTA file
fasta_sequences = SeqIO.parse(open(args.fasta),'fasta')

# open output file
for record in fasta_sequences:
    if record.description not in delete_seqs:
    	# print(record.id)
    	print(record.format("fasta"))
