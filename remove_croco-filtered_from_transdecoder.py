#! /usr/bin/env python3

# load required modules
import sys
import os
import re
import argparse
from Bio import SeqIO

# create an ArgumentParser object ('parser') that will hold all the information necessary to parse the command line
parser = argparse.ArgumentParser(description="Remove CroCo-identified crosstalk (FASTA format) from Transdecoder FASTA file")

# add arguments
parser.add_argument( "trans", help="FASTA format; sequences in trans.fasta that are not in croco.fasta will be written to STDOUT" )
parser.add_argument( "croco", help="FASTA format; sequences in trans.fasta that are not in croco.fasta will be written to STDOUT" )

args = parser.parse_args()

croco_descriptions = []

# open and parse Transdecoder FASTA file
trans_sequences = SeqIO.parse(open(args.trans),'fasta')

# open and parse CroCo FASTA file; sequences in this file will be removed from the Transdecoder file
croco_sequences = SeqIO.parse(open(args.croco),'fasta')

# SeqIO does not implement comparisons between 2 sets of FASTA files
# so add seq descriptions to a list, fasta2_descriptions, then compare the two files that way
for record in croco_sequences:
	croco_descriptions.append(record.description)

# compare the two sets of sequences
for trans_record in trans_sequences:
	match = re.search("([A-Z].*)(_m\.\d+)*", trans_record.id)
	seq_base_name = match.group(1)
	
	if seq_base_name in croco_descriptions:
		continue
	else:
		print(trans_record.format("fasta"), end = '')

