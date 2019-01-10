#! /usr/bin/env python3

# load required modules
import sys
import os
import re
import argparse
from Bio import SeqIO

# create an ArgumentParser object ('parser') that will hold all the information necessary to parse the command line
parser = argparse.ArgumentParser(description="Selects the longest translation from Transdecoder output (FASTA format)")

# add arguments
parser.add_argument( "fasta", help="Transdecoder FASTA file" )

args = parser.parse_args()

longest_orfs        = {}
longest_orf_lengths = {}

# open and parse the Transdecoder FASTA sequences
fasta_sequences = SeqIO.parse(open(args.fasta),'fasta')

# open and parse FASTA file
with open(args.fasta, 'r') as fasta_file:

	for line in fasta_file:
		line = line.rstrip()

		# is this a header line?
		if re.search("^>", line):
			line = line[1:]

			a = line.split()

			# get the sequence ID and remove the suffix added by Transdecoder
			id = a[0]
			id_basename_re = "(.*)\|.*$"
			match = re.search(id_basename_re, id)
			id = match.group(1)
			# print(id)

			# get the length of this ORF
			b = a[6].split(":")
			orf_length = int(b[1])
			# print(orf_length)

			# check whether this id exists in the dictionary of orf lengths
			# if it does, compare the lengths and keep the longest of the two
			# if it doesn't, store in in longest_orfs dictionary
			if id in longest_orf_lengths:
				# print(line)
				if orf_length > int(longest_orf_lengths[id]):
					# store the length of this orf
					longest_orf_lengths[id] = orf_length
					# store the fasta header (description) of this orf
					longest_orfs[id] = line
				else:
					continue
			else:
				# store the length of this orf
				longest_orf_lengths[id] = orf_length
				# store the fasta header (description) of this orf
				longest_orfs[id] = line

		# continue if not a header line (sequence line)
		else:
			continue


# put the sequence descriptions of the longest orfs in a list for more easy access
sequence_descriptions = list(longest_orfs.values())

# loop over the FASTA sequences, printing the longest ORF for each Transdecoder input sequence
for record in fasta_sequences:
	if record.description in sequence_descriptions:
		print(record.format("fasta"), end = '')


