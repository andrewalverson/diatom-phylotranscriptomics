#! /usr/bin/env python3

# load required modules
import sys
import os
import re
import argparse
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

# create an ArgumentParser object ('parser') that will hold all the information necessary to parse the command line
parser = argparse.ArgumentParser(description="Extract CDS sequences from EMBL flat file")

# add arguments
parser.add_argument( "embl", help="name of EMBL flat file" )

args = parser.parse_args()

# open and parse the EMBL file
for record in SeqIO.parse(open(args.embl, "r"), "embl"):
	for feature in (record.features):
		# print(dir(feature))
		# print(feature.qualifiers)
		if feature.type == "CDS":
			# print(feature.qualifiers)
			# locus_tag = feature.qualifiers['locus_tag'][0]
			note = feature.qualifiers['note']
			id_match = "ID:(.*).CDS\d+"
			id = re.search(id_match, note[1]).group(1)
			# print(id)
			print( '>Pseudonitzschiamultistr_' + id.replace("_V1.4", "-V1.4", 1))
			print(feature.extract(record.seq))
			# print(dir(feature_seq))


