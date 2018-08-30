#! /usr/bin/env python3

# this script opens BLAST output, which has been filtered ahead of time to include only the query sequences that need to be reverse complemented
# a FASTA file, potentially containing a mix of properly and improperly oriented sequences, is opened, and reverse complements output for
# all the queries in the BLAST output file


# load required modules
import sys
import os
import re
import argparse
import csv
from Bio import SeqIO

# create an ArgumentParser object ('parser') that will hold all the information necessary to parse the command line
parser = argparse.ArgumentParser()

# add arguments
parser.add_argument( "blast", help="BLAST outfmt 6 with queries to be reverse complemented" )
parser.add_argument( "fasta", help="FASTA file containing sequences to be reverse complemented" )

args = parser.parse_args()

complement_these = []

# open and parse BLAST output
with open(args.blast, 'r') as blast_file:
    for line in csv.reader(blast_file, delimiter = '\t' ):
        complement_these.append(line[0])
blast_file.close()


# open and parse FASTA file
fasta_sequences = SeqIO.parse(open(args.fasta),'fasta')
for record in fasta_sequences:
    if record.id in complement_these:
        record = record.reverse_complement(id=True, name=True, description=True)
        print(record.format("fasta"), end = '')
    else:
        print(record.format("fasta"), end = '')

fasta_sequences.close()



