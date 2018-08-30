#! /usr/bin/env python3

# load required modules
import sys
import os
import re
import argparse

# create an ArgumentParser object ('parser') that will hold all the information necessary to parse the command line
parser = argparse.ArgumentParser(description="Show redundant query sequences in BLASTP output")

# add arguments
parser.add_argument( "blastp",  help="BLASTP of pep to translated CDS" )
parser.add_argument( "query_fasta",  help="FASTA file of the query sequences" )

args = parser.parse_args()

# a list of all the query sequences (first field in the BLAST output)
queries = []

# open and parse the BLASTP output (outfmt 6)
with open(args.blastp, 'r') as blp:
  for line in blp:
    line = line.rstrip()
    a = line.split()
    queries.append(a[0])
    
# open and parse FASTA file
fasta_sequences = SeqIO.parse(open(args.query_fasta),'fasta')
for record in fasta_sequences:
    if record.id in queries:
    	continue
    else:
        print(record.id)

