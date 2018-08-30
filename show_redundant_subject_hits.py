#! /usr/bin/env python3

# load required modules
import sys
import os
import re
import argparse

# create an ArgumentParser object ('parser') that will hold all the information necessary to parse the command line
parser = argparse.ArgumentParser(description="Show redundant subject hits in BLASTP output")

# add arguments
parser.add_argument( "blastp",  help="BLASTP of pep to translated CDS" )

args = parser.parse_args()

# key = header in CDS file, value = corresponding header in pep file
subjects = []

# open and parse the BLASTP output (outfmt 6)
with open(args.blastp, 'r') as blp:
  for line in blp:
    line = line.rstrip()
    a = line.split()
    subjects.append(a[1])
    

for i in subjects:
    if subjects.count(i) > 1:
      print(i)

