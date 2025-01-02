#! /usr/bin/env python3

# load required modules
import argparse
import csv
import re
import os

# create an ArgumentParser object ('parser') that will hold all the information necessary to parse the command line
parser = argparse.ArgumentParser(description="Extracts simple summary data from batch BUSCO outfile")

# add arguments
parser.add_argument( "b", help="Summary outfile from batch BUSCO search" )

args = parser.parse_args()

print(','.join(['tip','complete','fragmented','missing','searched']))

# open and parse the GenBank taxonomy strings
with open(args.b, 'r') as t:
	for line in t:
		line = line.strip()
		f = re.compile("Input file is")
		file_match = re.search(f, line)
		if(file_match):
			path = line.split()
			head_tail = os.path.split(path[-1])
			tip = head_tail[-1].replace(".pep", "")
			print(tip, end = ',')

		b = re.compile("\|C:")
		busco_match = re.search(b, line)
		if(busco_match):
			c = re.search("C:(.*?)%", line)
			#print(line)
			complete = c.group(1)

			f = re.search("F:(.*?)%", line)
			fragmented = f.group(1)

			m = re.search("M:(.*?)%", line)
			missing = m.group(1)

			s = re.search("n:(\d+)", line)
			searched = s.group(1)

			print(','.join([complete,fragmented,missing,searched]))