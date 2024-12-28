#! /usr/bin/env python3

# load required modules
import sys
import os
import csv
import argparse

# create an ArgumentParser object ('parser') that will hold all the information necessary to parse the command line
parser = argparse.ArgumentParser(description="Parse output of `fasta_alignment_statistics.py`, copying alignment files into directory bins according to some cutoff value")

# add arguments
parser.add_argument( "stats",   help="Output of `fasta_alignment_statistics.py`" )
parser.add_argument( "measure", help="[0] alignment, [1] length_with_gaps, [2] num_seqs, [3] min, [4] max, [5] mean, [6] stdev, [7] median, [8] median_of_upper_half", type=int)
parser.add_argument( "cutoff",  help="Cutoff for binning measure into two groups", type=int )

args = parser.parse_args()

alignment_dir = '/storage/aja/orthofinder/orthogroups/orthogroup_alignments_occupancy_10/AA_alignments_trees_round1/'
raxml_dir     = '/storage/aja/orthofinder/orthogroups/orthogroup_alignments_occupancy_10/AA_alignments_trees_round1/raxml/'
fasttree_dir  = '/storage/aja/orthofinder/orthogroups/orthogroup_alignments_occupancy_10/AA_alignments_trees_round1/fasttree/'

# [0] alignment, [1] length_with_gaps, [2] num_seqs, [3] min, [4] max, [5] mean, [6] stdev, [7] median, [8] median_of_upper_half

# open and parse stats file
with open(args.stats, 'r') as stats_file:
	for line in csv.reader(stats_file, delimiter = ','):
		# print(line[int(args.measure)])
		if(int(line[args.measure]) > args.cutoff):
			os.rename(alignment_dir+line[0], fasttree_dir+line[0])
			# print(line[args.measure], fasttree_dir+line[0])
		else:
			os.rename(alignment_dir+line[0], raxml_dir+line[0])
			# print(line[args.measure], raxml_dir+line[0])

stats_file.close()

