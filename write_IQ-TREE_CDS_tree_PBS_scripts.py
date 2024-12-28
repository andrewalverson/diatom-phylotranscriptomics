#! /usr/bin/env python3

import sys
import os
import re
import argparse
import statistics
from Bio import SeqIO


def get_args():
	# create an ArgumentParser object ('parser') that will hold all the information necessary to parse the command line
	parser = argparse.ArgumentParser()

	# add arguments
	parser.add_argument( "fasta", help="FASTA codon alignment" )
	parser.add_argument( "--queue", "-q", help="AHPCC queue", default="q06h32c" )
	parser.add_argument( "--wall",  "-w", help="walltime", type=int, default=6 )

	return parser.parse_args()


def parse_filename():

	# parse name of input file in order to set prefix for output file names
	a = args.fasta.split("_")
	
	return a[0]


def get_alignment_length():

	# alignment length
	alignment_length = 0

	# initialize a counter
	counter = 0

	# open and parse FASTA file, storing sequence lengths
	fasta_sequences = SeqIO.parse(open(args.fasta),'fasta')
	for record in fasta_sequences:
		if counter > 0:
			continue

		else:
		    alignment_length = (len(record.seq))
		    counter += 1
	fasta_sequences.close()

	return alignment_length


def write_pbs(prefix):

	# create and open PBS output file
	pbs = open(prefix + '_iqtree.pbs', "w")

	pbs.write('#PBS -N ' + prefix + '.iqtree-cds' + '\n')
	pbs.write('#PBS -q ' + args.queue + '\n')
	pbs.write('#PBS -j oe ' + '\n')
	pbs.write('#PBS -o ' + prefix + '.iqtree-cds.$PBS_JOBID' + '\n')
	pbs.write('#PBS -l nodes=1:ppn=16' + '\n')
	pbs.write('#PBS -l walltime=' + str(args.wall) + ':00:00' + '\n')
	pbs.write('\n')

	pbs.write('cd $PBS_O_WORKDIR' + '\n')
	pbs.write('\n')

	pbs.write('# run IQ-TREE' + '\n')
	pbs.write('/home/aja/local/bin/iqtree -s ' + args.fasta + ' -pre ' + prefix + '_iqtree' + ' -spp ' + prefix +  '_iqtree.partition' + ' -st DNA -nt AUTO -bb 1000 -alrt 1000 -m TESTMERGE --runs 5 > ' +  prefix + '_iqtree.err' + '\n')

	pbs.close()


def write_partitions(prefix, length):

	# create and open partition output file
	part = open(prefix + '_iqtree.partition', "w")

	part.write('#nexus\n')
	part.write('begin sets;\n')
	part.write('    ' + 'charset c1 = 1-' + str(length-2) + '\\3' + ';\n')
	part.write('    ' + 'charset c2 = 2-' + str(length-1) + '\\3' + ';\n')
	part.write('    ' + 'charset c3 = 3-' + str(length)   + '\\3' + ';\n')
	part.write('end;\n')

	part.close()

def main():
	prefix = parse_filename()
	length = get_alignment_length()
	write_pbs(prefix)
	write_partitions(prefix, length)


# get the arguments before calling main
args = get_args()


# execute the program by calling main
if __name__ == "__main__":
	main()

