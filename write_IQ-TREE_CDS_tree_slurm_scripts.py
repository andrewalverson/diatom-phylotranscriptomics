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
	parser.add_argument( "--queue", "-q", help="Pinnacle queue", default="comp72" )
	parser.add_argument( "--nt", "-n", help="Number of threads", default="AUTO" )
	parser.add_argument( "--cpu", "-c", help="Slurm CPU request", default="4" )

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


def get_cpus_time():
	num_cpus = 0
	time     = 0

	if args.queue == 'comp06':
		# num_cpus = 32
		time = 6
	elif args.queue == 'comp72':
		# num_cpus = 32
		time = 72
	elif args.queue == 'c1427':
		num_cpus = 24
		time = 3600
	elif args.queue == 'c1421':
		num_cpus = 24
		time = 3600
	elif args.queue == 'himem06':
		# num_cpus = 24
		time = 6
	elif args.queue == 'himem72':
		# num_cpus = 24
		time = 72

	return(time, num_cpus)


def write_slurm(prefix):

	time, num_cpu = get_cpus_time()

	# create and open slurm output file
	slurm = open(prefix + '_iqtree.slurm', "w")

	slurm.write('#!/bin/bash' + '\n')
	slurm.write('#SBATCH --job-name=' + prefix + '.iqtree-cds.$PBS_JOBID' + '\n')

	if args.queue == 'c1427':
		slurm.write('#SBATCH --partition condo' + '\n')
		slurm.write('#SBATCH --constraint aja&0gpu&768gb' + '\n')
		slurm.write('#SBATCH --qos condo')
		slurm.write('#SBATCH --tasks-per-node=24' + '\n')
		slurm.write('#SBATCH --time=3600:00:00' + '\n')
	elif args.queue == 'c1421':
		slurm.write('#SBATCH --partition condo' + '\n')
		slurm.write('#SBATCH --constraint aja&0gpu&192gb' + '\n')
		slurm.write('#SBATCH --qos condo')
		slurm.write('#SBATCH --tasks-per-node=24' + '\n')
		slurm.write('#SBATCH --time=3600:00:00' + '\n')
	else:
		slurm.write('#SBATCH --partition' + ' ' + args.queue + '\n')
		slurm.write('#SBATCH --nodes=1' + '\n')
		slurm.write('#SBATCH --tasks-per-node=' + str(args.cpu) + '\n')
		slurm.write('#SBATCH --time=' + str(time) + ':00:00' + '\n')
	
	slurm.write('\n')
	slurm.write('module purge' + '\n')
	slurm.write('module load intel/18.0.1 impi/18.0.1 mkl/18.0.1' + '\n')
	slurm.write('module load iqtree/1.6.12' + '\n')
	slurm.write('\n')

	slurm.write('cd $SLURM_SUBMIT_DIR' + '\n')
	slurm.write('\n')

	slurm.write('# run IQ-TREE' + '\n')
	slurm.write('iqtree -s ' + args.fasta + ' -pre ' + prefix + '_iqtree' + ' -spp ' + prefix +  '_iqtree.partition' + ' -st DNA -nt ' + args.nt + ' -bb 1000 -alrt 1000 -m TESTMERGE --runs 5 >> ' +  prefix + '_iqtree.err' + '\n')

	slurm.close()


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
	write_slurm(prefix)
	write_partitions(prefix, length)


# get the arguments before calling main
args = get_args()


# execute the program by calling main
if __name__ == "__main__":
	main()

