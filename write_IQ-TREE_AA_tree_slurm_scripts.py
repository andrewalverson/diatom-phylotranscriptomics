#! /usr/bin/env python3

import sys
import os
import re
import argparse
import statistics
import random
from Bio import SeqIO


def get_args():
	# create an ArgumentParser object ('parser') that will hold all the information necessary to parse the command line
	parser = argparse.ArgumentParser()

	# add arguments
	parser.add_argument( "fasta", help="FASTA AA alignment" )
	parser.add_argument( "--queue", "-q", help="Pinnacle queue", default="comp72" )

	return parser.parse_args()


def parse_filename():

	# parse name of input file in order to set prefix for output file names
	# a = args.fasta.split("_")
	
	# needed to change this for the ortholog trees
	# a = args.fasta[:-3] + '_alignment'
	a = args.fasta.split("_trimal_final.fa")
	return a[0]

def get_cpus_time():
	num_cpus = 0
	time     = 0

	if args.queue == 'comp06':
		num_cpus = 24
		time = 6
	elif args.queue == 'comp72':
		num_cpus = 24
		time = 72
	elif args.queue == 'c1427':
		num_cpus = 24
		time = 3600
	elif args.queue == 'c1421':
		num_cpus = 24
		time = 3600
	elif args.queue == 'himem06':
		num_cpus = 24
		time = 6
	elif args.queue == 'himem72':
		num_cpus = 24
		time = 72
	elif queue == 'tres72':
		num_cpus = 32
		time = 72

	return(time, num_cpus)


def write_slurm(prefix):

	queue_list = ['comp72','himem72', 'c1421', 'c1427', 'tres72', 'c1427', 'tres72']
 
	queue, time, num_cpus = True, True, True

	queue = random.choice(queue_list)
	time, num_cpus = get_cpus_time(queue)


	# create and open PBS output file
	slurm = open(prefix + '.slurm', "w")

	slurm.write('#!/bin/bash' + '\n')
	slurm.write('#SBATCH --job-name=' + prefix + '-AA' + '\n')

	if args.queue == 'c1427':
		slurm.write('#SBATCH --partition condo' + '\n')
		slurm.write('#SBATCH --constraint aja&0gpu&768gb' + '\n')
		slurm.write('#SBATCH --nodes=1' + '\n')
		slurm.write('#SBATCH --qos condo')
		slurm.write('#SBATCH --tasks-per-node=24' + '\n')
		slurm.write('#SBATCH --time=3600:00:00' + '\n')
	elif args.queue == 'c1421':
		slurm.write('#SBATCH --partition condo' + '\n')
		slurm.write('#SBATCH --constraint aja&0gpu&192gb' + '\n')
		slurm.write('#SBATCH --qos condo')
		slurm.write('#SBATCH --nodes=1' + '\n')
		slurm.write('#SBATCH --tasks-per-node=24' + '\n')
		slurm.write('#SBATCH --time=3600:00:00' + '\n')
	else:
		slurm.write('#SBATCH --partition' + ' ' + args.queue + '\n')
		slurm.write('#SBATCH --nodes=1' + '\n')
		slurm.write('#SBATCH --tasks-per-node=' + str(num_cpu) + '\n')
		slurm.write('#SBATCH --time=' + str(time) + ':00:00' + '\n')
	
	slurm.write('\n')

	slurm.write('module purge' + '\n')
	slurm.write('module load intel/18.0.1 impi/18.0.1 mkl/18.0.1' + '\n')
	slurm.write('module load iqtree/1.6.12' + '\n')
	slurm.write('\n')
	slurm.write('cd $SLURM_SUBMIT_DIR' + '\n')
	slurm.write('\n')

	slurm.write('# run IQ-TREE' + '\n')


	# ortholog trees
	slurm.write('iqtree -s ' + args.fasta + ' -pre ' + prefix + ' -st AA -nt AUTO -bb 1000 -m TEST --runs 10 >> ' +  prefix + '.err' + '\n')

	# homolog trees
	# slurm.write('iqtree -s ' + args.fasta + ' -pre ' + prefix + ' -st AA -nt AUTO -bb 1000 -alrt 1000 -m TEST --runs 5 >> ' +  prefix + '.err' + '\n')

	slurm.close()


def main():
	prefix = parse_filename()
	write_slurm(prefix)


# get the arguments before calling main
args = get_args()


# execute the program by calling main
if __name__ == "__main__":
	main()

