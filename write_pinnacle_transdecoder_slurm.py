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
	parser.add_argument( "fasta", help='FASTA, output from transrate-filtering step: *_nuclear.good_transcripts.fa')
	parser.add_argument( "--queue", "-q", help="Pinnacle queue (comp06, comp72, c1307 [768 GB], c1308 [192 GB], himem06, himem72)", default="comp72" )

	return parser.parse_args()


def get_cpus_time():
	num_cpus = 0
	time     = 0

	if args.queue == 'comp72':
		num_cpus = 16
		time = 72
	elif args.queue == 'c1421':
		num_cpus = 12
		time = 3600
	elif args.queue == 'c1427':
		num_cpus = 12
		time = 3600
	elif args.queue == 'himem06':
		num_cpus = 24
		time = 6
	elif args.queue == 'himem72':
		num_cpus = 12
		time = 72

	return(time, num_cpus)


def print_slurm():

	time, num_cpu = get_cpus_time()

	job = args.fasta.replace('_nuclear.fa', '.transdecoder')
	left  = args.fasta.replace('_nuclear.good_transcripts.fa', '_trimmed_filtered_1.fq.gz')
	right = args.fasta.replace('_nuclear.good_transcripts.fa', '_trimmed_filtered_2.fq.gz')

	print('#!/bin/bash' + '\n')
	print('#SBATCH --job-name=' + job)
	print('#SBATCH --mail-type=ALL')
	print('#SBATCH --mail-user=aja@uark.edu')

	if args.queue == 'c1427':
		print('#SBATCH --partition condo')
		print('#SBATCH --constraint aja&0gpu&768gb')
		print('#SBATCH --qos condo')
		print('#SBATCH --nodes=1')
		print('#SBATCH --tasks-per-node=24')
		print('#SBATCH --time=3600:00:00')
	elif args.queue == 'c1421':
		print('#SBATCH --partition condo')
		print('#SBATCH --constraint aja&0gpu&192gb')
		print('#SBATCH --qos condo')
		print('#SBATCH --nodes=1')
		print('#SBATCH --tasks-per-node=1')
		print('#SBATCH --time=3600:00:00')
	else:
		print('#SBATCH --partition' + ' ' + args.queue)
		print('#SBATCH --qos comp')
		print('#SBATCH --nodes=1')
		print('#SBATCH --tasks-per-node=' + str(num_cpu))
		print('#SBATCH --time=' + str(time) + ':00:00')
	
	print()

	print('export OMP_NUM_THREADS=' + str(num_cpu))
	print('module load python/2.7.13-anaconda')
	print('module load salmon/0.14.1')

	print()
	print('cd $SLURM_SUBMIT_DIR')
	print()

	print('python2 /home/aja/local/src/yang-2021/phylogenomic_dataset_construction/scripts/corset_wrapper.py', \
	 args.fasta, left, right, str(num_cpu), '.', 'salmon')


def main():
	print_slurm()


# get the arguments before calling main
args = get_args()


# execute the program by calling main
if __name__ == "__main__":
	main()
