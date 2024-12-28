#! /usr/bin/env python3

import sys
import os
import re
import argparse
import random

def get_args():
	# create an ArgumentParser object ('parser') that will hold all the information necessary to parse the command line
	parser = argparse.ArgumentParser()

	# add arguments
	parser.add_argument( "fasta", help='FASTA, output from transrate-filtering step: *_nuclear.good_transcripts.fa')

	return parser.parse_args()


def get_cpu_time(queue):
	num_cpu = 0
	time     = 0

	if queue == 'comp72':
		num_cpu = 8
		time = 72
	elif queue == 'c1421':
		num_cpu = 12
		time = 3600
	elif queue == 'c1427':
		num_cpu = 12
		time = 3600
	elif queue == 'himem06':
		num_cpu = 12
		time = 6
	elif queue == 'himem72':
		num_cpu = 12
		time = 72

	return(time, num_cpu)


def print_slurm():

	# queue = random.choice(['c1307', 'c1308'])
	queue = 'himem06'
	time, num_cpu = get_cpu_time(queue)

	job   = args.fasta.replace('_nuclear.good_transcripts.short_name.filtered_transcripts.fa', '.corset')
	left  = args.fasta.replace('_nuclear.good_transcripts.short_name.filtered_transcripts.fa', '_trimmed_filtered_1.fq')
	right = args.fasta.replace('_nuclear.good_transcripts.short_name.filtered_transcripts.fa', '_trimmed_filtered_2.fq')

	salmon_index = args.fasta.replace('_nuclear.good_transcripts.short_name.filtered_transcripts.fa', '_salmon_index')
	salmon_quant  = args.fasta.replace('_nuclear.good_transcripts.short_name.filtered_transcripts.fa', '_salmon_quant')
	eq_classes = salmon_quant + '/aux_info/eq_classes.txt'
	corset_out = args.fasta.replace('_nuclear.good_transcripts.short_name.filtered_transcripts.fa', '_salmon_quant-clusters.txt')

	print('#!/bin/bash' + '\n')
	print('#SBATCH --job-name=' + job)

	if queue == 'c1427':
		print('#SBATCH --partition condo')
		print('#SBATCH --constraint aja&0gpu&768gb')
		print('#SBATCH --qos condo')
		print('#SBATCH --nodes=1')
		print('#SBATCH --tasks-per-node=' + str(num_cpu))
		print('#SBATCH --time=3600:00:00')
	elif queue == 'c1421':
		print('#SBATCH --partition condo')
		print('#SBATCH --constraint aja&0gpu&192gb')
		print('#SBATCH --qos condo')
		print('#SBATCH --nodes=1')
		print('#SBATCH --tasks-per-node=' + str(num_cpu))
		print('#SBATCH --time=3600:00:00')
	else:
		print('#SBATCH --partition' + ' ' + queue)
		print('#SBATCH --nodes=1')
		print('#SBATCH --tasks-per-node=' + str(num_cpu))
		print('#SBATCH --time=' + str(time) + ':00:00')
	

	print('#SBATCH -o', job + '.%j.out')
	print('#SBATCH -e', job + '.%j.err')


	print()
	print('export PATH=$PATH:/home/aja/local/bin')
	print('export PATH=$PATH:/home/aja/local/src/yang-2021/phylogenomic_dataset_construction/scripts')

	print()
	print('export OMP_NUM_THREADS=' + str(num_cpu))

	print()
	print('module load python/2.7.13-anaconda')
	# print('module load salmon/0.14.1')
	print('module load salmon/1.4.0')

	print()
	print('cd $SLURM_SUBMIT_DIR')

	# print()
	# print('# build salmon index')
	# print('salmon index --index', salmon_index, '--transcripts', args.fasta)

	# print()
	# print('# quantify transcripts with salmon')
	# print('salmon quant --index', salmon_index, '--libType A', '--dumpEq', '-1', left, '-2', right, '--output', salmon_quant)

	# print()
	# print('# run Corset')
	# print('gunzip', eq_classes + '.gz')
	# print('corset -i salmon_eq_classes', eq_classes, '-l 5 -m 5 -x 25 -f true', '-p', salmon_quant)

	print()
	print('# run filter_corset_output.py')
	print('python2 /home/aja/local/src/yang-2021/phylogenomic_dataset_construction/scripts/filter_corset_output.py', \
	 args.fasta, corset_out, '.')

def main():
	print_slurm()


# get the arguments before calling main
args = get_args()


# execute the program by calling main
if __name__ == "__main__":
	main()
