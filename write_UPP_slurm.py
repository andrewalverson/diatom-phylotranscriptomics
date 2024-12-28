#! /usr/bin/env python3

import random
import re
import argparse
import statistics
import numpy
from Bio import SeqIO


def get_args():

	# create an ArgumentParser object ('parser') that will hold all the information necessary to parse the command line
	parser = argparse.ArgumentParser(description = "write a PBS job script to run UPP with the '-M' setting")

	# add arguments
	parser.add_argument( "fasta", help="unaligned FASTA file" )
	parser.add_argument( "-m", "--median_full_length", action='store_true', default=False, help="use UPP's '-M' argument (default=False)")


	return(parser.parse_args())


def get_M_length():

	# compute median sequence length for this dataset
	lengths = [] # list of sequence lengths

	# open and parse FASTA file, storing sequence lengths
	fasta_sequences = SeqIO.parse(open(args.fasta),'fasta')
	for record in fasta_sequences:
	    lengths.append(len(record.seq))
	fasta_sequences.close()

	return(int(numpy.quantile(lengths, 0.75)), len(lengths))


def get_cpus_time(queue):

	if queue == 'comp01':
		num_cpus = 4
		time = 1
	elif queue == 'comp06':
		num_cpus = 4
		time = 6
	elif queue == 'comp72':
		num_cpus = 4
		time = 72
	elif queue == 'c1421':
		num_cpus = 4
		time = 3600
	elif queue == 'c1427':
		num_cpus = 4
		time = 3600
	elif queue == 'himem06':
		num_cpus = 4
		time = 6
	elif queue == 'himem72':
		num_cpus = 4
		time = 72
	elif queue == 'tres06':
		num_cpus = 4
		time = 6
	elif queue == 'tres72':
		num_cpus = 4
		time = 72

	return(time, num_cpus)


def print_slurm():

	# set input directory
	# indir = '/storage/aja/orthofinder/orthogroup_fasta_occupancy_20_no_stops'
	indir = '/scratch/$SLURM_JOB_ID/tmp'

	# set output directory
	# outdir = '/storage/aja/orthofinder/gene_trees_round_1/upp'
	outdir = '$SLURM_SUBMIT_DIR'

	# get median sequence length and number of sequences in the alignment
	M_length, num_seqs = get_M_length()

	# queue_gt300 = ['comp72', 'comp72', 'comp72', 'himem72', 'himem72', 'himem72', 'tres72', 'tres72', 'tres72']
	# queue_lt300 = ['comp06', 'comp06', 'comp72', 'comp72', 'comp72', 'himem06', 'himem06', 'himem72', 'himem72', 'himem72', 'tres72', 'tres72', 'tres72']
	queue_gt300 = ['comp01', 'tres06']
	queue_lt300 = ['comp01', 'tres06']
	# queue_list = ['himem72', 'c1421', 'c1427']

	queue, time, num_cpus = True, True, True

	if(int(num_seqs) < 300):
		queue = random.choice(queue_lt300)
		time, num_cpus = get_cpus_time(queue)
	else:
		queue = random.choice(queue_gt300)
		time, num_cpus = get_cpus_time(queue)

	# parse name of input file in order to set prefix for output file names
	# had to change this for ortholog files from Yang's RT output
	# a = args.fasta.split(".")
	# prefix = a[0]

	prefix = args.fasta[:-3]

	job = prefix + '.upp'

	print('#!/bin/bash' + '\n')
	print('#SBATCH --job-name=' + job)

	if queue == 'c1421':
		time = '3600'
		print('#SBATCH --partition condo')
		print('#SBATCH --constraint aja&0gpu&192gb')
		print('#SBATCH --qos condo')
		print('#SBATCH --nodes=1')
		print('#SBATCH --tasks-per-node=' + str(num_cpus))
		print('#SBATCH --time=' + str(time) + ':00:00')
	elif queue == 'c1427':
		time = '3600'
		print('#SBATCH --partition condo')
		print('#SBATCH --constraint aja&0gpu&768gb')
		print('#SBATCH --qos condo')
		print('#SBATCH --nodes=1')
		print('#SBATCH --tasks-per-node=' + str(num_cpus))
		print('#SBATCH --time=' + str(time) + ':00:00')
	else:
		print('#SBATCH --partition', queue)
		print('#SBATCH --nodes=1')
		print('#SBATCH --tasks-per-node=' + str(num_cpus))
		print('#SBATCH --time=' + str(time) + ':00:00')
	
	print()

	print('export OMP_NUM_THREADS=' + str(num_cpus))
	print()

	print('module purge')
	print('module load gcc/8.3.1  python/3.7.3-anaconda')
	print('source /share/apps/bin/conda-3.7.3.sh')
	print('conda activate bioconda3-el7')

	print()

	# need to put blast module after conda activate to override internal version
	# PASTA 1.8.6 installed in conda but module needed for run_pasta
	print('module load blast/2.11.0+bin PASTA/1.8.6')
	print()
	
	# note: mafft in bioconda3 are newer than module versions

	# create and make a temp directory in /scratch
	print('export TMPDIR=/scratch/$SLURM_JOB_ID/tmp')
	print('mkdir $TMPDIR')
	print('cp', args.fasta, '/scratch/$SLURM_JOB_ID/')
	print()

	print('# move into the temp directory')
	print('cd $TMPDIR')
	print()

	print('# run UPP')
	

	# with '-M'
	if args.median_full_length:

		print('run_upp.py', '-m amino', '-x', str(num_cpus), '-s', '/scratch/$SLURM_JOB_ID/' + args.fasta, '-d', '$TMPDIR', '-o', prefix, '-M', M_length)

	# with defaults
	else:
		print('run_upp.py', '-m amino', '-x', str(num_cpus), '-s', '/scratch/$SLURM_JOB_ID/' + args.fasta, '-d', '$TMPDIR', '-o', prefix)

	print()
	print('# copy output files')
	print('cp $TMPDIR/*_alignment* $SLURM_SUBMIT_DIR')
	print('cp $TMPDIR/*_insertion_columns.txt $SLURM_SUBMIT_DIR')


def main():
	print_slurm()


# get the arguments before calling main
args = get_args()


# execute the program by calling main
if __name__ == "__main__":
	main()
