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
	parser.add_argument( "--queue", "-q", help="Pinnacle queue (comp06, comp72, c1307 [768 GB], c1308 [192 GB], himem06, himem72)", default="comp72" )
	parser.add_argument( "--jobID", "-j", help='job name', default='aja')
	parser.add_argument( "fasta", help='trinity assembly FASTA file name')

	return parser.parse_args()


def get_cpus_time():
	num_cpus = 0
	time     = 0

	if args.queue == 'comp72':
		num_cpus = 32
		time = 72
	elif args.queue == 'c1421':
		num_cpus = 24
		time = 3600
	elif args.queue == 'c1427':
		num_cpus = 24
		time = 3600
	elif args.queue == 'himem06':
		num_cpus = 24
		time = 6
	elif args.queue == 'himem72':
		num_cpus = 24
		time = 72

	return(time, num_cpus)


def parse_names():
	input_directory = args.fasta + '.transdecoder_dir'
	output_file = args.fasta + '.pfam.domtblout'

	return(input_directory, output_file)

def print_slurm():

	time, num_cpus = get_cpus_time()
	input_directory, output_file = parse_names()

	# alverson_lab and NIES-3581
	pfam   = args.fasta.replace('largest_cluster_transcripts.fa', 'fa.hmmer.pfam.domtblout')
	uniref = args.fasta.replace('largest_cluster_transcripts.fa', 'fa.vs.uniref.blastp')


	print('#!/bin/bash' + '\n')
	print('#SBATCH --job-name=' + args.fasta + '.hmmer.pfam')
	print('#SBATCH --mail-type=ALL')
	print('#SBATCH --mail-user=aja@uark.edu')

	if args.queue == 'c1427':
		print('#SBATCH --partition condo')
		print('#SBATCH --constraint aja&0gpu&768gb')
		print('#SBATCH --qos condo')
		print('#SBATCH --nodes=1')
		print('#SBATCH --tasks-per-node=' + str(num_cpus))
		print('#SBATCH --time=3600:00:00')
	elif args.queue == 'c1421':
		print('#SBATCH --partition condo')
		print('#SBATCH --constraint aja&0gpu&192gb')
		print('#SBATCH --qos condo')
		print('#SBATCH --nodes=1')
		print('#SBATCH --tasks-per-node=' + str(num_cpus))
		print('#SBATCH --time=3600:00:00')
	else:
		print('#SBATCH --partition' + ' ' + args.queue)
		print('#SBATCH --nodes=1')
		print('#SBATCH --tasks-per-node=' + str(num_cpus))
		print('#SBATCH --time=' + str(time) + ':00:00')
	
	print()

	print('module load intel/18.0.1 impi/18.0.1 mkl/18.0.1')
	print('module load perl/5.10.1')
	print('module load blast/2.9.0+')
	print('module load hmmer/3.1b2')
	print('module load transdecoder/2.0.1')
	print('module load cdhit/4.6.8')

	print()
	print('cd $SLURM_SUBMIT_DIR')
	print()


	print('# TransDecoder step 2.2: HMMR search search long ORFs to Pfam database')
	print('/share/apps/bioinformatics/hmmer/3.1b2/bin/hmmscan --cpu', str(num_cpus), '--domtblout', pfam + '.out',  \
		'/storage/aja/pfam_db/Pfam-A.hmm', input_directory + '/' + 'longest_orfs.pep')

	# print('# TransDecoder step 3: Use Transdecoder to translate nuclear contigs using homology search information')
	# print('/share/apps/bioinformatics/transdecoder/5.5.0/TransDecoder.Predict -t', args.input_fasta,  '--retain_pfam_hits', output_file,  \
	#	'--retain_blastp_hits uniref.blastp --single_best_only')

	# print('# Run CD-HIT to remove redundant/overlapping contigs')
	# print('/share/apps/bioinformatics/cdhit/v4.6.8-2017-1208/cd-hit -i L503_nuclear.fa.transdecoder.pep -o L503_cd-hit.pep \
	#	-c 0.99 -n 5 -d 1000 -M 64000 -T 0')




def main():
	print_slurm()

# get the arguments before calling main
args = get_args()


# execute the program by calling main
if __name__ == "__main__":
	main()

