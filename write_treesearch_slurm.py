#! /usr/bin/env python3

import argparse
import random
import subprocess

def get_args():
	# create an ArgumentParser object ('parser') that will hold all the information necessary to parse the command line
	parser = argparse.ArgumentParser()

	# add arguments
	parser.add_argument( "fasta", help="Name of input FASTA file" )

	return parser.parse_args()


def fasta_count():
	command = "grep " + "'>' " + args.fasta + ' | wc -l'
	count = subprocess.check_output(command, shell=True, universal_newlines=True)

	return count


def get_cpus_time(queue):

	if queue == 'comp06':
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


def write_slurm(seq_count):

	job_type = ''

	# parse name of input file in order to set prefix for output file names
	a = args.fasta.split(".")
	file_prefix = a[0]

	if int(seq_count) < 1000:
		job_type = 'iqtree'
	else:
		job_type = 'fasttree'

	'''
	queue_list = ['comp06', 'comp06', 'comp06', 'comp72', 'comp72', 'comp72', 'himem06', 'himem06', 'himem06', 'himem72', 'himem72', 'himem72', 
	'tres06', 'tres06', 'tres06', 'tres72', 'tres72', 'tres72', 'c1421', 'c1427', 'c1421', 'c1427']

	queue_list = ['comp06', 'tres06', 'comp72', 'tres72', 'comp06', 'tres06', 'comp72', 'tres72', 'himem72', 'himem06', 'c1427']
	queue_list = ['c1421', 'c1427', 'tres72']
	'''

	queue_list = ['c1427']


	queue = random.choice(queue_list)
	time, num_cpus = get_cpus_time(queue)

	job = file_prefix + '_' + job_type
	outfile = job + '.slurm'

	out = open(outfile, 'w')

	out.write('#!/bin/bash' + '\n')
	out.write('#SBATCH --job-name=' + job + '\n')

	if queue == 'c1421':
		out.write('#SBATCH --partition condo' + '\n')
		out.write('#SBATCH --constraint aja&0gpu&192gb' + '\n')
		out.write('#SBATCH --qos condo' + '\n')
		out.write('#SBATCH --nodes=1' + '\n')
		out.write('#SBATCH --tasks-per-node=' + str(num_cpus) + '\n')
		out.write('#SBATCH --time=' + str(time) + ':00:00' + '\n')
	elif queue == 'c1427':
		out.write('#SBATCH --partition condo' + '\n')
		out.write('#SBATCH --constraint aja&0gpu&768gb' + '\n')
		out.write('#SBATCH --qos condo' + '\n')
		out.write('#SBATCH --nodes=1' + '\n')
		out.write('#SBATCH --tasks-per-node=' + str(num_cpus) + '\n')
		out.write('#SBATCH --time=' + str(time) + ':00:00' + '\n')
	else:
		out.write('#SBATCH --partition' + ' ' + queue + '\n')
		out.write('#SBATCH --nodes=1' + '\n')
		out.write('#SBATCH --tasks-per-node=' + str(num_cpus) + '\n')
		out.write('#SBATCH --time=' + str(time) + ':00:00' + '\n')
	out.write('\n')
	
	out.write('export OMP_NUM_THREADS=' + str(num_cpus) + '\n')
	out.write('\n')

	out.write('module load intel/18.0.1 impi/18.0.1 mkl/18.0.1' + '\n')
	out.write('module load iqtree/1.6.12' + '\n')
	out.write('module load fasttree/2.1.10' + '\n')
	out.write('\n')

	out.write('cd $SLURM_SUBMIT_DIR' + '\n')
	out.write('\n')

	if job_type == 'fasttree':
		out.write('# run FASTTREE' + '\n')
		out.write('# -pseudo for gappy alignments, -gamma for good brlens, -nosupport for topology only' + '\n')
		out.write('FastTreeMP -wag -gamma -pseudo -nosupport ' + args.fasta + ' > ' + file_prefix + '.' + job_type + '\n')
	else:
		out.write('# run IQTREE' + '\n')
		# version with '-safe'
		out.write('iqtree -s ' + args.fasta + ' -pre ' + file_prefix + ' -nt ' + str(num_cpus) + ' -safe -st AA -nt AUTO -m TEST --runs 1 > ' + job + '.err')
		# out.write('iqtree -s ' + args.fasta + ' -pre ' + file_prefix + ' -nt ' + str(num_cpus) + ' -st AA -nt AUTO -m TEST --runs 1 > ' + job + '.err')

	out.write('\n')


def main():
	seq_count = fasta_count()
	write_slurm(seq_count)

# get the arguments before calling main
args = get_args()


# execute the program by calling main
if __name__ == "__main__":
	main()

