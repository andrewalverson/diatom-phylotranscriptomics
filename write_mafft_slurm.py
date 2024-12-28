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


def get_cpus_time():

	#queue_list = ['comp06', 'comp06', 'comp06', 'himem06', 'himem06', 'himem72', 'himem72', 'tres06', 'tres06', 'tres72', 'tres72', 'tres72']
	queue_list = ['tres72']
	queue = random.choice(queue_list)

	num_cpus = 0
	time     = 0

	if queue == 'comp06':
		num_cpus = 4
		time = 6
	elif queue == 'tres06':
		num_cpus = 4
		time = 6
	elif queue == 'comp72':
		num_cpus = 4
		time = 18
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
	elif queue == 'tres72':
		num_cpus = 4
		time = 72

	return(queue, time, num_cpus)


def fasta_count():
	command = "grep " + "'>' " + args.fasta + ' | wc -l'
	count = subprocess.check_output(command, shell=True, universal_newlines=True)
	return count


def print_slurm(seq_count, queue, time, num_cpus):

	# parse name of input file in order to set prefix for output file names
	a = args.fasta.split(".")
	prefix = a[0]

	job = prefix + '-mafft'
	
	print('#!/bin/bash' + '\n')
	print('#SBATCH --job-name=' + job)

	if queue == 'c1427':
		print('#SBATCH --partition condo')
		print('#SBATCH --constraint aja&0gpu&768gb')
		print('#SBATCH --qos condo')
		print('#SBATCH --nodes=1')
		print('#SBATCH --tasks-per-node=' + str(num_cpus))
		print('#SBATCH --time=' + str(time) + ':00:00')
	elif queue == 'c1421':
		print('#SBATCH --partition condo')
		print('#SBATCH --constraint aja&0gpu&192gb')
		print('#SBATCH --qos condo')
		print('#SBATCH --nodes=1')
		print('#SBATCH --tasks-per-node=' + str(num_cpus))
		print('#SBATCH --time=' + str(time) + ':00:00')
	else:
		print('#SBATCH --partition' + ' ' + queue)
		print('#SBATCH --nodes=1')
		print('#SBATCH --tasks-per-node=' + str(num_cpus))
		print('#SBATCH --time=' + str(time) + ':00:00')
	
	print()
	print('export OMP_NUM_THREADS=' + str(num_cpus))
	print()

	print('module load intel/18.0.1 impi/18.0.1 mkl/18.0.1')

	print('module load python/3.7.3-anaconda')
	print()
	print('conda init bash')
	print('source activate ~/.conda/envs/mafft')
	print()

	print('cd $SLURM_SUBMIT_DIR')
	print()

	print('# run MAFFT')
	if(int(seq_count) > 200):
		print('mafft --anysymbol --maxiterate 1000 --thread ' + str(num_cpus) + ' ' + args.fasta + ' > ' + prefix + '.mafft')
	else:
		print('mafft --anysymbol --genafpair --maxiterate 1000 --thread ' + str(num_cpus) + ' ' + args.fasta + ' > ' + prefix + '.mafft')


def main():
	queue, time, num_cpus = get_cpus_time()
	seq_count = fasta_count()
	print_slurm(seq_count, queue, time, num_cpus)


# get the arguments before calling main
args = get_args()


# execute the program by calling main
if __name__ == "__main__":
	main()

