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
		num_cpus = 16
		time = 6
	elif queue == 'comp72':
		num_cpus = 16
		time = 72
	elif queue == 'c1421':
		num_cpus = 24
		time = 3600
	elif queue == 'c1427':
		num_cpus = 24
		time = 3600
	elif queue == 'himem06':
		num_cpus = 12
		time = 6
	elif queue == 'himem72':
		num_cpus = 12
		time = 72
	elif queue == 'tres06':
		num_cpus = 16
		time = 6
	elif queue == 'tres72':
		num_cpus = 16
		time = 72

	return(time, num_cpus)


def write_slurm(seq_count):

	a = args.fasta.split(".")
	file_prefix = a[0]

	queue_lt100 = ['tres06', 'comp06']
	queue_gt100 = ['comp72', 'comp72', 'himem72', 'tres72']

	'''
	queue_lt100 = ['tres06' 'tres72', 'c1421', 'c1427']
	queue_gt100 = ['comp72', 'himem72', 'tres72', 'tres72', 'c1421', 'c1421', 'c1427', 'c1427']
	'''

	if int(seq_count) > 100:
		queue = random.choice(queue_gt100)
	else:
		queue = random.choice(queue_lt100)
	
	time, num_cpus = get_cpus_time(queue)

	out_prefix = file_prefix + '_iqtree'
	outfile    = out_prefix + '.slurm'

	out = open(outfile, 'w')

	out.write('#!/bin/bash' + '\n')
	out.write('#SBATCH --job-name=' + out_prefix + '\n')

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
	elif queue == 'tres06' or queue == 'tres72' or queue == 'razr06' or queue == 'razr72':
		out.write('#SBATCH --partition' + ' ' + queue + '\n')
		out.write('#SBATCH --qos tres' + '\n')
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
	out.write('\n')

	out.write('cd $SLURM_SUBMIT_DIR' + '\n')
	out.write('\n')

	out.write('# run IQTREE' + '\n')
	out.write('iqtree -s ' + args.fasta + ' -pre ' + file_prefix + ' -nt ' + str(num_cpus) + ' -st AA -mset WAG,LG --runs 5 -bb 1000 -wbt > ' + out_prefix + '.err')
	out.write('\n')


def main():
	seq_count = fasta_count()
	write_slurm(seq_count)

# get the arguments before calling main
args = get_args()


# execute the program by calling main
if __name__ == "__main__":
	main()

