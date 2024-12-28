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

	# queue_list = ['aja', 'onenode16core', 'mem768GB32core', 'mem768GB32core', 'mem768GB32core', 'med12core', 'tiny12core', 'tiny16core']
	queue_list = ['aja', 'onenode16core', 'mem512GB64core', 'mem768GB32core', 'med12core', 'tiny12core', 'tiny16core', 'med16core']
	queue = random.choice(queue_list)

	num_cpus = 6
	time     = 6

	if queue == 'aja':
		num_cpus = 4
		time = 72
	elif queue == 'onenode16core':
		num_cpus = 4
		time = 72
	elif queue == 'mem768GB32core':
		num_cpus = 4
		time = 72
	elif queue == 'mem512GB64core':
		num_cpus = 4
		time = 72
	elif queue == 'med12core':
		num_cpus = 4
		time = 72
	elif queue == 'serial12core':
		num_cpus = 4
		time = 72
	elif queue == 'mem96GB12core':
		num_cpus = 4
		time = 72
	elif queue == 'tiny12core':
		num_cpus = 4
		time = 6
	elif queue == 'tiny16core':
		num_cpus = 4
		time = 6
	elif queue == 'student12core':
		num_cpus = 4
		time = 6

	return(queue, time, num_cpus)


def fasta_count():
	command = "grep " + "'>' " + args.fasta + ' | wc -l'
	count = subprocess.check_output(command, shell=True, universal_newlines=True)
	return count


def print_pbs(seq_count, queue, time, num_cpus):

	# parse name of input file in order to set prefix for output file names
	a = args.fasta.split(".")
	prefix = a[0]

	job = prefix + '-mafft'
	
	
	# print top of PBS script
	print('#PBS -N', prefix+'.mafft')
	print('#PBS -q', queue)
	print('#PBS -j oe')
	print('#PBS -o',  prefix+'.mafft.$PBS_JOBID')
	print('#PBS -l nodes=1:ppn=' + str(num_cpus))
	print('#PBS -l walltime=' + str(time) + ':00:00')
	print()

	print('export OMP_NUM_THREADS=' + str(num_cpus))
	print()

	print('# load modules')
	print('module load gcc/7.2.1 java/sunjdk_1.8.0 python/3.7.3-anaconda-razor')
	print('module load mafft')
	print()

	# print('conda init bash')
	# print('source activate ~/.conda/envs/mafft-razor')
	# print()

	print('cd $PBS_O_WORKDIR')
	print()

	print('# run MAFFT')
	if(int(seq_count) > 200):
		print('mafft --anysymbol --maxiterate 1000 --thread ' + str(num_cpus) + ' ' + args.fasta + ' > ' + prefix + '.mafft')
	else:
		print('mafft --anysymbol --genafpair --maxiterate 1000 --thread ' + str(num_cpus) + ' ' + args.fasta + ' > ' + prefix + '.mafft')


def main():
	queue, time, num_cpus = get_cpus_time()
	seq_count = fasta_count()
	print_pbs(seq_count, queue, time, num_cpus)


# get the arguments before calling main
args = get_args()


# execute the program by calling main
if __name__ == "__main__":
	main()

