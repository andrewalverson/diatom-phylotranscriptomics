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

	if queue == 'aja':
		num_cpus = 16
		time = 2400
	elif queue == 'tiny16core':
		num_cpus = 16
		time = 6
	elif queue == 'med16core':
		num_cpus = 16
		time = 72
	elif queue == 'med12core':
		num_cpus = 12
		time = 72
	elif queue == 'onenode16core':
		num_cpus = 16
		time = 168
	elif queue == 'mem768GB32core':
		num_cpus = 16
		time = 72
	elif queue == 'mem512GB64core':
		num_cpus = 32
		time = 72

	return(time, num_cpus)


def write_pbs(seq_count):

	job_type = ''

	# parse name of input file in order to set prefix for output file names
	a = args.fasta.split(".")
	file_prefix = a[0]

	if int(seq_count) < 1000:
		job_type = 'iqtree'
	else:
		job_type = 'fasttree'

	queue_list = ['aja', 'onenode16core', 'mem768GB32core', 'mem768GB32core', 'mem512GB64core', 'mem512GB64core', 'mem512GB64core', 'tiny16core', 'tiny16core', 'med16core', 'med16core', 'med12core']
	queue = random.choice(queue_list)
	time, num_cpus = get_cpus_time(queue)

	job = file_prefix + '_' + job_type
	outfile = job + '.pbs'

	out = open(outfile, 'w')

	# print top of PBS script
	out.write('#PBS -N ' + job + '\n')
	out.write('#PBS -q ' + queue + '\n')
	out.write('#PBS -j oe' + '\n')
	out.write('#PBS -o ' +  job + '.$PBS_JOBID' + '\n')
	out.write('#PBS -l nodes=1:ppn=' + str(num_cpus) + '\n')
	out.write('#PBS -l walltime=' + str(time) + ':00:00' + '\n')
	out.write('\n')

	out.write('export OMP_NUM_THREADS=' + str(num_cpus) + '\n')
	out.write('\n')

	out.write('# load modules' + '\n')
	out.write('# module purge' + '\n')
	out.write('module load gcc/7.2.1 java/sunjdk_1.8.0 python/3.7.3-anaconda-razor' + '\n')
	out.write('module load fasttree/2.1.10\n')
	out.write('module load iqtree/1.6.12\n')
	out.write('\n')


	out.write('cd $PBS_O_WORKDIR' + '\n')
	out.write('\n')

	if job_type == 'fasttree':
		out.write('# run FASTTREE' + '\n')
		out.write('# -pseudo for gappy alignments, -gamma for good brlens, -nosupport for topology only' + '\n')
		out.write('FastTreeMP -wag -gamma -pseudo -nosupport ' + args.fasta + ' > ' + file_prefix + '.' + job_type + '\n')
	else:
		out.write('# run IQTREE' + '\n')
		out.write('iqtree -s ' + args.fasta + ' -pre ' + file_prefix + ' -nt ' + str(num_cpus) + ' -st AA -nt AUTO -m TEST --runs 1 > ' + job + '.err')

	out.write('\n')


def main():
	seq_count = fasta_count()
	write_pbs(seq_count)


# get the arguments before calling main
args = get_args()


# execute the program by calling main
if __name__ == "__main__":
	main()

