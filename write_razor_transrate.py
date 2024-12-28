#! /usr/bin/env python3

# load required modules
import argparse

def get_args():
	# create an ArgumentParser object ('parser') that will hold all the information necessary to parse the command line
	parser = argparse.ArgumentParser(description='this script writes a PBS script to run transrate')

	# add positional arguments
	parser.add_argument( "fasta", help='trinity assembly FASTA file name (transdecoder --assembly')

	# add optional arguments
	parser.add_argument( "--queue", "-q", help='razor queue', default='aja')
	parser.add_argument( "--wall",  "-w", help='walltime', type=int, default=72)
	parser.add_argument( "--cpus",  "-c", help='number of cpus', type=int, default='16')

	return parser.parse_args()


def parse_names():
	output = args.fasta + '.transrate'

	return(input_directory, output_file)

def print_pbs():

	job = args.fasta.replace('_nuclear.fa', '.transrate')
	left  = args.fasta.replace('_nuclear.fa', '_trimmed_filtered_1.fq.gz')
	right = args.fasta.replace('_nuclear.fa', '_trimmed_filtered_2.fq.gz')
	outdir = args.fasta + '.transrate'
	base = args.fasta.split('.fa')
	basename = base[0]

	print('#PBS -N', job)
	print('#PBS -q', args.queue)
	print('#PBS -o',  job + '.$PBS_JOBID.out')
	print('#PBS -l nodes=1:ppn=' + str(args.cpus))
	print('#PBS -l walltime=' + str(args.wall) + ':00:00')
	print()

	print('export PATH=$PATH:/home/aja/local/bin')
	print('export PATH=$PATH:/home/aja/local/src/yang-2021/phylogenomic_dataset_construction/scripts')

	print('export OMP_NUM_THREADS=' + str(args.cpus))
	print()

	print('cd $PBS_O_WORKDIR')
	print()

	print('module purge')
	print('module load python/2.7.13-anaconda')
	print('module load blast/2.2.29+')

	print()

	print('# Run transrate')
	print('transrate --assembly', args.fasta, '--left', left, '--right', right, '--threads', args.cpus, \
		'--reference /storage/aja/diatom_protein_dbs/phaeoAA.fa', '--output', outdir)
	
	print()
	print('# Run filter_transcripts_transrate.py')
	print('python2 /home/aja/local/src/yang-2021/phylogenomic_dataset_construction/scripts/filter_transcripts_transrate.py', args.fasta, outdir + '/' + basename + '/' + \
		'contigs.csv . >', basename + '.filtersum')


def main():
	print_pbs()

# get the arguments before calling main
args = get_args()

# execute the program by calling main
if __name__ == "__main__":
	main()


