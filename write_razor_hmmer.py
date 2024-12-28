#! /usr/bin/env python3

# load required modules
import argparse

def get_args():
	# create an ArgumentParser object ('parser') that will hold all the information necessary to parse the command line
	parser = argparse.ArgumentParser(description='this script parses output from `check_for_stops.pl`, identifies out-of-frame CDS sequences, then uses BLASTX to search the problematic CDS sequence against the corresponding AA sequence and replace the out-of-frame with the in-frame CDS sequence')

	# add positional arguments
	parser.add_argument( "fasta", help='trinity assembly FASTA file name')

	# add optional arguments
	parser.add_argument( "--jobID", "-j", help='job name', default='aja')
	parser.add_argument( "--queue", "-q", help='razor queue', default='aja')
	parser.add_argument( "--wall",  "-w", help='walltime', type=int, default=72)
	parser.add_argument( "--cpus",  "-c", help='number of cpus', type=int, default='16')

	return parser.parse_args()


def parse_names():
	input_directory = args.fasta + '.transdecoder_dir'
	output_file = args.fasta + '.pfam.domtblout'

	return(input_directory, output_file)

def print_pbs():

	input_directory, output_file = parse_names()
	jobID = args.fasta + '.hmmer.pfam.domtblout'



	print('#PBS -N', jobID)
	print('#PBS -q', args.queue)
	print('#PBS -o',  jobID)
	print('#PBS -l nodes=1:ppn=' + str(args.cpus))
	print('#PBS -l walltime=' + str(args.wall) + ':00:00')
	print()

	print('cd $PBS_O_WORKDIR')
	print()

	print('module load blast/2.6.0+')
	print('module load perl/5.10.1')
	print('module load hmmer/3.1b2')
	print('module load transdecoder/2.0.1')
	print('module load cdhit/4.6.8')

	print()

	print('# TransDecoder step 2.2: HMMR search search long ORFs to Pfam database')
	print('hmmscan --cpu', args.cpus, '--domtblout', '$PBS_JOBID' + '.out',  \
		'/storage/aja/pfam_db/Pfam-A.hmm', input_directory + '/' + 'longest_orfs.pep')


def main():
	print_pbs()

# get the arguments before calling main
args = get_args()

# execute the program by calling main
if __name__ == "__main__":
	main()


