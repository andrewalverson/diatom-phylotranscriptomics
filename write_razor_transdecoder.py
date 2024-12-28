#! /usr/bin/env python3

# load required modules
import argparse

def get_args():
	# create an ArgumentParser object ('parser') that will hold all the information necessary to parse the command line
	parser = argparse.ArgumentParser(description='Runs the final TransDecoder prediction')

	# add positional arguments
	parser.add_argument( "fasta", help='trinity assembly FASTA file name')

	# add optional arguments
	parser.add_argument( "--queue", "-q", help='razor queue', default='aja')
	parser.add_argument( "--wall",  "-w", help='walltime', type=int, default=72)
	parser.add_argument( "--cpus",  "-c", help='number of cpus', type=int, default='16')

	return parser.parse_args()


def print_pbs():

	job = args.fasta.replace('largest_cluster_transcripts.fa', 'transdecoder')
	input_directory = args.fasta + '.transdecoder_dir'


	# mmetsp
	# pfam   = args.fasta.replace('good_transcripts.short_name.filtered_transcripts.largest_cluster_transcripts.fa', 'fa.hmmer.pfam.domtblout')
	# uniref = args.fasta.replace('good_transcripts.short_name.filtered_transcripts.largest_cluster_transcripts.fa', 'fa.vs.uniref.blastp')

	# alverson_lab and NIES-3581
	pfam   = args.fasta.replace('largest_cluster_transcripts.fa', 'fa.hmmer.pfam.domtblout')
	uniref = args.fasta.replace('largest_cluster_transcripts.fa', 'fa.vs.uniref.blastp')

	print('#PBS -N', job)
	print('#PBS -q', args.queue)
	print('#PBS -o',  job + '.out')
	print('#PBS -l nodes=1:ppn=' + str(args.cpus))
	print('#PBS -l walltime=' + str(args.wall) + ':00:00')
	print()

	print('cd $PBS_O_WORKDIR')
	print()


	print('module load python/3.7.3-anaconda-razor')
	print('module load blast/2.3.0+')
	print('module load hmmer/3.1b2')

	print()
	print('conda init bash')
	print('source activate /home/aja/.conda/envs/transdecoder')


	print()
	print('# TransDecoder step 1: extract long orfs')
	print('TransDecoder.LongOrfs -t', args.fasta)

	print()
	print('# TransDecoder step 2.1: BLAST search to uniref')
	print('/share/apps/bioinformatics/blast/ncbi-blast-2.6.0+/bin/blastp -outfmt 6 -max_target_seqs 1 -evalue 1e-3 -db /storage/aja/swissprot_db/uniprot_sprot.fasta', \
		'-num_threads 32 -query', input_directory + '/' + 'longest_orfs.pep', '>', uniref)

	print()
	print('# TransDecoder step 2.2: HMMR search search long ORFs to Pfam database')
	print('/share/apps/bioinformatics/hmmer/hmmer-3.1b2/src/hmmscan --cpu', args.cpus, '--domtblout', pfam,  \
		'/storage/aja/pfam_db/Pfam-A.hmm', input_directory + '/' + 'longest_orfs.pep')

	print()
	print('# TransDecoder step 3: Use Transdecoder to translate nuclear contigs using homology search information')
	print('TransDecoder.Predict --single_best_only -t', args.fasta, \
		'--retain_pfam_hits', pfam, '--retain_blastp_hits', uniref)
	

def main():
	print_pbs()

# get the arguments before calling main
args = get_args()

# execute the program by calling main
if __name__ == "__main__":
	main()


