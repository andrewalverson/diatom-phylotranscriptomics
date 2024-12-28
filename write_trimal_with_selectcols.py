#! /usr/bin/env python3

# load required modules
import argparse
import os, sys
import re
import statistics
import collections
from Bio import SeqIO


def get_args():
	# create an ArgumentParser object ('parser') that will hold all the information necessary to parse the command line
	parser = argparse.ArgumentParser(description="This script makes a global list of masked sites (from UPP) and trimmed sites (from trimal) to trim from an alignment")

	# add arguments
	parser.add_argument( "alignment", help="original FASTA alignment" )

	# parse the arguments
	return parser.parse_args()


def get_alignment_length(alignment_path):
	# fasta_sequences  = list(SeqIO.parse(open(alignment_path + args.alignment),'fasta'))
	fasta_sequences  = list(SeqIO.parse(open(args.alignment),'fasta'))
	alignment_length = len(fasta_sequences[0].seq)
	return alignment_length


def parse_trimal(base_filename, alignment_length, trimal_path):

	# name of trimal column mapping file
	# trimal = base_filename + '_alignment.sites'
	trimal = base_filename + '_trimal_gt0.1.sites'

	if os.path.exists(trimal):

		# with open(trimal_path + trimal, 'r') as trimal_mapping:
		with open(trimal, 'r') as trimal_mapping:
			untrimmed_columns = trimal_mapping.read()
			untrimmed_columns = untrimmed_columns.rstrip()

		# remove trimal line header
		untrimmed_columns = re.sub(r'#ColumnsMap\s+', r'', untrimmed_columns)
		# print(untrimmed_columns)
		# print(alignment_length)

		# covert untrimmed_columns from string to list
		untrimmed_columns_list = untrimmed_columns.split(', ')
		# for i in untrimmed_columns_list:
		# 	print(i)

		# define list to hold the trimmed columns
		trimmed_columns_list = []
		
		# from the list of *retained* columns from trimal, get the inverse list of *trimmed* columns
		for i in range(0, alignment_length):
			if str(i) in untrimmed_columns_list:
				continue
			else:
				trimmed_columns_list.append(str(i))
		return coords_to_ranges(trimmed_columns_list)
	else:
		sys.exit()


def coords_to_ranges(coords):
	starts = collections.OrderedDict()
	ends   = collections.OrderedDict()
	for idx, page in enumerate(coords):
		section = int(page) - int(idx)
		starts.setdefault(section, page)
		ends[section] = page
	page_parts = []
	for section, start in starts.items():
		end = ends[section]
		if start == end:
			page_parts.append("{0}".format(start))
		else:
			page_parts.append("{0}-{1}".format(start, end))
	# print(','.join(page_parts))
	return(','.join(page_parts))


def parse_upp_masked_columns(base_filename, alignment_path):

	# path to alignments and masked sites
	# alignment_path = '/storage/aja/orthofinder/orthogroups/orthogroup_alignments_occupancy_10/AA_alignments_correct_headers/'  # round 1
	# alignment_path   = '/storage/aja/orthofinder/orthogroups/orthogroup_alignments_occupancy_10_round2/UPP_M_setting/'           # round 2
	
	upp_masked = base_filename + '_insertion_columns.txt'

	# if os.path.exists(alignment_path + upp_masked):
	if os.path.exists(upp_masked):

		# with open(alignment_path + upp_masked, 'r') as masked:
		with open(upp_masked, 'r') as masked:
			masked_columns = masked.read()
		
		# return comma separated list of columns masked by UPP
		return(masked_columns.rstrip())
	else:
		sys.exit()


def print_command(base_filename, trimmed_columns, masked_columns, alignment_path, out_path):

	# compile trimmed_columns string, accounting for possibility that one or both of masked_columns and trimmed_columns might be empty
	if masked_columns:
		if trimmed_columns:
			selectcols = masked_columns + ',' + trimmed_columns
		else:
			selectcols = masked_columns
	else:
		if trimmed_columns:
			selectcols = trimmed_columns	
		else:
			selectcols = ''

	# set walltime and number of cpu's
	# wall    = '06'
	# num_cpu = '1'

	# print top of PBS script
	# print('#PBS -N', base_filename+'.trimal-AA')
	# print('#PBS -q q06h32c')
	# print('#PBS -j oe')
	# print('#PBS -o',  base_filename+'.trimal-AA.$PBS_JOBID')
	# print('#PBS -l nodes=1:ppn=' + num_cpu)
	# print('#PBS -l walltime=' + wall + ':00:00')
	# print()

	# print('cd $PBS_O_WORKDIR')
	# print()

	# print trimal command
	# print('# run trimal')
	# print('/share/apps/bioinformatics/trimal/1.4rev22/trimal', '-selectcols', '{ ' + selectcols + ' }', '-keepheader', '-colnumbering', '-in', alignment_path + args.alignment, '-out', out_path + args.alignment, '>', out_path + base_filename + '_trimal.sites', '2>', out_path + base_filename + '_alignment.err')
	print('/share/apps/bioinformatics/trimal/1.4rev22/trimal', '-selectcols', '{ ' + selectcols + ' }', '-keepheader', '-colnumbering', '-in', args.alignment, '-out', base_filename + '_trimal_final.fa', '>', base_filename + '_trimal_final.sites', '2>', base_filename + '_trimal_final.err')


def main():
	# parse name of alignment file in order to set the base filename for output files
	# a = args.alignment.split("_")
	a = args.alignment.split("_alignment.fasta")
	base_filename = a[0]

	# print(base_filename + '\n')

	# set paths to input and output files (ROUND 1)
	# alignment_path = '/storage/aja/orthofinder/orthogroups/orthogroup_alignments_occupancy_10/AA_alignments_correct_headers/' # round 1
	# trimal_path    = '/storage/aja/orthofinder/orthogroups/orthogroup_alignments_occupancy_10/AA_alignments_trimal_gappyout/' # round 1
	# out_path       = '/storage/aja/orthofinder/orthogroups/orthogroup_alignments_occupancy_10/AA_alignments_masked_gappyout/' # round 1

	# set paths to input and output files (ROUND 2)
	# alignment_path = '/storage/aja/orthofinder/orthogroups/orthogroup_alignments_occupancy_10_round2/AA_alignments_untrimmed/'       # round 2
	# trimal_path    = '/storage/aja/orthofinder/orthogroups/orthogroup_alignments_occupancy_10_round2/AA_alignments_trimal_gappyout/' # round 2
	# out_path       = '/storage/aja/orthofinder/orthogroups/orthogroup_alignments_occupancy_10_round2/AA_alignments_masked_gappyout/' # round 2

	# set paths to input and output files (ROUND 2)
	alignment_path = '/storage/aja/orthofinder/gene_trees_round_4_ortholog_selection_and_trees/RT_ortholog_realignments/' # round 4
	trimal_path    = '/storage/aja/orthofinder/gene_trees_round_4_ortholog_selection_and_trees/RT_ortholog_realignments/' # round 4
	out_path       = '/storage/aja/orthofinder/gene_trees_round_4_ortholog_selection_and_trees/RT_ortholog_realignments/' # round 4

	alignment_length = get_alignment_length(alignment_path)
	trimmed_columns  = parse_trimal(base_filename, alignment_length, trimal_path)
	masked_columns   = parse_upp_masked_columns(base_filename, alignment_path)
	print_command(base_filename, trimmed_columns, masked_columns, alignment_path, out_path)


# get the command-line arguments
args = get_args()

# execute the program by calling main
if __name__ == '__main__':
	main()

