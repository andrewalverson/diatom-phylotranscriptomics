#! /usr/bin/env python3

# load required modules
import argparse
import os
import sys
import re
import statistics
import collections
from Bio import SeqIO


def get_args():
	# create an ArgumentParser object ('parser') that will hold all the information necessary to parse the command line
	parser = argparse.ArgumentParser(description="This script makes parses the #ColumnsMap list of amino acid sites output by trimal to produce a list of codon sites to for trimming from a corresponding codon alignment")

	# add arguments
	parser.add_argument( "cds", help="FASTA alignment (codons)" )

	# parse the arguments
	return parser.parse_args()


def get_alignment_length(alignment_path, alignment):
	fasta_sequences  = list(SeqIO.parse(open(alignment_path + alignment),'fasta'))
	alignment_length = len(fasta_sequences[0].seq)
	return alignment_length


def compare_aa_cds_alignments(aa_alignment_length, cds_alignment_length, aa_alignment, cds_alignment):
	if  cds_alignment_length != (3 * aa_alignment_length):
		error_message = aa_alignment + ' (' + str(aa_alignment_length) + ') and ' + cds_alignment + ' (' + str(cds_alignment_length) + ') lengths do not correspond'
		sys.exit(error_message)


def parse_trimal(base_filename, aa_alignment_length, trimal_path):

	# name of trimal column mapping file (pal2nal)
	# trimal = base_filename + '_trimal.sites'

	# name of trimal column mapping file (translatorX)
	# trimal = base_filename + '_trimal_aa.sites'

	# name of trimal column mapping file (translatorX -- without internal stop codons)
	trimal = base_filename + '_in_frame.aa_ali_trimal.sites'

	if os.path.exists(trimal_path + trimal):

		with open(trimal_path + trimal, 'r') as trimal_mapping:
			untrimmed_aa_columns = trimal_mapping.read()
			untrimmed_aa_columns = untrimmed_aa_columns.rstrip()

		# remove trimal line header
		untrimmed_aa_columns = re.sub(r'#ColumnsMap\s+', r'', untrimmed_aa_columns)
		# print(untrimmed_columns)
		# print(alignment_length)

		# covert untrimmed_columns from string to list
		untrimmed_aa_columns_list = untrimmed_aa_columns.split(', ')
		# for i in untrimmed_columns_list:
		# 	print(i)

		# define list to hold the trimmed columns
		trimmed_aa_columns_list = []
		
		# from the list of *retained* columns from trimal, get the inverse list of *trimmed* columns
		for i in range(0, aa_alignment_length):
			if str(i) in untrimmed_aa_columns_list:
				continue
			else:
				trimmed_aa_columns_list.append(str(i))
		

		return trimmed_aa_columns_list


def aa2cds_columns(aa_length, trimmed_aa_columns):
	# a list of codon columns to trim
	trimmed_cds_columns = []

	# print(coords_to_ranges(trimmed_aa_columns))
	# loop over list of AA columns to trim, convert AA --> CDS column coords, then add those to trimmed_cds_columns[]
	for i in range(0, aa_length):
		if str(i) in trimmed_aa_columns:
			trimmed_cds_columns.extend([i*3, (i*3)+1, (i*3)+2])

	# convert list of coordinates to ranges, and return this to main()
	return coords_to_ranges(trimmed_cds_columns)


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



def write_PBS(base_filename, cds_selectcols, alignment_path, out_path):

	# set walltime and number of cpu's
	wall    = '06'
	num_cpu = '1'

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
	trimmed_cds = base_filename + '_alignment_trimmed.cds'

	print('/share/apps/bioinformatics/trimal/1.4rev22/trimal', '-selectcols', '{ ' + cds_selectcols + ' }', '-keepheader', '-colnumbering', '-in', alignment_path + args.cds, '-out', out_path + trimmed_cds, '>', out_path + base_filename + '_alignment_trimmed.sites', '2>', out_path + base_filename + '_alignment_trimmed.err')


def main():
	# parse name of alignment file in order to set the base filename for output files
	a = args.cds.split("_")
	base_filename = a[0]

	"""
	# name of untrimmed AA alignment (for original pal2nal alignments)
	aa_alignment = base_filename + '_alignment.fasta'

	set paths to input and output files (for pal2nal alignments)
	aa_alignment_path  = '/storage/aja/orthofinder/orthogroups/orthogroup_alignments_occupancy_10_round2/AA_alignments_untrimmed/'
	cds_alignment_path = '/storage/aja/orthofinder/orthogroups/orthogroup_alignments_occupancy_10_round2/CDS_aligned_untrimmed/'
	trimal_path        = '/storage/aja/orthofinder/orthogroups/orthogroup_alignments_occupancy_10_round2/AA_alignments_masked_gappyout/'
	out_path           = '/storage/aja/orthofinder/orthogroups/orthogroup_alignments_occupancy_10_round2/CDS_trimmed/'

	# name of untrimmed AA alignment (for translatorX alignments -- with internal stops version)
	aa_alignment = base_filename + '_in_frame.aa_ali.fasta'

	set paths to input and output files (for translatorX alignments)
	aa_alignment_path  = '/storage/aja/orthofinder/orthogroups/orthogroup_alignments_occupancy_10_round2/CDS_unaligned/with_internal_stops/'
	cds_alignment_path = '/storage/aja/orthofinder/orthogroups/orthogroup_alignments_occupancy_10_round2/CDS_unaligned/with_internal_stops/'
	trimal_path        = '/storage/aja/orthofinder/orthogroups/orthogroup_alignments_occupancy_10_round2/CDS_unaligned/with_internal_stops/'
	out_path           = '/storage/aja/orthofinder/orthogroups/orthogroup_alignments_occupancy_10_round2/CDS_unaligned/with_internal_stops/'

	# name of untrimmed AA alignment (for failed pal2nal caused by mismatches in number of sequences between AA and CDS files, which traces back to UPP removing some sequences)
	aa_alignment = base_filename + '_alignment_subsetted.fasta'

	# set paths to input and output files (for failed pal2nal alignments)
	aa_alignment_path  = '/storage/aja/orthofinder/orthogroups/orthogroup_alignments_occupancy_10_round2/pal2nal/pal2nal_failed/'
	cds_alignment_path = '/storage/aja/orthofinder/orthogroups/orthogroup_alignments_occupancy_10_round2/pal2nal/pal2nal_failed/'
	trimal_path        = '/storage/aja/orthofinder/orthogroups/orthogroup_alignments_occupancy_10_round2/AA_alignments_masked_gappyout/'
	out_path           = '/storage/aja/orthofinder/orthogroups/orthogroup_alignments_occupancy_10_round2/pal2nal/pal2nal_failed/'
	"""

	# name of untrimmed AA alignment (for translatorX alignments -- without internal stops version))
	aa_alignment = base_filename + '_in_frame.aa_ali.fasta'

	# set paths to input and output files (for translatorX alignments)
	aa_alignment_path  = '/storage/aja/orthofinder/orthogroups/orthogroup_alignments_occupancy_10_round2/CDS_unaligned/without_internal_stops/translatorX/'
	cds_alignment_path = '/storage/aja/orthofinder/orthogroups/orthogroup_alignments_occupancy_10_round2/CDS_unaligned/without_internal_stops/translatorX/'
	trimal_path        = '/storage/aja/orthofinder/orthogroups/orthogroup_alignments_occupancy_10_round2/CDS_unaligned/without_internal_stops/translatorX/'
	out_path           = '/storage/aja/orthofinder/orthogroups/orthogroup_alignments_occupancy_10_round2/CDS_unaligned/without_internal_stops/translatorX/'


	aa_alignment_length  = get_alignment_length(aa_alignment_path, aa_alignment)
	cds_alignment_length = get_alignment_length(cds_alignment_path, args.cds)
	compare_aa_cds_alignments(aa_alignment_length, cds_alignment_length, aa_alignment, args.cds)
	trimmed_aa_columns   = parse_trimal(base_filename, aa_alignment_length, trimal_path)
	cds_selectcols       = aa2cds_columns(aa_alignment_length, trimmed_aa_columns)
	write_PBS(base_filename, cds_selectcols, cds_alignment_path, out_path)


# get the command-line arguments
args = get_args()

# execute the program by calling main
if __name__ == '__main__':
	main()

