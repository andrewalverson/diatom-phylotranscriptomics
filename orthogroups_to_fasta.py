#! /usr/bin/env python3

# load required modules
import sys
import os
import re
import argparse
from Bio import SeqIO

# create an ArgumentParser object ('parser') that will hold all the information necessary to parse the command line
parser = argparse.ArgumentParser(description="Create FASTA files from OrthoFinder's 'Orthogroups.csv' output file; FASTA files are written to the working directory")

# add arguments
parser.add_argument( "ortho", help="'Orthogroups.csv' file" )
parser.add_argument( "fasta", help="Directory containing FASTA files for all species in the OrthoFinder analysis" )
parser.add_argument( "-s", "--fasta_suffix", default='fa', help="file suffix for individual FASTA files in the 'fasta' directory (default = 'fa')")
# parser.add_argument( "outdir", help="Path to output directory" )

args = parser.parse_args()

# exit if Orthogroups.csv or the specified directories don't exist
os.path.isdir(args.fasta)  or sys.exit("The directory", args.fasta, "does not exist – exiting.")
os.path.isfile(args.ortho) or sys.exit("The file", args.ortho, "does not exist – exiting.")
# os.path.isdir(args.outdir) or sys.exit("The directory", args.outdir, "does not exist – exiting.")

# check for trailing '/' on the path to the RNAseq assemblies, add if not there
if not args.fasta.endswith('/'):
	args.fasta += '/'
# if not args.outdir.endswith('/'):
# 	args.outdir += '/'

fasta_dictionary = {} # holds the full list of species names, parsed from file names of the FASTA files; key = species name, value = assembly FASTA file
taxon_numbers    = {} # key = taxon number (order/index as listed in Orthogroups.csv), value = taxon name

# read in all the contents of the FASTA directory
fasta_files = os.listdir(args.fasta)

# loop over all assembly fasta files in the directory, parsing and storing the species name, as gleaned from the filename
# key = species name (from filename), value = assembly file name for this species
for fi in fasta_files:
	filename_re = '(.*)\.' + args.fasta_suffix + '$' # files with a FASTA suffix
	match = re.search(filename_re, fi) # only consider *.fa[sta] files
	if match:
		fasta_dictionary[(match.group(1))] = match.group(0)
	else:
		continue

# read and parse the orthogroups file
with open(args.ortho, 'r') as og:

	# loop over each line in the file: first line = taxon IDs, following lines = orthogroups
	for line in og:
		line = line.rstrip()

		# header line begsin with a tab
		if re.search("^\t", line):
			taxon_list = line.split()
			
			# loop over list of taxa; key = taxon number, value = taxon name
			for i in range(len(taxon_list)):
				taxon_numbers[i] = taxon_list[i]

		# else it's an orthogroup
		else: 
			# split on tab to isolate orthogroup ID (first element), with each subsequent element a csv list of all the genes for a given species in this orthogroup
			species = line.split('\t')

			# set name for output FASTA file -> orthogroup_name.fa
			# output_filename =  args.ortho + species.pop(0) + '.fa'
			output_filename =  species.pop(0) + '.fa'

			# open the output file for writing
			output_handle = open(output_filename, 'w')

			# loop over all the groups/species, where each group/species is either:
				# csv list of sequences for a given species (matches ', ')
				# single sequence for a given species (doesn't match ', ')
				# an empty string, i.e., no seqs for this species (matches '')
			# loop over the indices so I can index to get the RNAseq filename for this species from fasta_dictionary{}
			for j in range(len(species)):
				# some species will have no sequences in an orthogroup – skip these
				if species[j] == '':
					continue
				# this species has >=1 sequences in this orthogroup
				else:
					# open and parse the RNAseq assembly for this species, store in a dictionary for random record access
					path_to_assembly = args.fasta + fasta_dictionary[taxon_numbers[j]]
					assembly_dict = SeqIO.to_dict(SeqIO.parse(path_to_assembly,'fasta'))

					# if there's a match to ' ,', that means there's >1 sequence for this species - split genes into a list
					if re.search(', ', species[j]):
						# split all the genes for this species into a list called 'genes_for_this_species'
						genes_for_this_species = species[j].split(', ')

						# loop over all the genes (in this orthogroup) for this species
						for paralog in genes_for_this_species:

							# check whether the desired sequence is in the file – print to STDOUT if it's not
							if paralog in assembly_dict:
								# write the sequence to a file
								SeqIO.write(assembly_dict[paralog], output_handle, "fasta")
							else:
								print(output_filename, paralog)

					# else there's just one sequence for this species, so no comma to split on
					else:
						# check whether the desired sequence is in the file – print to STDOUT if it's not
						if species[j] in assembly_dict:
						# write the sequence to a file
							SeqIO.write(assembly_dict[species[j]], output_handle, "fasta")
						else:
							print(output_filename, species[j])




			output_handle.close()

og.close()

