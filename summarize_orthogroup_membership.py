#! /usr/bin/env python3

# load required modules
import sys
import os
import re
import csv
import argparse

# create an ArgumentParser object ('parser') that will hold all the information necessary to parse the command line
parser = argparse.ArgumentParser(description="Summarize ingroup/outgroup membership in OrthoFinder-based orthogroups from the 'Orthogroups.GeneCount.csv' file")

# add arguments
parser.add_argument( "ogs", help="OrthoFinder-based orthogroups in the 'Orthogroups.GeneCount.csv' file" )
parser.add_argument( "out", help="File containing list (one species per row) of outgroup species" )
parser.add_argument( "req", help="Report statistics for this species (e.g., a required outgroup)" )

args = parser.parse_args()

# list containing names of outgroup taxa
outgroups = []

# get list of outgroup taxa
with open(args.out, 'r') as outgroup_file:
	for line in outgroup_file:
		line = line.rstrip()
		outgroups.append(line)

# total number of outgroups
num_outgroups = len(outgroups)
# print("num_outgroups:", num_outgroups)

# print header line
s = "\t"
print(s.join(['Orthogroup', 'Unique_diatoms', 'Unique_outgroups', 'Total_species', 'Diatom_occupancy', 'Taxon_occupancy', args.req, 'Outgroup_count', 'Diatom_count', 'Total_gene_count', 'NumTaxGT100copies', 'NumTaxGT500copies', 'NumTaxGT1000copies' ]))

# dictionary to hold column numbers of outgroup taxa; key = index, value = taxon name
outgroup_indices = {}

# open and orthogroups
with open(args.ogs, 'r') as orthogroups:

	# loop over each line in the file: first line = taxon IDs, following lines = orthogroups
	for line in orthogroups:
		line = line.rstrip()

		# header line begsin with a tab
		if re.search("^\t", line):
			taxa = line.split()

			# total number of taxa (ingroup _and_ outgroup) in the matrix
			# subtract one to account for the "Total" column
			num_taxa = len(taxa) - 1
			# print("num_taxa:", num_taxa)			

			# get column numbers for outgroup taxa
			for i in outgroups:
				o = taxa.index(i) + 1  # adding +1 because the column/index for each species is one greater in the orthogroup lines -- want those to be identical
				outgroup_indices[o] = i
			
			# get column number for required outgroup/special taxon (args.req)
			p = taxa.index(args.req) + 1  # adding +1 because the column/index for each species is one greater in the orthogroup lines -- want those to be identical
			outgroup_indices[p] = args.req

			# print(taxa.index(args.req), outgroup_indices[taxa.index(args.req)])
			# for num, name in outgroup_indices.items():
			# 	print(name, "\t", num)

		# not the header line, so it's an orthogroup line
		else:

			# dictionary to hold number of copies per species per orthogroup; key = taxon number, value = num copies in this orthogroup
			copy_number_per_species = {}
			
			# for storing counts of number of taxa and genes per orthogroup
			diatom_unique   = 0 # number of unique diatom species in this orthogroup
			diatom_total    = 0 # total number of diatom sequences in this orthogroup
			total_species   = 0 # number of unique species (ingroup + outgroup) in this orthogroup
			outgroup_unique = 0 # number of unique outgroup species in this orthogroup
			outgroup_total  = 0 # total number of outgroup sequences in this orthogroup
			req_count       = 0 # total number of sequences in this orthogroup for the required outgroup taxon
			gene_count      = 0 # total number of sequences in this orthogroup (i.e., the orthogroup size)

			# split the line
			counts = line.split()

			# loop over all fields in this line, where a field is the gene count for each species
			for j in range(len(counts)):

				# first index = orthogroup ID -- print it
				if j == 0:
					print(counts[j], end='\t')
					continue
				
				# last index = total gene count -- record it
				elif j == (len(counts)-1):
					gene_count = counts[j]
					# print(gene_count)

				# else its the ortho count for some species
				else:
					# keep track of total copy count for this taxon in copy_number_per_species dictionary
					copy_number_per_species[j] = int(counts[j])

					# is this index (i.e., column number, where each column is a species) in the outgroup_indices dictionary?
					# i.e., is this an outgroup column?
					if(j in outgroup_indices):

						# does this particular outgroup taxon have a sequence in this orthogroup?
						if int(counts[j]) > 0:
							outgroup_unique += 1
							outgroup_total  += int(counts[j])
							total_species   += 1

						# is this outgroup sequence in the required/special outgroup taxon?
						if outgroup_indices[j] == args.req:
							# record the number of sequences in this orthogroup for the required/special taxon (= args.req)
							req_count += int(counts[j])
					# not an outgroup, so this is a diatom
					else:
						# does this diatom taxon have a sequence in this orthogroup?
						if int(counts[j]) > 0:
							diatom_unique += 1
							diatom_total  += int(counts[j])
							total_species += 1


			# calculate number of species with sequence/gene counts greater than 100, 500, and 1000
			gt100  = 0
			gt500  = 0
			gt1000 = 0

			for tax_num, copy_number in copy_number_per_species.items():
				if copy_number > 100:
					gt100 += 1
				elif copy_number > 500:
					gt500 += 1
				elif copy_number > 1000:
					gt1000 += 1

			# check whether my gene count matches the OrthoFinder count for this orthogroup
			# if they match, print output
			#['Orthogroup', 'Unique_diatoms', 'Unique_outgroups', 'Total_species', 'Diatom_occupancy', 'Taxon_occupancy', 'Bolidomonas_count', 'Outgroup_count', 'Diatom_count', 'Total_gene_count', 'NumTaxGT100copies', 'NumTaxGT500copies', 'NumTaxGT1000copies' ]
			if int(gene_count) == (diatom_total + outgroup_total):
				print(diatom_unique,   end='\t') # total number of unique diatom taxa
				print(outgroup_unique, end='\t') # total number of unique outgroups
				print(total_species, end='\t')   # total number of unique species (diatoms + outgroups)
				print( "%.2f" % (diatom_unique/(num_taxa - num_outgroups)),   end='\t') # taxon occupancy (diatoms)
				print( "%.2f" % ((diatom_unique + outgroup_unique)/num_taxa), end='\t') # taxon occupancy (total)
				print(req_count, end='\t') # total number of sequences in the required outgroup taxon (args.req)
				print(outgroup_total, end='\t') # total number of outgroup sequences
				print(diatom_total, end='\t') # total number of diatom sequences
				print(gene_count, end='\t') # total number of sequences in this orthogroup
				print(gt100, end='\t')  # total number of species with >100  sequences in this orthogroup
				print(gt500, end='\t')  # total number of species with >500  sequences in this orthogroup
				print(gt1000) # total number of species with >1000 sequences in this orthogroup

			# if my total gene count doesn't match OrthoFinders, thrown an error
			else:
				print("OrthoFinder count: ", gene_count)
				print("My count: ", diatom_total+outgroup_total)
				sys.exit()


