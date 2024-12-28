#! /usr/bin/env python3

# load required modules
import argparse
import re

# create an ArgumentParser object ('parser') that will hold all the information necessary to parse the command line
parser = argparse.ArgumentParser(description="OUTPUT: number of sequences, min length, max length, mean length, median length, median of upper half of sequence lengths")

# add arguments
parser.add_argument( "fasta", help="FASTA file" )

args = parser.parse_args()

genome_prefixes = ['Aureococcusanoph', 'Cyclotella_cryptica_CCMP332', 'Ectocarpussil', \
'Fistuliferasol', 'Fragilariopsiscyl_jgi', 'Nannochloropsisgad', 'NITZ4', 'Phaeodactylum_tricornutum_Phatr3', \
'Psammoneisjap', 'Pseudonitzschiamultiser_jgi', 'Pseudonitzschiamultistr_PSNMU-V1.4', 'Semro1', \
'Thalassiosira_oceanica_CCMP1005_mRNA', 'Thalassiosirapseu']


# parse FASTA file
with open(args.fasta,'r') as infile:
	for line in infile:
		line = line.rstrip()

		if(re.search('^>', line)):
			# remove '>'
			line = line[1:]

			# Trinity transcripts
			if(re.search('TRINITY', line)):
				a = line.split('_TRINITY')
				print(line, a[0])

			# genome taxon
			else:
				for i in genome_prefixes:
					if(re.search(i, line)):
						b = line.split(i)
						print(line, i)
						break
		else:
			continue


'''
# parse key file
with open(args.key,'r') as keyfile:
	for line in keyfile:
		line = line.rstrip()
		genome_prefixes.append(line)
'''
