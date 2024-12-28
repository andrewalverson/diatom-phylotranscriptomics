#! /usr/bin/env python3

# load required modules
import argparse
import os,sys,re
from Bio import SeqIO

def get_args():
	# create an ArgumentParser object ('parser') that will hold all the information necessary to parse the command line
	parser = argparse.ArgumentParser(description="this script reads a FASTA file and prints the taxon portion of each FASTA description; for the diatom-287 dataset")

	# add arguments
	parser.add_argument( "fasta", help="Concatenated FASTA file of sequences from diatom-287 dataset" )

	parser.add_argument( "-m", "--astral_mapping", help="output an ASTRAL-Pro mapping file (default=false)", action='store_true' )


	# parse the arguments
	return parser.parse_args()


def get_name(label):
	"""Get taxonID from a tip label"""

	# get rid of single quotes surrounding taxon name
	label = label.replace("'","")

	# split labels as they appear in our dataset

	rnaseq_re = re.compile("(\S+)_TRINITY_\S+$")
	aureo     = re.compile("Aureococcusanoph")
	cryptica  = re.compile("Cyclotella_cryptica_CCMP332")
	ecto      = re.compile("Ectocarpussil")
	fist      = re.compile("Fistuliferasol")
	fracy     = re.compile("Fragilariopsiscyl_jgi")
	nanno     = re.compile("Nannochloropsisgad")
	nitz4     = re.compile("NITZ4")
	phaeo     = re.compile("Phaeodactylum_tricornutum_Phatr3")
	psamm     = re.compile("Psammoneisjap")
	multiser  = re.compile("Pseudonitzschiamultiser_jgi")
	multistr  = re.compile("Pseudonitzschiamultistr_PSNMU-V1.4")
	semi      = re.compile("Semro1")
	oceanica  = re.compile("Thalassiosira_oceanica_CCMP1005_mRNA")
	thaps     = re.compile("Thalassiosirapseu")
	bolido    = re.compile("Bolidomonas_pacifica_CCMP1866")

	rnaseq_match   = re.match(rnaseq_re, label)
	aureo_match    = re.match(aureo, label)
	cryptica_match = re.match(cryptica, label)
	ecto_match     = re.match(ecto, label)
	fist_match     = re.match(fist, label)
	fracy_match    = re.match(fracy, label)
	nanno_match    = re.match(nanno, label)
	nitz4_match    = re.match(nitz4, label)
	phaeo_match    = re.match(phaeo, label)
	psamm_match    = re.match(psamm, label)
	multiser_match = re.match(multiser, label)
	multistr_match = re.match(multistr, label)
	semi_match     = re.match(semi, label)
	oceanica_match = re.match(oceanica, label)
	thaps_match    = re.match(thaps, label)
	bolido_match   = re.match(bolido, label)

	if rnaseq_match:
		# print rnaseq_match.group(1)
		return rnaseq_match.group(1)
	elif aureo_match:
		return "Aureococcusanoph"
	elif cryptica_match:
		return "Cyclotella_cryptica_CCMP332"
	elif ecto_match:
		return "Ectocarpussil"
	elif fist_match:
		return "Fistuliferasol"
	elif fracy_match:
		return "Fragilariopsiscyl_jgi"
	elif nanno_match:
		return "Nannochloropsisgad"
	elif nitz4_match:
		return "NITZ4"
	elif phaeo_match:
		return "Phaeodactylum_tricornutum_Phatr3"
	elif psamm_match:
		return "Psammoneisjap"
	elif multiser_match:
		return "Pseudonitzschiamultiser_jgi"
	elif multistr_match:
		return "Pseudonitzschiamultistr_PSNMU-V1.4"
	elif semi_match:
		return "Semro1"
	elif oceanica_match:
		return "Thalassiosira_oceanica_CCMP1005_mRNA"
	elif thaps_match:
		return "Thalassiosirapseu"
	elif bolido_match:
		return "Bolidomonas_pacifica_CCMP1866"
	else:
		# this shouldn't happen
		print("other:", label.split("_")[0], label)
		# return label.split("_")[0]
		

def main():
	fasta_sequences = SeqIO.parse(open(args.fasta),'fasta')

	for record in fasta_sequences:
		taxon = get_name(record.description)
		if(args.astral_mapping):
			print(' '.join([record.description,taxon]))
		else:
			print(taxon)

# get the command-line arguments
args = get_args()

# execute the program by calling main
if __name__ == '__main__':
	main()
