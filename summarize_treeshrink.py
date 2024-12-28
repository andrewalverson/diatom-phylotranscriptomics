#! /usr/bin/env python3

import argparse
import re
import sys,os
import csv

# create an ArgumentParser object ('parser') that will hold all the information necessary to parse the command line
parser = argparse.ArgumentParser(description = "summarize the tips removed by TreeShrink" )

# add arguments
parser.add_argument( "ending", help="file ending for summary files" )
parser.add_argument( "directory", help="working directory" )

args = parser.parse_args()

def get_taxon(name):
	
	rnaseq_re = re.compile("(\S+)_TRINITY")
	aureo     = re.compile("Aureococcusanoph")
	cryptica  = re.compile("CCMP332_Cyclotella_cryptica")
	ecto  = re.compile("Ectocarpussil")
	fist  = re.compile("Fistuliferasol")
	fracy = re.compile("Fragilariopsiscyl_jgi")
	nanno = re.compile("Nannochloropsisgad")
	nitz4 = re.compile("NITZ4")
	phaeo = re.compile("Phaeodactylum_tricornutum_Phatr3")
	psamm    = re.compile("Psammoneisjap")
	multiser = re.compile("Pseudonitzschiamultiser_jgi")
	multistr = re.compile("Pseudonitzschiamultistr_PSNMU-V1.4")
	semi     = re.compile("Semro1")
	oceanica = re.compile("Thalassiosira_oceanica_CCMP1005_mRNA")
	thaps    = re.compile("Thalassiosirapseu")
	bolido   = re.compile("Bolidomonas_pacifica_CCMP1866")

	rnaseq_match   = re.match(rnaseq_re, name)
	aureo_match    = re.match(aureo, name)
	cryptica_match = re.match(cryptica, name)
	ecto_match     = re.match(ecto, name)
	fist_match     = re.match(fist, name)
	fracy_match    = re.match(fracy, name)
	nanno_match    = re.match(nanno, name)
	nitz4_match    = re.match(nitz4, name)
	phaeo_match    = re.match(phaeo, name)
	psamm_match    = re.match(psamm, name)
	multiser_match = re.match(multiser, name)
	multistr_match = re.match(multistr, name)
	semi_match     = re.match(semi, name)
	oceanica_match = re.match(oceanica, name)
	thaps_match    = re.match(thaps, name)
	bolido_match   = re.match(bolido, name)

	if rnaseq_match:
		return rnaseq_match.group(1)
	elif aureo_match:
		return "OUTGROUP_Aureococcus"
	elif cryptica_match:
		return "Cyclotella_cryptica_CCMP332"
	elif ecto_match:
		return "OUTGROUP_Ectocarpus"
	elif fist_match:
		return "Fistuliferasol"
	elif fracy_match:
		return "Fragilariopsiscyl_jgi"
	elif nanno_match:
		return "OUTGROUP_Nannochloropsisgad"
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
		return "OUTGROUP_Bolidomonas_pacifica_CCMP1866"

def get_front_names(tips):
	return [get_taxon(i) for i in tips]

outfile = open("treeshrink_summary", "w")
outfile.write("taxonID\ttips_removed\n")

DIR = args.directory + "/"
DICT = {}

for i in os.listdir(DIR):
	if i[-len(args.ending):] == args.ending in i:
		if os.path.getsize(i) > 1:
			print(i)
			with open(DIR+i, "r") as infile:
				values = infile.read().split('\t')
			names = list(map(get_taxon, values))
			for taxon in names:
				if taxon is None:
					continue
				if taxon not in DICT:
					DICT[taxon] = 0
				DICT[taxon] += 1

for taxon in DICT:
	outfile.write(taxon + '\t' + str(DICT[taxon]) + '\n')
outfile.close()

