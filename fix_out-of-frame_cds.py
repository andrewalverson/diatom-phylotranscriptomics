#! /usr/bin/env python3

# this script does the following: (1) parses output from `check_for_stops.pl`, (2) identifies out-of-frame CDS sequences, (3) uses BLASTX 
# to search the problematic CDS sequence against the corresponding AA sequence, (4) extract the in-frame CDS sequence, and (5) replace
# the out-of-frame with the in-frame CDS sequence


# load required modules
import os
import sys
import argparse
import csv
from Bio import SeqIO
from collections import defaultdict


def get_args():
	# create an ArgumentParser object ('parser') that will hold all the information necessary to parse the command line
	parser = argparse.ArgumentParser(description='this script parses output from `check_for_stops.pl`, identifies out-of-frame CDS sequences, then uses BLASTX to search the problematic CDS sequence against the corresponding AA sequence and replace the out-of-frame with the in-frame CDS sequence')

	# add positional arguments
	parser.add_argument( "stops", help='output from `check_for_stops.pl`')

	# add optional arguments
	parser.add_argument( "--id", help='percent ID to consider CDS and AA sequence a match', type=float, default=98)

	return parser.parse_args()


def parse_stops(corrected_cds):
	# open and parse the check_for_stops output
	with open(args.stops, 'r') as stops_out:
		for line in csv.reader(stops_out, delimiter = '\t' ):
			
			# extract the CDS sequence and write to file
			cds_header, cds_seq = get_sequence_from_multifasta(line[0], line[1])

			# set filename for writing out this CDS sequence
			cds_fasta = line[1] + "." + "cds"

			# write out the CDS sequence
			cds_out = open(cds_fasta, 'w')
			cds_out.write(">" + cds_header + "\n")
			cds_out.write(str(cds_seq) + "\n")
			cds_out.close()

			# set name of the file that contains all the corresponding AA sequences
			aa_file  = line[0].replace("cds", "fasta")

			# get the AA sequence and FASTA header
			aa_header, aa_seq   = get_sequence_from_multifasta(aa_file, line[1])

			# set filename for writing out this AA sequence
			aa_fasta = line[1] + "." + "prot"

			# write out the AA sequence
			aa_out = open(aa_fasta, 'w')
			aa_out.write(">" + aa_header + "\n")
			aa_out.write(str(aa_seq) + "\n")
			aa_out.close()

			# use BLASTX to search the problematic CDS sequence against this datbase
			blastx_output = aa_header + ".blastx"
			blastx = "blastx -query " + cds_fasta + " -subject " + aa_fasta + " -max_hsps 1 -outfmt '6 length nident mismatch gaps qstart qend'" + " > " + blastx_output
			os.system(blastx)

			# parse the BLAST output
			parsed_blast_return = parse_blast(blastx_output, cds_seq)

			# if return value is a float, then percent identity between CDS and AA was less than args.id, so print error message
			if isinstance(parsed_blast_return, float):
				error_message = ', '.join([args.stops, cds_header, "AA and CDS percent identity: " + str(round(parsed_blast_return, 2))])
				print(error_message, file=sys.stderr)

			# else it's the in-frame CDS sequence, store it
			else:
				corrected_cds[cds_header] = parsed_blast_return
				# print(">" + cds_header)
				# print(parsed_blast_return)

			# remove BLASTX output and AA and CDS sequence files
			os.remove(blastx_output)
			os.remove(aa_fasta)
			os.remove(cds_fasta)

	return corrected_cds, line[0]


def substitute_in_frame_cds(corrected_cds, seq_file):
	# open and parse CDS sequence file
	cds_sequences = SeqIO.parse(open(seq_file),'fasta')

	for record in cds_sequences:
		if record.description in corrected_cds:
			print(">" + record.description)
			print(corrected_cds.get(record.description))
		else:
			print(record.format("fasta"), end='')


def parse_blast(blastout, cds_seq):
	# open and parse the check_for_stops output
	with open(blastout, 'r') as blast:
		for line in csv.reader(blast, delimiter = '\t' ):

			# convert all numbers in the BLAST output to integers
			line = [ int(x) for x in line ]

			# name the BLAST output fields
			length, nident, mismatch, gaps, qstart, qend = line

			# sum the number of identical matches and number of gaps
			identity_with_gaps = 100 * ((nident + gaps) / length)

			# if percent ID is high enough, extract the CDS
			if identity_with_gaps >= args.id:
				return extract_in_frame_cds(qstart, qend, cds_seq)
			else:
				return float(identity_with_gaps)


def extract_in_frame_cds(begin, end, seq):
	if(begin < end):
		# return seq[int(begin) - 1:int(end)]
		return seq[begin - 1:end]

	else:
		return seq[int(begin) - 1:int(end)].reverse_complement()


def get_sequence_from_multifasta(seq_file, header):
	# parse fasta file into a dictionay
	fasta_dict = SeqIO.index(seq_file, "fasta")

	# retrieve sequence and remove gaps
	# seq = str(fasta_dict[header].seq)
	# seq = seq.replace("-", "")

	# return the FASTA header and sequence
	return fasta_dict[header].description, fasta_dict[header].seq
	

def main():
	corrected_cds = defaultdict(dict)
	corrected_cds, cds_fasta = parse_stops(corrected_cds)
	substitute_in_frame_cds(corrected_cds, cds_fasta)


# get the arguments before calling main
args = get_args()

# execute the program by calling main
if __name__ == "__main__":
	main()
















