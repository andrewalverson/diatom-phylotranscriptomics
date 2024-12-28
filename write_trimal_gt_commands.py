#! /usr/bin/env python3

import argparse
import random
import subprocess


def get_args():
	# create an ArgumentParser object ('parser') that will hold all the information necessary to parse the command line
	parser = argparse.ArgumentParser(description = "Print trimal commands using the '-gt' option")

	# add arguments
	parser.add_argument( "fasta", help="Name of input FASTA file" )
	parser.add_argument( "-g", "--gt", help="Value for '-gt' (default = 0.1)", default = '0.1')

	return parser.parse_args()


def print_trimal():

	# parse name of input file in order to set prefix for output file names
	a = args.fasta.split(".")
	prefix = a[0]

	outfile = prefix + '.trimal'
	
	print('trimal -in', args.fasta, '-out', outfile, '-gt 0.1')


def main():
	print_trimal()


# get the arguments before calling main
args = get_args()


# execute the program by calling main
if __name__ == "__main__":
	main()

