#! /usr/bin/env python3

# load required modules
import argparse
import re
import statistics
from Bio import SeqIO

# create an ArgumentParser object ('parser') that will hold all the information necessary to parse the command line
parser = argparse.ArgumentParser(description="This script makes a global list of masked sites (from UPP) and trimmed sites (from trimal) to make a global list of sites to trim from an alignment")

# add arguments
parser.add_argument( "alignment", help="original FASTA alignment" )

args = parser.parse_args()

# parse name of alignment file in order to set prefix for output file names
a = args.alignment.split("_")
prefix = a[0]

# set paths to input and output files
alignment_path = '/storage/aja/orthofinder/orthogroups/orthogroup_alignments_occupancy_10/AA_alignments_correct_headers/'
trimal_path    = '/storage/aja/orthofinder/orthogroups/orthogroup_alignments_occupancy_10/AA_alignments_trimal/'
out_path       = '/storage/aja/orthofinder/orthogroups/orthogroup_alignments_occupancy_10/AA_alignments_trimmed/'

# name of trimal column mapping file
trimal = prefix + '_trimal.sites'

# name of UPP masked sites file
upp_masked = prefix + '_insertion_columns.txt'

# name of trimal error file
err = prefix + '_alignment.err'

# set walltime and number of cpu's
wall    = '06'
num_cpu = '1'


alignment_length = 0
trimmed_columns = []

# get alignment length
fasta_sequences  = list(SeqIO.parse(open(args.alignment),'fasta'))
alignment_length = len(fasta_sequences[0].seq)
# print(alignment_length)

# parse trimal mapping
with open(trimal, 'r') as trimal_mapping:
	untrimmed_columns = trimal_mapping.read()
	untrimmed_columns = untrimmed_columns.rstrip()

# remove trimal line header
untrimmed_columns = re.sub(r'#ColumnsMap\s+', r'', untrimmed_columns)
# print(untrimmed_columns)

# covert untrimmed_columns from string to list
untrimmed_columns_list = untrimmed_columns.split(', ')
# for i in untrimmed_columns_list:
# 	print(i)

# from the list of *retained* columns from trimal, get the inverse list of *trimmed* columns
for i in range(0, alignment_length):
	if str(i) in untrimmed_columns_list:
		continue
	else:
		trimmed_columns.append(str(i))

# parse sites masked by UPP
with open(upp_masked, 'r') as masked:
	masked_columns = masked.read()
	masked_columns = masked_columns.rstrip()



# print top of PBS script
print('#PBS -N', prefix+'.trimal-AA')
print('#PBS -q q06h32c')
print('#PBS -j oe')
print('#PBS -o',  prefix+'.trimal-AA.$PBS_JOBID')
print('#PBS -l nodes=1:ppn=' + num_cpu)
print('#PBS -l walltime=' + wall + ':00:00')
print()

print('cd /scratch/$PBS_JOBID')
print()

print('# run trimal')
print('/share/apps/bioinformatics/trimal/1.4rev22/trimal', '-selectcols { ' + masked_columns + ',' + ','.join(trimmed_columns) + ' }', '-keepheader', '-colnumbering', '-in', alignment_path + args.alignment, '-out', out_path + args.alignment, '>', out_path + trimal, '2>', out_path + err)


