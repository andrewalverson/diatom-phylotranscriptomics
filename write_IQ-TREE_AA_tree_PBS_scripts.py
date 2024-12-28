#! /usr/bin/env python3

import sys
import os
import re
import argparse
import statistics
from Bio import SeqIO

# create an ArgumentParser object ('parser') that will hold all the information necessary to parse the command line
parser = argparse.ArgumentParser()

# add arguments
parser.add_argument( "fasta", help="FASTA alignment" )

args = parser.parse_args()

# set walltime and number of cpu's
wall    = '168'
num_cpu = '8'

# set input directory
indir  = '/storage/aja/orthofinder/orthogroups/gene_trees_round2/aa_iqtree'

# set output directory
outdir = '/storage/aja/orthofinder/orthogroups/gene_trees_round2/aa_iqtree'

# parse name of input file in order to set prefix for output file names
a = args.fasta.split(".")
prefix = a[0]

# print top of PBS script
print('#PBS -N', prefix + '.IQTREE-AA')
print('#PBS -q onenode16core')
print('#PBS -j oe')
print('#PBS -o',  prefix + '.IQTREE-AA.$PBS_JOBID')
print('#PBS -l nodes=1:ppn=' + num_cpu)
print('#PBS -l walltime=' + wall + ':00:00')
print()

print('cd $PBS_O_WORKDIR')
print()

print('module load intel/18.0.1 impi/18.0.1 mkl/18.0.1')
print()

print('# run IQ-TREE')
print('/home/aja/local/bin/iqtree -s', '/'.join([indir, args.fasta]), '-pre', prefix, '-st AA -nt AUTO -bb 1000 -alrt 1000 -m TEST --runs 5 >', '/'.join([indir, prefix + '.err']))

# print('/share/apps/bioinformatics/iqtree/iqtree-1.6.9-Linux/bin/iqtree', '-m TEST --runs 100 -nt AUTO -bb 1000 -alrt 1000 -st AA', '-s', '/'.join([indir, args.fasta]), '-pre', prefix, '>', prefix + '.err')

