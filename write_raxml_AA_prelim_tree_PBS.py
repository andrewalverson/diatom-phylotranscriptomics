#! /usr/bin/env python3

import sys
import os
import argparse
import random

# create an ArgumentParser object ('parser') that will hold all the information necessary to parse the command line
parser = argparse.ArgumentParser()

# add arguments
parser.add_argument( "alignment", help="FASTA file" )
parser.add_argument( "-q", "--queue", help="name of Razor queue", default="mem512GB64core" )

args = parser.parse_args()

# set walltime and number of cpu's
wall    = 72
num_cpu = 1

# parse name of input file in order to set prefix for output file names
a = args.alignment.split(".")
base_filename = a[0]

# print top of PBS script
print('#PBS -N', base_filename+'.raxml')
print('#PBS -q', args.queue)
print('#PBS -j oe')
print('#PBS -o',  base_filename+'.raxml.$PBS_JOBID')
print('#PBS -l nodes=1:ppn=' + str(num_cpu))
print('#PBS -l walltime=' + str(wall) + ':00:00')
print()

print('cd $PBS_O_WORKDIR')
print('cp', args.alignment, '/scratch/$PBS_JOBID')
print('cd /scratch/$PBS_JOBID')
print()

print('# load modules')
print('module load intel/12.0.0 openmpi/1.5.1 raxml/8.2.11')
print()

print('# run RAxML')
print('raxmlHPC-SSE3', '-s', args.alignment, '-f d', '-m PROTCATWAG', '-n', base_filename, '-p', str(random.randint(100,100000)))
# print('raxmlHPC-AVX-PTHREADS', '-s', args.alignment, '-f d', '-m PROTCATWAG', '-n', base_filename, '-T', num_cpu, '-p', str(random.randint(100,100000)))
print('cp /scratch/$PBS_JOBID/RAxML_bestTree.' + base_filename, '$PBS_O_WORKDIR/' + base_filename + '_raxml.tre' )
