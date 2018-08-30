#! /usr/bin/env python3

# this script parses the output of `summarize_orthogroup_membership.py` to select orthogroups with taxon occupancy >= a specified
# threshold, then makes symbolic links to those files

# load required modules
import sys
import os
import re
import argparse

# create an ArgumentParser object ('parser') that will hold all the information necessary to parse the command line
parser = argparse.ArgumentParser(description="Parse output of `summarize_orthogroup_membership.py` to create symbolic links to FASTA files of orthogroups with a minimum specified taxon occupancy")

# add arguments
parser.add_argument( "orthogroups_summary", help="output of summarize_orthogroup_membership.py" )
parser.add_argument( "directory", help="Path to orthogroup FASTA files" )
parser.add_argument( "--min_diatom_occ", help="Link to orthogroups with diatom taxon occupancy >= proportion <min_diatom_occ>", default=0.1, type=float )

args = parser.parse_args()

# check user input
# exit if the directory of FASTA files doesn't exist
os.path.isdir(args.directory) or sys.exit("The file", args.directory, "does not exist – exiting.")

# check whether the orthogroup summary file exists
os.path.isfile(args.orthogroups_summary)  or sys.exit("The file", args.orthogroups_summary, "does not exist – exiting.")

# check for trailing '/' on the path to the orthogroup FASTA files – add one if it's missing
if not args.directory.endswith('/'):
    args.directory += '/'

# open and parse orthogroup summary file
with open(args.orthogroups_summary, 'r') as orthogroups:
    # skip first line
    next(orthogroups)
    
    for line in orthogroups:
        line = line.rstrip()

        a = line.split()

        if float(a[4]) >= 0.1:
            link = "ln -s " + args.directory + a[0] + ".fa" + "\n"
            #print(link)
            os.system(link)
        else:
            continue
