#! /share/apps/Python/3.5.2/bin/python3.5

import sys
import os
import re

# get directory contents
folder = os.listdir(".")

# loop over files
for fi in folder:
    # only need to consider one of the two files per sample, so select "contam"
    if re.search("contam\.fasta", fi):
		# split filename on the underscore
    	parts = fi.split("_")
    	cat_command = "cat" + " " + fi + " " + parts[0] + "_dubious.fasta" + " > " + parts[0] + "_croco_removed.fasta"
    	# print(cat_command)
    	os.system(cat_command)

    else:
    	continue
