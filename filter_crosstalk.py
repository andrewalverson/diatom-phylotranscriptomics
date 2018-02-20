#! /usr/bin/env python3

# this script is run in a directory that includes all of the assemblies that were multiplexed and run on the same Illumina lane
# it generates all the pairwise CD-HIT commands, which are designed to remove index crosstalk/barcode bleeding contigs from each sample

# about CD-HIT-EST-2D output files
# "the output are two files: a fasta file of sequences in db2 that are not similar to db1 and a text file that lists similar sequences between db1 & db2"

# load required modules
import sys
import os
import re
from Bio import SeqIO

min_seq_length = 300

# get current working directory
cwd = os.getcwd()
    
# read in all the contents of the directory
directory_listing = os.listdir(cwd)

# loop over all the files in this directory, considering only *.fa files
for file1 in directory_listing:
    from collections import defaultdict
    # delete_seqs = defaultdict(dict)
    delete_seqs = defaultdict(list)
    filename_re = '([RL]\d+)(_.*)\.fa$' # e.g., R25_ssu_diatom_rrna.fa
   
    if re.search( filename_re, file1 ):
        for file2 in directory_listing:
            if re.search( filename_re, file2 ) and ( file1 != file2 ):
                outfile_pt1 = re.search( filename_re, file1 ).group(1)
                outfile_pt2 = re.search( filename_re, file2 ).group(1)
                outfile = outfile_pt1 + '-' + outfile_pt2

                # generate and run CD-HIT-EST-2D command (nucleotide sequences)
                cd_hit = "cd-hit-est-2d" + " -c 1" + " -T 0" + " -d 300" + " -aS 1" + " -n 10" + " -p 1" + " -i " + file2 + " -i2 " + file1 + " -o " + outfile + "> /dev/null"
                os.system( cd_hit )
                os.remove( outfile )

                # open and parse *.clstr file
                clstr_file = open( outfile + ".clstr", 'r' )
                for line in clstr_file:
                    line = line.rstrip()

                    # find only lines that begin with 1-9, which are the smaller sub-matches
                    if re.search( '^([1-9]\d*)\s+', line ):
                        f = line.split()
                        delete_seqs[file1].append(f[2][1:-3])
                clstr_file.close()

    else:
        continue
    
    # set output filename
    filtered_fasta = re.search( filename_re, file1 ).group(1) + re.search( filename_re, file1 ).group(2) + '.filtered.fa'

    # open and parse file1 (FASTA file that we're filtering contaminants from)
    fasta_sequences = SeqIO.parse(open(file1),'fasta')

    # print(delete_seqs[file1])

    # open output file
    with open( filtered_fasta, "w" ) as f:
        for acc in fasta_sequences:
            if acc.id not in delete_seqs[file1] and len(acc.seq) > min_seq_length:
                SeqIO.write( [acc], f, "fasta" )


# the *.clstr files, which are formatted like so:

# >Cluster 0
# 0	1980nt, >TRINITY_DN1_c0_g1_i1_R15_Attheya_longicornis... *
# 1	480nt, >TRINITY_DN1_c0_g1_i1_R15_Attheya_longicornis_created... at 1:480:961:1440/+/100.00%

# Remember that I have specified the potential "contaminating" species as db1, so this species --
   # not the one we're looking for contaminants in -- is shown as the longer reference sequence in this output
# Thus, the clusters are defined based on sequences in the "contaminating" db1 file, not the target species
# So in this output, we want to remove '>TRINITY_DN1_c0_g1_i1_R15_Attheya_longicornis_created' from R12, which was specified
   # as db2 in this example



