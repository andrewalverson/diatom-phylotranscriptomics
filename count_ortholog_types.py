#! /usr/bin/env python3

# load required modules
import argparse
import os,sys,re

inclades_with_subtrees = {}
total_inclades = {}

subtree_count  = 0
ortho_count    = 0
unrooted_count = 0
inclades_without_subtrees_count = 0
inclades_with_subtree_count    = 0

files = os.listdir()

# loop over *all* tree files
for f in files:
	# skip possible system files or non-tree files that might be in the directory
	if not (f.startswith("OG0")):
		continue
	# skip unrooted trees
	elif(re.search(".*unrooted.*", f)):
		unrooted_count += 1
	else:
		# match the "ortho*.tre" files, keeping the the *.inclade[0-9] part
		subtree_match = re.search("((^OG.*)\.inclade[1-9])(\.ortho[1-9]\.tre)", f)
		if subtree_match:
			# each subtree file is unique – count them here
			subtree_count += 1
			# print(subtree_match.group(0))
			# print(subtree_match.group(1))
			# store the base inclade prefix, i.e. store 'OG0008776_1_1.inclade1' for 'OG0008776_1_1.inclade1.ortho1.tre')
			# these are not unique because there can be multiple subtrees per inclade, so store in dictionary
			inclades_with_subtrees[subtree_match.group(1)] = 1
			# inclades.append(subtree_match.group(1))
			# print(inclade_match.group(1))
		# this is an inclade tree
		else:
			total_inclades[f] = 1


inclades_without_subtrees_count = len(total_inclades.keys()) - len(inclades_with_subtrees.keys())
ortho_count                     = subtree_count + unrooted_count


# print('total inclades without extracted subtrees):', inclade_with_subtree_count)
# print('unique orthogroups (subtrees + inclades with no extracted subtrees):', ortho_count + unrooted_count)

print('total inclades (above and below taxon occupancy threshold):', len(total_inclades.keys())) # verified correct
print('num trees below taxon occupancy threshold:', inclades_without_subtrees_count) # verified correct
print('num trees above taxon occupancy threshold:', subtree_count) # verified correct
print('num unrooted trees:', unrooted_count) # verified correct
print('total num ortholog trees above taxon occupancy threshold:', ortho_count)
# print('num inclades with subtrees :', len(inclades_with_subtrees.keys())) # verified correct


'''
# iterate again over all files, including inclades and subtrees
for i in files:
	# skip possible system files or non-tree files that might be in the directory
	if not (f.startswith("OG0")):
		continue
	# skip unrooted trees
	elif(re.search(".*unrooted.*", f)):
		ortho_count += 1
	# count the total number of unique inclades and subtrees
	# counting this way skips inclades that were stored in inclades{} **because**
	# they had a match to a subtree
	else:
		ortho_count += 1
		# print(i)
'''