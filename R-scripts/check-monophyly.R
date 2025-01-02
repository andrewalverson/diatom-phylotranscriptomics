library(ape)
library(castor)
library(readxl)
library(tidyverse)
library(here)

# get taxonomy metadata
metadata <- read_excel("voucher-list.xlsx")

# get tip labels for Thalassiosirales
thals <- metadata %>% filter(order == 'Thalassiosirales') %>% 
  select(label)

# get tip labels for raphid pennates
raphids <- metadata %>% filter(group == 'raphid') %>% 
  select(label)

# get tip labels for pennates (araphid + raphid)
# note that some araphids are probably misclassified based on their GenBank taxonomy,
# so this might not be a good monophyly test
pennates <- metadata %>% filter(class == 'Bacillariophyceae') %>% 
  select(label)

# specify the directory of gene trees and the file endings
trees_dir <- file.path('RT_ortholog_trees_final')
trees_gt <- dir(path=trees_dir, pattern='.treefile')

# function to check monophyly
check_monophyly <- function(file) {
  # this line removes the 'treefile' suffix from the tree filename
  OG_no <- sub('.treefile$', '\\1', perl=TRUE, x=file)

  # read in tree
  tree <- read.tree(paste(trees_dir, file, sep='/'))

  # create list of tips in gene trees
  list1 <- thals$label[thals$label %in% tree$tip.label]
  list2 <- raphids$label[raphids$label %in% tree$tip.label]
  list3 <- pennates$label[pennates$label %in% tree$tip.label]
  
  # check if Thalassiosirales are monophyletic
  thals_monophyly <- is_monophyletic(tree=tree, focal_tips=list1)
  # check if raphids are monophyletic
  raphid_monophyly <- is_monophyletic(tree=tree, focal_tips=list2)
  # check if pennates are monophyletic
  pennate_monophyly <- is_monophyletic(tree=tree, focal_tips=list3)
  
  return(c(file, thals_monophyly, raphid_monophyly, pennate_monophyly))
}

# loop over all gene trees and return a data frame with results
monophyly <- lapply(trees_gt, check_monophyly)
monophyly <- data.frame(matrix(unlist(monophyly), 
                                       nrow=(length(monophyly)),
                                       byrow=T))
colnames(monophyly) <- c('OG', 'Thalassiosirales', 'Raphids', 'Pennates')

# count number of gene trees where Thalassiosirales is monophyletic
monophyly %>% 
  count(Thalassiosirales)

# count number of gene trees where raphids are monophyletic
monophyly %>% 
  count(Raphids)

# count number of gene trees where raphids are monophyletic
monophyly %>% 
  count(Pennates)

# get gene trees with monophyletic Thalassiosirales
thal_trees <- monophyly %>% 
  filter(Thalassiosirales == TRUE)

# get gene trees with monophyletic Thalassiosirales *and* monophyletic raphids
thal_raphid_trees <- thal_trees %>%  
  filter(Raphids == TRUE)

# export names of trees with monophyletic raphids and Thalassiosirales
write_csv(tibble(thal_raphid_trees$OG),
          file = 'ortholog_sets_for_species_trees/thal_raphid_trees.csv', col_names = F)

# compare thal-raphid-monophyly trees against taxon-occupancy cutoffs
tax0.4 <- read.csv(file = "ortholog_sets_for_species_trees/tax0.4_trees.csv", header = F)
names(tax0.4) <- c('OG')

tax0.6 <- read.csv(file = "ortholog_sets_for_species_trees/tax0.6_trees.csv", header = F)
names(tax0.6) <- c('OG')

tax0.7 <- read.csv(file = "ortholog_sets_for_species_trees/tax0.7_trees.csv", header = F)
names(tax0.7) <- c('OG')

tax0.75 <- read.csv(file = "ortholog_sets_for_species_trees/tax0.75_trees.csv", header = F)
names(tax0.75) <- c('OG')

tax0.8  <- read.csv(file = "ortholog_sets_for_species_trees/tax0.8_trees.csv", header = F)
names(tax0.8) <- c('OG')

# dplyr's intersect() and setdiff() commands work with dataframes, so
# make a dataframe of the thal-raphid tree names, then change the column header
# back to "OG"
trt <- as.data.frame(thal_raphid_trees$OG)
trt <- rename(trt, 'OG' = 'thal_raphid_trees$OG')

# find intersection of thal-raphid trees and ones of different taxon occupancy
# substitute in different values
intersect(trt, tax0.4)
setdiff(trt, tax0.4)

# write out the full set of gene trees with monophyly of thals and raphids
write_csv(as.data.frame(thal_raphid_trees$OG),
          file = 'ortholog_sets_for_species_trees/thal_raphid_tax0.4_trees.csv', col_names = F)

a <- intersect(trt, tax0.6)

# write out gene trees with monophyly of thals and raphids *and* taxon occupancy >= 0.6
write_csv(a,
          file = 'ortholog_sets_for_species_trees/thal_raphid_tax0.6_trees.csv', col_names = F)

