library(here)
library(readxl)
library(tidyverse)
library(ggrepel)
library(viridis)
library(treespace)
library(ape)
library(ggtree)
library(treeio)
library(TreeDist)
library(TreeTools)
library(protoclust)

# load species trees metadata
metadata <- read_excel("species_trees/species_trees_metadata.xlsx")

# load trees
astral_pro_tax0.4                  <- read.tree('species_trees/tree_files/astral-pro-tax0.4.tre')
astral_pro_tax0.6                  <- read.tree('species_trees/tree_files/astral-pro-tax0.6.tre')
astral_pro_tax0.75                 <- read.tree('species_trees/tree_files/astral-pro-tax0.75.tre')
astral_pro_tax0.8                  <- read.tree('species_trees/tree_files/astral-pro-tax0.8.tre')
astral_pro_tax0.75_topPI           <- read.tree('species_trees/tree_files/astral-pro-tax0.75-top-PI.tre')
astral_pro_tax0.8_topPI            <- read.tree('species_trees/tree_files/astral-pro-tax0.8-top-PI.tre')

astral_weighted_0.4                <- read.tree('species_trees/tree_files/astral-weighted-0.4.tre')
astral_weighted_0.6                <- read.tree('species_trees/tree_files/astral-weighted-0.6.tre')
astral_weighted_0.7                <- read.tree('species_trees/tree_files/astral-weighted-0.7.tre')
astral_weighted_0.75               <- read.tree('species_trees/tree_files/astral-weighted-0.75.tre')
astral_weighted_0.8                <- read.tree('species_trees/tree_files/astral-weighted-0.8.tre')
astral_weighted_0.75_topPI         <- read.tree('species_trees/tree_files/astral-weighted-0.75-top-PI.tre')
astral_weighted_0.8_topPI          <- read.tree('species_trees/tree_files/astral-weighted-0.8-top-PI.tre')

astral_weighted_thal_raphid        <- read.tree('species_trees/tree_files/astral-weighted-thal-raphid.tre')
astral_weighted_thal_raphid_tax0.6 <- read.tree('species_trees/tree_files/astral-weighted-thal-raphid-tax0.6.tre')

concat_tax0.6_symtest              <- read.tree('species_trees/tree_files/concat-tax0.6-symtest.treefile')
concat_tax0.7_symtest              <- read.tree('species_trees/tree_files/concat-tax0.7-symtest.treefile')
concat_tax0.75_symtest             <- read.tree('species_trees/tree_files/concat-tax0.75-symtest.treefile')
concat_tax0.8_symtest              <- read.tree('species_trees/tree_files/concat-tax0.8-symtest.treefile')
concat_tax0.75_topPI_symtest       <- read.tree('species_trees/tree_files/concat-tax0.75-top-PI-symtest.treefile')
concat_tax0.8_topPI_symtest        <- read.tree('species_trees/tree_files/concat-tax0.8-top-PI-symtest.treefile')

pmsf_tax0.75_topPI                <- read.tree('species_trees/tree_files/tax-0.75-top-PI-pmsf.treefile')
pmsf_tax0.8_topPI                 <-  read.tree('species_trees/tree_files/tax-0.8-top-PI-pmsf.treefile')
pmsf_tax0.6_thal_raphid           <-  read.tree('species_trees/tree_files/thal-raphid-tax0.6-pmsf.treefile')


# root trees
astral_pro_tax0.4 <- drop.tip(astral_pro_tax0.4, c("Ectocarpussil","Nannochloropsisgad"))
astral_pro_tax0.4 <- root.phylo(phy = astral_pro_tax0.4,
                                outgroup = "Aureococcusanoph",
                                resolve.root = TRUE,
                                edgelabel = TRUE)
astral_pro_tax0.6 <- drop.tip(astral_pro_tax0.6, c("Ectocarpussil","Nannochloropsisgad"))
astral_pro_tax0.6 <- root.phylo(phy = astral_pro_tax0.6,
                                outgroup = "Aureococcusanoph",
                                resolve.root = TRUE,
                                edgelabel = TRUE)
astral_pro_tax0.75 <- drop.tip(astral_pro_tax0.75, c("Ectocarpussil","Nannochloropsisgad"))
astral_pro_tax0.75 <- root.phylo(phy = astral_pro_tax0.75,
                                 outgroup = "Aureococcusanoph",
                                 resolve.root = TRUE,
                                 edgelabel = TRUE)
astral_pro_tax0.75_topPI <- drop.tip(astral_pro_tax0.75_topPI, c("Ectocarpussil","Nannochloropsisgad"))
astral_pro_tax0.75_topPI <- root.phylo(phy = astral_pro_tax0.75_topPI,
                                       outgroup = "Aureococcusanoph",
                                       resolve.root = TRUE,
                                       edgelabel = TRUE)
astral_pro_tax0.8 <- drop.tip(astral_pro_tax0.8_topPI, c("Ectocarpussil","Nannochloropsisgad"))
astral_pro_tax0.8 <- root.phylo(phy = astral_pro_tax0.8,
                                outgroup = "Aureococcusanoph",
                                resolve.root = TRUE,
                                edgelabel = TRUE)
astral_pro_tax0.8_topPI <- drop.tip(astral_pro_tax0.8_topPI, c("Ectocarpussil","Nannochloropsisgad"))
astral_pro_tax0.8_topPI <- root.phylo(phy = astral_pro_tax0.8_topPI,
                                      outgroup = "Aureococcusanoph",
                                      resolve.root = TRUE,
                                      edgelabel = TRUE)
astral_weighted_0.4 <- root.phylo(phy = astral_weighted_0.4,
                                  outgroup = "Aureococcusanoph",
                                  resolve.root = TRUE,
                                  edgelabel = TRUE)
astral_weighted_0.6 <- root.phylo(phy = astral_weighted_0.6,
                                  outgroup = "Aureococcusanoph",
                                  resolve.root = TRUE,
                                  edgelabel = TRUE)
astral_weighted_0.7 <- root.phylo(phy = astral_weighted_0.7,
                                  outgroup = "Aureococcusanoph",
                                  resolve.root = TRUE,
                                  edgelabel = TRUE)
astral_weighted_0.75 <- root.phylo(phy = astral_weighted_0.75,
                                   outgroup = "Aureococcusanoph",
                                   resolve.root = TRUE,
                                   edgelabel = TRUE)
astral_weighted_0.8 <- root.phylo(phy = astral_weighted_0.8,
                                  outgroup = "Aureococcusanoph",
                                  resolve.root = TRUE,
                                  edgelabel = TRUE)
astral_weighted_0.75_topPI <- root.phylo(phy = astral_weighted_0.75_topPI,
                                         outgroup = "Aureococcusanoph",
                                         resolve.root = TRUE,
                                         edgelabel = TRUE)
astral_weighted_0.8_topPI <- root.phylo(phy = astral_weighted_0.8_topPI,
                                        outgroup = "Aureococcusanoph",
                                        resolve.root = TRUE,
                                        edgelabel = TRUE)
astral_weighted_thal_raphid <- root.phylo(phy = astral_weighted_thal_raphid,
                                          outgroup = "Aureococcusanoph",
                                          resolve.root = TRUE,
                                          edgelabel = TRUE)
astral_weighted_thal_raphid_tax0.6 <- root.phylo(phy = astral_weighted_thal_raphid_tax0.6,
                                        outgroup = "Aureococcusanoph",
                                        resolve.root = TRUE,
                                        edgelabel = TRUE)
concat_tax0.6_symtest <- root.phylo(phy = concat_tax0.6_symtest,
                                    outgroup = "Aureococcusanoph",
                                    resolve.root = TRUE,
                                    edgelabel = TRUE)
concat_tax0.7_symtest <- root.phylo(phy = concat_tax0.7_symtest,
                                    outgroup = "Aureococcusanoph",
                                    resolve.root = TRUE,
                                    edgelabel = TRUE)
concat_tax0.75_symtest <- root.phylo(phy = concat_tax0.75_symtest,
                                     outgroup = "Aureococcusanoph",
                                     resolve.root = TRUE,
                                     edgelabel = TRUE)
concat_tax0.8_symtest <- root.phylo(phy = concat_tax0.8_symtest,
                                    outgroup = "Aureococcusanoph",
                                    resolve.root = TRUE,
                                    edgelabel = TRUE)
concat_tax0.75_topPI_symtest <- root.phylo(phy = concat_tax0.75_topPI_symtest,
                                           outgroup = "Aureococcusanoph",
                                           resolve.root = TRUE,
                                           edgelabel = TRUE)
concat_tax0.8_topPI_symtest <- root.phylo(phy = concat_tax0.8_topPI_symtest,
                                          outgroup = "Aureococcusanoph",
                                          resolve.root = TRUE,
                                          edgelabel = TRUE)
pmsf_tax0.75_topPI <- root.phylo(phy = pmsf_tax0.75_topPI,
                                  outgroup = "Aureococcusanoph",
                                  resolve.root = TRUE,
                                  edgelabel = TRUE)
pmsf_tax0.8_topPI <- root.phylo(phy = pmsf_tax0.8_topPI,
                                  outgroup = "Aureococcusanoph",
                                  resolve.root = TRUE,
                                  edgelabel = TRUE)
pmsf_tax0.6_thal_raphid <- root.phylo(phy = pmsf_tax0.6_thal_raphid,
                                                         outgroup = "Aureococcusanoph",
                                                         resolve.root = TRUE,
                                                         edgelabel = TRUE)


# collect all trees
all_trees <- c(astral_pro_tax0.4, astral_pro_tax0.6, astral_pro_tax0.75, 
               astral_pro_tax0.8, astral_pro_tax0.75_topPI, astral_pro_tax0.8_topPI,
               astral_weighted_0.4, astral_weighted_0.6, astral_weighted_0.7, 
               astral_weighted_0.75, astral_weighted_0.8, astral_weighted_0.75_topPI,
               astral_weighted_0.8_topPI, astral_weighted_thal_raphid,
               astral_weighted_thal_raphid_tax0.6,
               concat_tax0.6_symtest, concat_tax0.7_symtest, concat_tax0.75_symtest, 
               concat_tax0.8_symtest, concat_tax0.75_topPI_symtest, concat_tax0.8_topPI_symtest, 
               pmsf_tax0.75_topPI, pmsf_tax0.8_topPI, pmsf_tax0.6_thal_raphid)


# collect all astral trees
all_astral <- c(astral_pro_tax0.4, astral_pro_tax0.6, astral_pro_tax0.75, 
                astral_pro_tax0.8, astral_pro_tax0.75_topPI, astral_pro_tax0.8_topPI,
                astral_weighted_0.4, astral_weighted_0.6, astral_weighted_0.7, 
                astral_weighted_0.75, astral_weighted_0.8, astral_weighted_0.75_topPI,
                astral_weighted_0.8_topPI, astral_weighted_thal_raphid,
                astral_weighted_thal_raphid_tax0.6)

# collect all concat trees
all_concat <- c(concat_tax0.6_symtest, concat_tax0.7_symtest, concat_tax0.75_symtest, 
                concat_tax0.8_symtest, concat_tax0.75_topPI_symtest, concat_tax0.8_topPI_symtest, 
                pmsf_tax0.75_topPI, pmsf_tax0.8_topPI, pmsf_tax0.6_thal_raphid)


# create multiPhylo objects
class(all_trees) <- 'multiPhylo'
class(all_astral) <- 'multiPhylo'
class(all_concat) <- 'multiPhylo'

# create tree names for the plot
names(all_trees)[1] <- paste0('Astral-Pro-Tax0.4')
names(all_trees)[2] <- paste0('Astral-Pro-Tax0.6')
names(all_trees)[3] <- paste0('Astral-Pro-Tax0.75')

names(all_trees)[4] <- paste0('Astral-Pro-Tax0.8')
names(all_trees)[5] <- paste0('Astral-Pro-Tax0.75-TopPI')
names(all_trees)[6] <- paste0('Astral-Pro-Tax0.8-TopPI')

names(all_trees)[7] <- paste0('Astral-Weighted-Tax0.4')
names(all_trees)[8] <- paste0('Astral-Weighted-Tax0.6')
names(all_trees)[9] <- paste0('Astral-Weighted-Tax0.7')

names(all_trees)[10] <- paste0('Astral-Weighted-Tax0.75')
names(all_trees)[11] <- paste0('Astral-Weighted-Tax0.8')
names(all_trees)[12] <- paste0('Astral-Weighted-Tax0.75-TopPI')
names(all_trees)[13] <- paste0('Astral-Weighted-Tax0.8-TopPI')

names(all_trees)[14] <- paste0('Astral-Weighted-Thal-Raphid')
names(all_trees)[15] <- paste0('Astral-Weighted-Thal-Raphid-Tax0.6')

names(all_trees)[16] <- paste0('Concat-Tax0.6')
names(all_trees)[17] <- paste0('Concat-Tax0.7')
names(all_trees)[18] <- paste0('Concat-Tax0.75')

names(all_trees)[19] <- paste0('Concat-Tax0.8')
names(all_trees)[20] <- paste0('Concat-Tax0.75-TopPI')
names(all_trees)[21] <- paste0('Concat-Tax0.8-TopPI')

names(all_trees)[22] <- paste0('PMSF-TopPI-Tax0.75')
names(all_trees)[23] <- paste0('PMSF-TopPI-Tax0.8')
names(all_trees)[24] <- paste0('PMSF-Thal-Raphid-Tax0.6')


# create vector of the tree types
analysis <- c(rep('ASTRAL-PRO',6), rep('ASTRAL-WEIGHTED',9), rep('CONCAT', 6), rep('PMSF', 3))

# create list of taxon occupancy
tax_occ <- c('0.4', '0.6', '0.75',
             '0.8', '0.75', '0.8',
             '0.4', '0.6', '0.7',
             '0.75', '0.8', '0.75', '0.8',
             '0.4', '0.6',
             '0.6', '0.7', '0.75',
             '0.8', '0.75', '0.8',
             '0.75', '0.8', '0.6')

a_name = c('Astral-Pro-Tax0.4', 'Astral-Pro-Tax0.6', 'Astral-Pro-Tax0.75',
           'Astral-Pro-Tax0.8', 'Astral-Pro-tax0.75-TopPI', 'Astral-Pro-Tax0.8-TopPI',
           'Astral-Weighted-Tax0.4', 'Astral-Weighted-Tax0.6', 'Astral-Weighted-Tax0.7',
           'Astral-Weighted-Tax0.75', 'Astral-Weighted-Tax0.8',
           'Astral-Weighted-Tax0.75-TopPI', 'Astral-Weighted-Tax0.8-TopPI', 
           'Astral-Weighted-Thal-Raphid', 'Astral-Weighted-Thal-Raphid-Tax0.6',
           'Concat-Tax0.6', 'Concat-Tax0.7', 'Concat-Tax0.75', 
           'Concat-Tax0.8', 'Concat-Tax0.75-TopPI', 'Concat-Tax0.8-TopPI', 
           'PMSF-TopPI-Tax0.75', 'PMSF-TopPI-Tax0.8', 'PMSF-Thal-Raphid-Tax0.6')

# use TreeDist to calculate a more accurate distance measure
distances <- ClusteringInfoDistance(all_trees)
mapping   <- cmdscale(distances, k = 12)
mapping   <- data.frame(mapping)

# add metadata columns to mapping data frame
mapping <- mapping %>% 
  mutate(occupancy = metadata$min.taxon.occupancy) %>% 
  mutate(tree.type = metadata$method) %>% 
  mutate(name = metadata$name) %>% 
  mutate(num.orthologs = metadata$num.orthologs)

# mapping <- mapping %>% 
#   mutate(occ = tax_occ) %>% 
#   mutate(atype = analysis) %>% 
#   mutate(label = a_name)

# manually set breaks
bb = c(100, 500, 1500, 2500)

dist <- ggplot(mapping, aes(x = X1, y = X2)) +
  geom_point(aes(color = occupancy, shape = tree.type, size = num.orthologs)) +
  scale_size_continuous(name = "Num orthologs", 
                       range = c(2, 12),
                       breaks = bb,
                       limits = c(100, 3000)) +
  scale_color_viridis() +
  geom_label_repel(aes(label = name),
                  size = 2,
                  box.padding = 0.35, 
                  point.padding = 0.5,
                  segment.color = 'grey50') +
  labs(color="Minimum occupancy",shape="Analysis") +
  theme_bw() +
  theme(axis.title = element_blank(), text = element_text(size = 14)) + 
  guides(shape = guide_legend(override.aes = list(size = 3)))

ggsave(dist, filename = "figs/trees-ordination.png",
       height = 5, width = 8)

ggsave(dist, filename = "figs/trees-ordination.pdf",
       height = 5, width = 8)

############################## TREE CLUSTERING ##############################

# conclusion: meh support for 10 clusters

# use TreeDist to cluster trees
possibleClusters <- 2:10

pamClusters <- lapply(possibleClusters, function(k) cluster::pam(distances, k = k))
pamSils <- vapply(pamClusters, function(pamCluster) {
  mean(cluster::silhouette(pamCluster)[, 3])
}, double(1))

bestPam <- which.max(pamSils)
pamSil <- pamSils[bestPam]
pamCluster <- pamClusters[[bestPam]]$cluster

hTree <- protoclust(distances)
hClusters <- lapply(possibleClusters, function(k) cutree(hTree, k = k))
hSils <- vapply(hClusters, function(hCluster) {
  mean(cluster::silhouette(hCluster, distances)[, 3])
}, double(1))


bestH <- which.max(hSils)
hSil <- hSils[bestH]
hCluster <- hClusters[[bestH]]

plot(pamSils ~ possibleClusters,
     xlab = "Number of clusters", ylab = "Silhouette coefficient",
     ylim = range(c(pamSils, hSils)))
points(hSils ~ possibleClusters, pch = 2)
legend("topleft", c("PAM", "Hierarchical"), pch = 1:2)


cluster <- hClusters[[10 - 1]]

class(hTree) <- "hclust"
par(mar = c(0, 0, 0, 0))
plot(hTree, labels = FALSE, main = "")
points(seq_along(trees), rep(1, length(trees)), pch = 16,
       col = spectrum[hTree$order])

par(mar = rep(0, 4))
plot(mapping,
     asp = 1, # Preserve aspect ratio - do not distort distances
     ann = FALSE, axes = FALSE, # Don't label axes: dimensions are meaningless
     pch = 16
)


# Robinson-Foulds distances
# Dscape <- treespace(all_trees[c(1:19)], method='RF', nf=2)
Dscape <- treespace(all_trees, method='RF', nf=2)
plotGrovesD3(Dscape$pco, groups=analysis[c(1:23)], treeNames=names(all_trees), symbol_var = tax_occ[c(1:23)],
             colors=c('#E69F00','#56B4E9','#009E73', '#CC79A7'))


############################## CONSENSUS TREES ##############################

# compute and write out consensus trees
consensus_astral <- consensus(all_astral, p=0.7)
consensus_concat <- consensus(all_concat, p=0.7)
consensus_all    <- consensus(all_trees, p=0.7)
consensus_all_strict    <- consensus(all_trees, p=1)

write.tree(consensus_astral, file = 'consensus-astral-0.7.tre')
write.tree(consensus_concat, file = 'consensus-concat-0.7.tre')
write.tree(consensus_all,    file = 'consensus-all-0.7.tre')
write.tree(consensus_all_strict,    file = 'consensus-all-strict.tre')

# read in metadata for tip labels
metadata <- read_excel("voucher-list.xlsx", sheet = "292-dataset")

# make columns with formatted labels
metadata <- metadata %>% 
  mutate(formatted.label = paste(metadata$final.genus,
                                 metadata$final.species,
                                 metadata$strain,
                                 sep = "_"))

# define "deep" outgroup clade
out.vec = c('Ectocarpussil', 'Nannochloropsisgad')

#### run the following code for each of the consensus trees ####
# tree <- consensus_astral
# tree <- consensus_concat
# tree <- consensus_all
tree <- consensus_all_strict

# drop distant outgroups
tree <- drop.tip(tree, out.vec)

# root the tree
rooted.tree <- root.phylo(phy = tree,
                          outgroup = "Aureococcusanoph",
                          resolve.root = TRUE,
                          edgelabel = TRUE)

# get vector of tips to drop
drop.dat <- metadata %>% filter(keep == "drop")
drop.vec <- drop.dat$label

# drop tips
rooted.tree <- drop.tip(rooted.tree, tip=drop.vec)

metadata_tips_only <- metadata[metadata$label %in% rooted.tree$tip.label,]

### group by morphological "group" == radial, polar, araphid, or raphid
# creates a list of lists, each one is the tips that belong to each group
major_group_info <- split(metadata_tips_only$label, metadata_tips_only$group)
major_groups <- groupOTU(rooted.tree, major_group_info)

# rename tips
major_groups.tree <- rename_taxa(major_groups, metadata_tips_only, key = label, value = formatted.label)

major_group_tree_fig <- ggtree(major_groups.tree, size = 0.4, aes(color=group), branch.length = "none") +
  scale_y_reverse() +
  ggplot2::xlim(0,40) +
  geom_tiplab(size = 2, family = "Helvetica") +
  # scale_x_continuous(limits=c(0,40)) +
  # geom_treescale(fontsize=4, linesize=0.2, x=1.5, y=-100) +
  theme_tree(legend.position = c(0.1, 0.1)) +
#  ggtitle("Astral consensus") +
#  ggtitle("Concatenated consensus") +
#  ggtitle("All consensus (70%)") +
  ggtitle("All consensus (strict)") +
  theme(plot.title = element_text(hjust = 0.5))

# ggsave(major_group_tree_fig, file='species_trees/consensus-astral.pdf', width = 8, height = 20, device = "pdf")
# ggsave(major_group_tree_fig, file='species_trees/consensus-concat.pdf', width = 8, height = 20, device = "pdf")
# ggsave(major_group_tree_fig, file='species_trees/consensus-all.pdf', width = 8, height = 20, device = "pdf")
# ggsave(major_group_tree_fig, file='species_trees/consensus-all-strict.pdf', width = 8, height = 20, device = "pdf")


# plot(consensus_astral, cex=0.5, node.depth=2)
# edgelabels(cex=0.5)
# consensus_astral <- castor::root_in_edge(consensus_astral, root_edge=1)
# consensus_astral <- ladderize(consensus_astral)
# write.tree(consensus_astral, file='consensus-50-all-astral.tree')
# 
# pdf(file='consensus-50-all-astral.pdf', height=8, width=6)
# plot(consensus_astral, cex=0.5, node.depth=2, no.margin=T)
# dev.off()
