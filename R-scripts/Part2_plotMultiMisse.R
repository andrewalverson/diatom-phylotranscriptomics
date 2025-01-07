
######################################################################################################################################
######################################################################################################################################
### Plots a densitree of a set of MiSSE runs:
######################################################################################################################################
######################################################################################################################################

source("Densi.mod.R")
library(hisse)

#Step 1: Get node and tip rates. For each tree, create edge rate vector, loop through and take average of rootward and tipward rates,
GetMapEdgeTreeList <- function(prefix="diatom_258*", rate.type="turnover"){
	runs <- system(paste("ls -1 ../FinalMiSSERunsFIXED/",  prefix, sep=""), intern=TRUE)
	tree.list <- as.list(1:length(runs))
	for(tree.index in 1:length(runs)){
		load(runs[tree.index])
		tmp.tree <- model.recons.misse[[1]]$phy
		tmp.rate <- c()
		rate.mat <- GetModelAveRates(model.recons.misse, type="both")
		tip.mat.modified <- rate.mat$tips
		tip.mat.modified <- tip.mat.modified[,-1]
		id <- 1:dim(tip.mat.modified)[1]
		tip.mat.modified <- cbind(id, tip.mat.modified)
		rate.total <- rbind(tip.mat.modified, rate.mat$nodes)
		for(row.index in 1:dim(tmp.tree$edge)[1]){
			tmp.rate <- c(tmp.rate, mean(c(rate.total[tmp.tree$edge[row.index,1], rate.type], rate.total[tmp.tree$edge[row.index,2], rate.type])))
		}
		tmp.tree$edge.rate <- tmp.rate
		tree.list[[tree.index]] <- tmp.tree
	}
	class(tree.list) <- "multiPhylo"
	return(tree.list)
}



#Step 2: Now pass this to Densi.plot function, and set the colors:
MakePlot <- function(trees, cex=0.01, rate.colors=NULL, consensus=NULL){
	if(is.null(rate.colors)){
		rate.colors <- c("blue", "red")
	}
	rate.colors <- colorRampPalette(rate.colors, space="Lab")(1001)
	
	for(tree.index in 1:length(trees)){
		rates.to.plot <- trees[[tree.index]]$edge.rate
		rate.lims <- range(rates.to.plot)
		lims.percentage.correction <- 0.001
		rate.lims[1] <- rate.lims[1] - lims.percentage.correction*abs(rate.lims[1])
		rate.lims[2] <- rate.lims[2] + lims.percentage.correction*abs(rate.lims[2])
		rate_normalized <- round((rate.lims[2] - rates.to.plot)/(rate.lims[2] - rate.lims[1]), 3)*1000
		rate.cols <- rev(colorRampPalette(rate.colors, space="Lab")(1001))
		trees[[tree.index]]$edge.colors <- rate.cols[rate_normalized]
	}
	densiTree(trees, cex=cex, type="cladogram", width=3, consensus=consensus)
}

turn.trees <- GetMapEdgeTreeList(rate.type="turnover")
netdiv.trees <- GetMapEdgeTreeList(rate.type="net.div")

trees.con <- read.tree("../InputTrees/thal-raphid-empirical-brlens.rooted.pruned.dated.tre")
trees.order <- ladderize(trees.con)

pdf("diatom_MultiMisseTURN.pdf", height=8, width=6)
MakePlot(turn.trees, cex=.2, consensus=trees.order, rate.colors=c("#3E049C", "#FCD225"))
dev.off()

pdf("diatom_MultiMisseNETDIV.pdf", height=8, width=6)
MakePlot(netdiv.trees, cex=.2, consensus=trees.order, rate.colors=c("#3E049C", "#FCD225"))
dev.off()



