
######################################################################################################################################
######################################################################################################################################
### densiTree -- Modified Version from phangorn's Densitree function
######################################################################################################################################
######################################################################################################################################

library(ape)

getAges <- function(x) {
  fun <- function(x) max(node.depth.edgelength(x))
  height <- NULL
  if (inherits(x, "phylo")) height <- fun(x)
  if (inherits(x, "multiPhylo")) {
	if (!is.null(attr(x, "TipLabel"))) {
	  x <- .uncompressTipLabel(x)
	  x <- unclass(x)
	  height <- vapply(x, fun, 0)
	}
	else {
	  x <- unclass(x)
	  height <- vapply(x, fun, 0)
	}
  }
  height
}


add_tiplabels <- function(xy, tip.label, direction, adj, font, srt = 0, cex = 1, col = 1, label.offset = 0) {
  direction <- match.arg(direction, c("rightwards", "leftwards",  "upwards",
	"downwards"))
  horizontal <- direction %in% c("rightwards", "leftwards")
  nTips <- length(tip.label)
  xx <- rep(1, nrow(xy))
  yy <- xy[, 2 ]
  if (direction == "leftwards" | direction == "downwards") xx <- xx * 0
  if (!horizontal) {
#    tmp <- yy
	yy <- xx
	xx <- xy[, 1]
  }
  MAXSTRING <- max(strwidth(tip.label, cex = cex))
  loy <- 0
  if (direction == "rightwards") lox <- label.offset + MAXSTRING * 1.05 * adj
  if (direction == "leftwards")
	lox <- -label.offset - MAXSTRING * 1.05 * (1 - adj)
  if (!horizontal) {
	psr <- par("usr")
	MAXSTRING <- MAXSTRING * 1.09 * (psr[4] - psr[3]) / (psr[2] - psr[1])
	loy <- label.offset + MAXSTRING * 1.05 * adj
	lox <- 0
	srt <- 90 + srt
	if (direction == "downwards") {
	  loy <- -loy
	  srt <- 180 + srt
	}
  }
  text(xx[1:nTips] + lox, yy[1:nTips] + loy, tip.label, adj = adj,
	font = font, srt = srt, cex = cex, col = col)
}



densiTree <- function(x, type = "cladogram", alpha = 1 / length(x), consensus = NULL, direction = "rightwards", optim = FALSE, scaleX = FALSE, col = 1, width = 1, lty = 1, cex = .8, font = 3, tip.color = 1, adj = 0, srt = 0, underscore = FALSE, label.offset = 0, scale.bar = TRUE, jitter = list(amount = 0, random = TRUE), ...) {
  
  if (!inherits(x, "multiPhylo")) stop("x must be of class multiPhylo")
  if (is.character(consensus)) {
	consensus <- stree(length(consensus), tip.label = consensus)
	consensus$edge.length <- rep(1.0, nrow(consensus$edge))
  }
  if (is.null(consensus)) {
	consensus <- tryCatch(consensus(x, p = .5),
						  error = function(e) unroot(midpoint(superTree(x))))
  }
  if (inherits(consensus, "multiPhylo")) consensus <- consensus[[1]]

  sort_tips <- function(x) {
	x <- reorder(x)
	nTip <- as.integer(length(x$tip.label))
	e2 <- x$edge[, 2]
	x$tip.label <- x$tip.label[e2[e2 <= nTip]]
	x$edge[e2 <= nTip, 2] <- as.integer(1L:nTip)
	x
  }

  type <- match.arg(type, c("phylogram", "cladogram"))
  direction <- match.arg(direction, c("rightwards", "leftwards",  "upwards",
	"downwards"))
  horizontal <- direction %in% c("rightwards", "leftwards")

  consensus <- reorder(consensus)
  nTip <- as.integer(length(consensus$tip.label))
  consensus <- sort_tips(consensus)
  consensus <- reorder(consensus, "postorder")

  maxBT <- max(getAges(x))
  if (scaleX) maxBT <- 1.0
  label <- rev(pretty(c(maxBT, 0)))
  maxBT <- max(label)
  xy <- plotPhyloCoor(consensus, direction = direction, ...)
  yy <- xy[, 2]

  tl <- which.max(nchar(consensus$tip.label))
  if(horizontal) pin1 <- par("pin")[1]
  else pin1 <- par("pin")[2]
  sw <- strwidth(consensus$tip.label[tl], "inch", cex = cex) / pin1  * 1.1

  if (direction == "rightwards") {
	plot.default(0, type = "n", xlim = c(0, 1.0 + sw), ylim = c(0, nTip + 1),
				 xlab = "", ylab = "", axes = FALSE, ...)
	if (scale.bar) axis(side = 1, at = seq(0, 1.0, length.out = length(label)),
						labels = label)
  }
  if (direction == "leftwards") {
	plot.default(0, type = "n", xlim = c(0 - sw, 1.0), ylim = c(0, nTip + 1),
				 xlab = "", ylab = "", axes = FALSE, ...)
	if (scale.bar) axis(side = 1, at = seq(0, 1.0, length.out = length(label)),
						labels = rev(label))
  }
  if (direction == "downwards") {
	plot.default(0, type = "n", xlim = c(0, nTip + 1), ylim = c(0 - sw, 1.0),
				 xlab = "", ylab = "", axes = FALSE, ...)
	if (scale.bar) axis(side = 2, at = seq(0, 1.0, length.out = length(label)),
						labels = rev(label))
  }
  if (direction == "upwards") {
	plot.default(0, type = "n", xlim = c(0, nTip + 1), ylim = c(0, 1.0 + sw),
				 xlab = "", ylab = "", axes = FALSE, ...)
	if (scale.bar) axis(side = 2, at = seq(0, 1.0, length.out = length(label)),
						labels = label)
  }
  tip_labels <- consensus$tip.label
  if (is.expression(consensus$tip.label))
	underscore <- TRUE
  if (!underscore)
	tip_labels <- gsub("_", " ", tip_labels)

  add_tiplabels(xy, tip_labels, direction, adj = adj, font = font, srt = srt,
	cex = cex, col = tip.color, label.offset = label.offset)

  col <- rep(col, length.out = length(x))
  tiporder <-  1:nTip
  names(tiporder) <- consensus$tip.label

  if (jitter$amount > 0) {
	if (jitter$random) jit <- runif(length(x), -jitter$amount, jitter$amount)
	else jit <- seq(-jitter$amount, jitter$amount, length = length(x))
  }

  for (treeindex in seq_along(x)) {
	tmp <- reorder(x[[treeindex]])

	tmp <- sort_tips(tmp)
	#    if(!compressed) tiporder <- match(tmp$tip.label, consensus$tip.label)
	xy <- plotPhyloCoor(tmp, tip.height = tiporder, direction = direction, ...)
	xx <- xy[, 1]
	yy <- xy[, 2]

	if (horizontal) {
	  if (scaleX) xx <- xx / max(xx)
	  else xx <- xx / maxBT
	  if (direction == "rightwards") xx <- xx + (1.0 - max(xx))
	  if (jitter$amount > 0) yy <- yy + jit[treeindex]
	}
	else {
	  if (scaleX) yy <- yy / max(yy)
	  else yy <- yy / maxBT
	  if (direction == "upwards") yy <- yy + (1.0 - max(yy))
	  if (jitter$amount > 0) xx <- xx + jit[treeindex]
	}
	e1 <- tmp$edge[, 1]
	if (type == "cladogram") cladogram.plot(tmp$edge, xx, yy, edge.color =
	adjustcolor(tmp$edge.color, alpha.f = alpha), edge.width = width,
		edge.lty = lty)
	if (type == "phylogram") {
	  Ntip <- min(e1) - 1L
	  Nnode <- tmp$Nnode
	  phylogram.plot(tmp$edge, Ntip, Nnode, xx, yy, horizontal, edge.color =
		adjustcolor(tmp$edge.color, alpha.f = alpha), edge.width = width,
		edge.lty = lty)
	}
  }
}






