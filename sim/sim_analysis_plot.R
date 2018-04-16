library(apeglm)
library(DESeq2)
library(edgeR)
library(limma)
library(ashr)
library(dplyr)
library(devtools)
library(REBayes)


### Analysis

load("eval/sim/bottomly/rdata/bottomly_sumexp.RData")
bottomly.se <- updateObject(bottomly)
levels(bottomly.se$strain) <- c("C","D")
bottomly.se$batch <- factor(bottomly.se$experiment.number)
bottomly.se <- DESeqDataSet(bottomly.se, ~batch + strain)
bottomly.se <- estimateSizeFactors(bottomly.se)
bottomly.se <- estimateDispersions(bottomly.se)

source("common/funcs.R")
source("common/plotfuncs.R")
source("common/eval.sim.R")
source("common/generateNB.R")

# from Pickrell
library(BiocParallel)
BPPARAM <- MulticoreParam(4)

load("eval/pickrell/rdata/meanDispPairs.RData")
is.pickrell <- TRUE

niter <- 10

pickrell5 <- bplapply(1:niter, function(i) {
  set.seed(i)
  simOut <- generateNB(meanDispPairs, is.pickrell, nper = 5)
  Y <- simOut$Y
  condition <- simOut$condition
  m <- nrow(Y); n <- ncol(Y)
  beta <- simOut$beta
  lfcs <- eval.sim(Y, condition, beta)
}, BPPARAM=BPPARAM)

pickrell10 <- bplapply(1:niter, function(i) {
  set.seed(i)
  simOut <- generateNB(meanDispPairs, is.pickrell, nper = 10)
  Y <- simOut$Y
  condition <- simOut$condition
  m <- nrow(Y); n <- ncol(Y)
  beta <- simOut$beta
  lfcs <- eval.sim(Y, condition, beta)
}, BPPARAM=BPPARAM)

# from Bottomly
meanDispPairs <- with(mcols(bottomly.se),
                      data.frame(x=baseMean, y=dispGeneEst))
meanDispPairs <- meanDispPairs[meanDispPairs$y > 1e-6,]
meanDispPairs <- meanDispPairs[!is.na(meanDispPairs$x),]
is.pickrell <- FALSE

bottomly5 <- bplapply(1:niter, function(i) {
  set.seed(i)
  simOut <- generateNB(meanDispPairs, is.pickrell, nper = 5)
  Y <- simOut$Y
  condition <- simOut$condition
  m <- nrow(Y); n <- ncol(Y)
  beta <- simOut$beta
  lfcs <- eval.sim(Y, condition, beta)
}, BPPARAM=BPPARAM)

bottomly10 <- bplapply(1:niter, function(i) {
  set.seed(i)
  simOut <- generateNB(meanDispPairs, is.pickrell, nper = 10)
  Y <- simOut$Y
  condition <- simOut$condition
  m <- nrow(Y); n <- ncol(Y)
  beta <- simOut$beta
  lfcs <- eval.sim(Y, condition, beta)
}, BPPARAM=BPPARAM)

save(pickrell5, bottomly5, pickrell10, bottomly10, file="eval/sim/rdata/res.rda")

load("eval/sim/rdata/res.rda")

### Plot 

library(ggplot2)
library(reshape)
library(cowplot)
library(rafalib)

png("eval/sim/outfigure/sim_p.png", width = 10, height = 10, units = "in", res = 300)
g1 <- plotOverLFC(pickrell5, error=FALSE) + 
			theme(axis.text = element_text(size=14),
				axis.title = element_text(size=22, face="bold"),
				legend.position="none", 
				axis.line = element_line(colour = "black"))
g2 <- catPlot(pickrell5, error=FALSE) + scale_y_continuous(breaks=1:50/50) + 
theme(axis.text = element_text(size=14),
				axis.title = element_text(size=22, face="bold"),
				legend.position="none",
				axis.line = element_line(colour = "black"))
g3 <- plotOverLFC(pickrell10, error=FALSE) + 
			theme(axis.text = element_text(size=14),
				axis.title = element_text(size=22, face="bold"),
				legend.position="none", 
				axis.line = element_line(colour = "black"))
g4 <- catPlot(pickrell10, error=FALSE) + scale_y_continuous(breaks=1:50/50) + 
theme(axis.text = element_text(size=14),
				axis.title = element_text(size=22, face="bold"),
				legend.position="none",
				axis.line = element_line(colour = "black"))
plot_grid(g1, g2, g3, g4, labels="auto", hjust = 0, vjust = 1, scale = c(1,1,1,1))
dev.off()


png("eval/sim/outfigure/sim_error_mean.png", width = 10, height = 10, units = "in", res = 300)
# Pickrell

g1 <- plotOverMean(pickrell5, error=FALSE) +  
			theme(axis.text = element_text(size=14),
				axis.title = element_text(size=22, face="bold"),
				legend.position="none",
				axis.line = element_line(colour = "black"))

g2 <- plotOverMean(pickrell10, error=FALSE) +  
			theme(axis.text = element_text(size=14),
				axis.title = element_text(size=22, face="bold"),
				legend.position="none",
				axis.line = element_line(colour = "black"))
# Bottomly
g3 <- plotOverMean(bottomly5, error=FALSE) + 
			theme(axis.text = element_text(size=14),
				axis.title = element_text(size=22, face="bold"),
				legend.position="none", 
				axis.line = element_line(colour = "black"))

g4 <- plotOverMean(bottomly10, error=FALSE) + 
			theme(axis.text = element_text(size=14),
				axis.title = element_text(size=22, face="bold"),
				legend.position="none", 
				axis.line = element_line(colour = "black"))				
plot_grid(g1, g2, g3, g4, labels="auto", hjust = 0, vjust = 1, scale = c(1,1,1,1))
dev.off()

png("eval/sim/outfigure/sim_b.png", width = 10, height = 10, units = "in", res = 300)
g1 <- plotOverLFC(bottomly5, error=FALSE) + 
			theme(axis.text = element_text(size=14),
				axis.title = element_text(size=22, face="bold"),
				legend.position="none", 
				axis.line = element_line(colour = "black"))
g2 <- catPlot(bottomly5, error=FALSE) + scale_y_continuous(breaks=1:50/50) + 
theme(axis.text = element_text(size=14),
				axis.title = element_text(size=22, face="bold"),
				legend.position="none",
				axis.line = element_line(colour = "black"))
g3 <- plotOverLFC(bottomly10, error=FALSE) + 
			theme(axis.text = element_text(size=14),
				axis.title = element_text(size=22, face="bold"),
				legend.position="none", 
				axis.line = element_line(colour = "black"))
g4 <- catPlot(bottomly10, error=FALSE) + scale_y_continuous(breaks=1:50/50) + 
theme(axis.text = element_text(size=14),
				axis.title = element_text(size=22, face="bold"),
				legend.position="none",
				axis.line = element_line(colour = "black"))
plot_grid(g1, g2, g3, g4, labels="auto", hjust = 0, vjust = 1, scale = c(1,1,1,1))
dev.off()


png("eval/sim/outfigure/sim_scatter_p_5.png", width = 5, height = 5, units = "in", res = 300)
plotLFCs(pickrell5[[1]]$lfcs, 6)
dev.off()

png("eval/sim/outfigure/sim_scatter_b_5.png", width = 5, height = 5, units = "in", res = 300)
plotLFCs(bottomly5[[1]]$lfcs, 3)
dev.off()


png("eval/sim/outfigure/sim_scatter_p_10.png", width = 5, height = 5, units = "in", res = 300)
plotLFCs(pickrell10[[1]]$lfcs, 6)
dev.off()

png("eval/sim/outfigure/sim_scatter_b_10.png", width = 5, height = 5, units = "in", res = 300)
plotLFCs(bottomly10[[1]]$lfcs, 3)
dev.off()
