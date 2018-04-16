library(rafalib)
library(reshape2)
library(ggplot2)
library(dplyr)
library(devtools)
library(cowplot)

source("common/plotfuncs.R")

load("highrep/rdata/analysis_yeast_3v3.RData")

res3 <- res
mle3 <- fullmle
mu3 <- full.dres$baseMean
keep3 <- mu3 > 1
mle3.keep <- as.vector(mle3[keep3])
mu3.keep <- as.vector(mu3[keep3])

# number of iterations
iters <- length(res3)

rm(res, fullmle, full.dres)

load("highrep/rdata/analysis_yeast_5v5.RData")

res5 <- res
mle5 <- fullmle
mu5 <- full.dres$baseMean
keep5 <- mu5 > 1
mle5.keep <- as.vector(mle5[keep5])
mu5.keep <- as.vector(mu5[keep5])

rm(res, fullmle, full.dres)

## Needs to study Mike's iCOBRA code ##
## load_all('iCOBRA')

# Plots 
g1 <- plotOverLFC(res3, keeprule = keep3, error=FALSE) + 
			theme(axis.text = element_text(size=14),
				axis.title = element_text(size=22, face="bold"),
				legend.position="none", 
				axis.line = element_line(colour = "black"))

g2 <- catPlot(res3, keeprule = keep3, error=FALSE) + 
		  scale_y_continuous(breaks=1:50/50) + 
			 theme(axis.text = element_text(size=14),
				axis.title = element_text(size=22, face="bold"),
				legend.position="none", 
				axis.line = element_line(colour = "black"))
g3 <- plotOverLFC(res5, keeprule = keep5, error=FALSE) + 
			theme(axis.text = element_text(size=14),
				axis.title = element_text(size=22, face="bold"),
				legend.position="none", 
				axis.line = element_line(colour = "black"))

g4 <- catPlot(res5, keeprule = keep5, error=FALSE) + 
		  scale_y_continuous(breaks=1:50/50) + 
			 theme(axis.text = element_text(size=14),
				axis.title = element_text(size=22, face="bold"),
				legend.position="none", 
				axis.line = element_line(colour = "black"))	
png("highrep/outfigure/yeast.png", 
    height = 10, width = 10, units = 'in', res = 300)
    plot_grid(g1, g2, g3, g4, labels="auto", hjust = 0, vjust = 1, scale = c(1,1,1,1))
dev.off()

png("highrep/outfigure/sup_est_tru_1iter_3.png", 
     height = 5, width = 5, units = 'in', res = 300)
	plotLFCs(lfcs = res3[[1]]$lfcs, lim = 1.5)
dev.off()

lfcs <- lapply(res3, function(x) x$lfcs)

png("highrep/outfigure/yeast_arrow_small_3.png", 
     height = 5, width = 6, units = 'in', res = 300)
	ArrowPlot(lfcs = lfcs, mu.obs = as.vector(mu3))
dev.off()

png("highrep/outfigure/yeast_overmean_3.png", 
    height = 5, width = 6, units = 'in', res = 300)
	plotOverMean(res3, keeprule = keep3, mu.obs = mu3.keep, error=FALSE)
dev.off()


png("highrep/outfigure/sup_est_tru_1iter_5.png", 
     height = 5, width = 5, units = 'in', res = 300)
	plotLFCs(lfcs = res5[[1]]$lfcs, lim = 1.5)
dev.off()

1500

rm(lfcs)
lfcs <- lapply(res5, function(x) x$lfcs)

png("highrep/outfigure/yeast_arrow_small_5.png", 
     height = 5, width = 6, units = 'in', res = 300)
	ArrowPlot(lfcs = lfcs, mu.obs = as.vector(mu5))
dev.off()

png("highrep/outfigure/yeast_overmean_5.png", 
    height = 5, width = 6, units = 'in', res = 300)
	plotOverMean(res5, keeprule = keep5, mu.obs = mu5.keep, error=FALSE)
dev.off()