library(ggplot2)
library(reshape)
library(cowplot)
library(rafalib)
source("mixology/plotfuncs_mixo.R")

### 075 vs 025
load(file = "mixology/rdata/075vs025.RData")
   keep <- finalres$baseMean > 1
   png("mixology/outfigure/mixology_hist.png", width = 5, height = 5, units = "in", res = 300)
   hist(finalres$lfcs$true, xlab = "True LFC", ylab = "", main = "")
   dev.off()
   png("mixology/outfigure/mixology_scatter.png", width = 5, height = 5, units = "in", res = 300)
   plot(x = finalres$lfcs$true, y = finalres$lfcs$apeglm, xlab = "True LFC", ylab = "Estimated LFC", main = "", col=rgb(0,0,0,.3), pch = ".")
   abline(h = c(0.75, 1), col = "red")
   abline(v = c(0.75, 1), col = "red")
   dev.off()
   g1 <- plotOverLFC(finalres, keeprule = keep, error=FALSE) + 
			theme(axis.text = element_text(size=14),
				axis.title = element_text(size=22, face="bold"),
				legend.position="none", 
				axis.line = element_line(colour = "black"))
   g2 <- plotOverEST(finalres, keeprule = keep, error=FALSE) + 
			theme(axis.text = element_text(size=14),
				axis.title = element_text(size=22, face="bold"),
				legend.position="none", 
				axis.line = element_line(colour = "black"))
   rm(finalres)
   load(file = "mixology/rdata/050vs025.RData")
   g3 <- plotOverLFC(finalres, keeprule = keep, error=FALSE) + 
			theme(axis.text = element_text(size=14),
				axis.title = element_text(size=22, face="bold"),
				legend.position="none", 
				axis.line = element_line(colour = "black"))
   g4 <- plotOverEST(finalres, keeprule = keep, error=FALSE) + 
			theme(axis.text = element_text(size=14),
				axis.title = element_text(size=22, face="bold"),
				legend.position="none", 
				axis.line = element_line(colour = "black"))
   png("mixology/outfigure/mixology.png", width = 10, height = 10, units = "in", res = 300)
   plot_grid(g1, g2, g3, g4, labels="auto", hjust = 0, vjust = 1, scale = c(1,1,1,1))
   dev.off()