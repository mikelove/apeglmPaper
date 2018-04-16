dir = "highrep/rdata"
setwd(dir)

library(DESeq2)
library(edgeR)
library(limma)
library(parallel)
library(apeglm)
source("common/funcs.R")
source("highrep/eval.yeast.R")
library(ashr)
library(REBayes)

# Here we have a DESeq data for all genes that has non-zero count (fulldata), 
# and the DESeq Results for the full data set (full.dres)
load("highrep/rdata/ddsdata_fulldres.RData")

# number of genes
ngene <-  nrow(fulldata)

# MLE on full data set
fullmle <- data.matrix(full.dres$lfcMLE)
rownames(fullmle) <- rownames(full.dres)

# Sampling and iterations
ind.mut <- which(tolower(ddsdata$condition)=="mut")    #44 samples
ind.wt <- which(tolower(ddsdata$condition)=="wt")      #42 samples

samplesize <- 3
iterations <- 100

set.seed(26)

options(mc.cores=12)
res <- mclapply(1:iterations, eval.yeast, fulldata = fulldata, ind.one = ind.mut, ind.two = ind.wt, sampsize = samplesize, MLE = full.dres$lfcMLE)
save.image(file = "analysis_yeast_3v3.RData")

sessionInfo()
