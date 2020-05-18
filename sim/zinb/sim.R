knitr::opts_chunk$set(cache=TRUE)
deprob <- 0.01/3
groupprob <- 0.2
name <- paste0("deprob_", deprob*3, "_group_", groupprob*100)
suppressPackageStartupMessages(library(splatter))

cddir <- "/proj/milovelab/zhu/projects/ape/lfcprior2/eval/code_release/"
source(paste0(cddir, "common/plotfuncs.R"))
source(paste0(cddir, "sim/loglikZINB.R"))
source(paste0(cddir, "sim/generateZINB.R"))

library(zinbwave)
library(apeglm)
library(ZIM)

sim <- generateZINB(deprob = deprob, groupprob = groupprob)

keep <- rowSums(counts(sim) >= 10) >= 5
table(keep)
zinb <- sim[keep,]
zinb$condition <- factor(zinb$Group)
nms <- c("counts", setdiff(assayNames(zinb), "counts"))
assays(zinb) <- assays(zinb)[nms]
# epsilon as recommended from Van den Berge and Perraudeau
zinb <- zinbwave(zinb, K=0, BPPARAM=SerialParam(), epsilon=1e12)

suppressPackageStartupMessages(library(DESeq2))
dds <- DESeqDataSet(zinb, design=~condition)
# arguments as recommended from Van den Berge and Perraudeau
dds <- DESeq(dds, test="LRT", reduced=~1,
             sfType="poscounts", minmu=1e-6, minRep=Inf)

## ----plotdisp------------------------------------------------------------
plotDispEsts(dds)

## ------------------------------------------------------------------------
wts <- assays(dds)[["weights"]]
ncts <- counts(dds, normalized=TRUE)
idx1 <- dds$condition == "Group1"
idx2 <- dds$condition == "Group2"
pc <- .1
simple.lfc <- log2(rowMeans(ncts[,idx2])+pc) -
  log2(rowMeans(ncts[,idx1])+pc)
wtd.lfc <- log2(rowSums(ncts[,idx2]*wts[,idx2])/rowSums(wts[,idx2])+pc) -
  log2(rowSums(ncts[,idx1]*wts[,idx1])/rowSums(wts[,idx1])+pc)

## ------------------------------------------------------------------------
res <- results(dds, name="condition_Group2_vs_Group1",　
               independentFiltering=FALSE)

## ------------------------------------------------------------------------
res.norm <- lfcShrink(dds, coef=2, type="normal")

## ------------------------------------------------------------------------
ape.nb <- lfcShrink(dds, coef=2, type="apeglm")

## ------------------------------------------------------------------------
Y <- counts(dds)
design <- model.matrix(design(dds), data=colData(dds))
disps <- dispersions(dds)
# combine dispersion and wts into a parameter matrix,
# which will be passed to apeglm
param <- cbind(disps, 1 - wts)
offset <- matrix(log(sizeFactors(dds)), nrow=nrow(dds),
                 ncol=ncol(dds), byrow=TRUE)
# need to put to natural log scale for apeglm
mle <- log(2) * cbind(res$log2FoldChange, res$lfcSE)
# run apeglm with a ZINB likelihood and zinbwave weights
# used to define the probability of an excess zero
fit <- apeglm(Y=Y, x=design, log.lik=logLikZINB, param=param,
              coef=2, mle=mle, offset=offset)
# need to put back to log2 scale
ape.zinb <- log2(exp(1)) * fit$map[,2]

## ------------------------------------------------------------------------

## ----plotlfc, fig.width=9, fig.height=7----------------------------------

png(paste0("/proj/milovelab/zhu/projects/ape/lfcprior2/eval/sim_ZINB/ZINB_MAE_", name, ".png"), height = 2000, width = 3000, res = 300)
par(mfrow=c(2,3), mar=c(2,3,2,1))
myplot(mcols(dds)$log2FC, simple.lfc, main="pseudocount")
myplot(mcols(dds)$log2FC, wtd.lfc, main="wtd pseudocount")
myplot(mcols(dds)$log2FC, res$log2FoldChange, main="MLE")
myplot(mcols(dds)$log2FC, res.norm$log2FoldChange, main="normal + wtd NB")
myplot(mcols(dds)$log2FC, ape.nb$log2FoldChange, main="apeglm + wtd NB")
myplot(mcols(dds)$log2FC, ape.zinb, main="apeglm + ZINB lik")
dev.off()
save.image(file = paste0("/proj/milovelab/zhu/projects/ape/lfcprior2/eval/sim_ZINB/", name, ".RData"))

## ------------------------------------------------------------------------
sessionInfo()

