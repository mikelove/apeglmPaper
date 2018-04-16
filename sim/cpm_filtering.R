suppressPackageStartupMessages(library(DESeq2))
library(edgeR)

load("bottomly_sumexp.RData")
bottomly.se <- updateObject(bottomly)
levels(bottomly.se$strain) <- c("C","D")
bottomly.se$batch <- factor(bottomly.se$experiment.number)
bottomly.se <- DESeqDataSet(bottomly.se, ~batch + strain)
system.time( bottomly.se <- DESeq(bottomly.se) )
res <- results(bottomly.se, alpha=.05)

cts <- counts(bottomly.se)
ncts <- counts(bottomly.se, normalized=TRUE)

keep <- matrix(NA, nrow=nrow(cts), ncol=100)
lfcs <- matrix(NA, nrow=nrow(cts), ncol=100)
for (i in 1:100) {
  set.seed(i)
  idx <- with(colData(bottomly.se), 
              c(sample(which(strain == "C" & batch == "4"), 1),
                sample(which(strain == "C" & batch == "6"), 1),
                sample(which(strain == "C" & batch == "7"), 1),
                sample(which(strain == "D" & batch == "4"), 1),
                sample(which(strain == "D" & batch == "6"), 1),
                sample(which(strain == "D" & batch == "7"), 1)))
  cts.sub <- cts[,idx]
  ncts.sub <- ncts[,idx]
  y <- DGEList(counts=cts.sub)
  y <- calcNormFactors(y)
  cpm <- cpm(y)
  L <- min(colSums(cts)/1e6)
  keep[,i] <- rowSums(cpm > 10/L) >= 3
  # the mean of per batch LFCs
  lfcs[,i] <- rowMeans(cbind(
    log2(ncts.sub[,4] + .125) - log2(ncts.sub[,1] + .125),
    log2(ncts.sub[,5] + .125) - log2(ncts.sub[,2] + .125),
    log2(ncts.sub[,6] + .125) - log2(ncts.sub[,3] + .125)))    
}

padj <- res$padj
padj[is.na(padj)] <- 1

rm.ncts <- rowMeans(ncts+.125)
filtered.sig <- rowMeans(keep) < .5 & padj < .05 & rm.ncts > 10
sum(filtered.sig)

wfs <- which(filtered.sig)
sign.correct <- rowMeans(sign(lfcs[wfs,]) == sign(res$log2FoldChange[wfs]))
mean(sign.correct)

library(rafalib)

bigpar()
plot(rm.ncts, rowMeans(keep), log="x", xlim=c(5,50),
     xlab="mean of normalized counts", ylab="how often passing CPM filter",
     col=ifelse(filtered.sig, rgb(1,0,0,.5), rgb(0,0,0,.5)),
     pch=20, cex=1.5)

##

bigpar()
hist(sign.correct, xlab="how often sign of LFC is correct",
     ylab="number of genes", main="")

##

bigpar()
plot(rowMeans(keep), -log10(res$pvalue),
     col=ifelse(filtered.sig, rgb(1,0,0,.5), rgb(0,0,0,.5)),
     ylim=c(0,50), xlim=c(0,.55))
pts.idx <- identify(rowMeans(keep), -log10(res$pvalue))
pts.idx <- c(11838, 16061, 24854, 36672)
which(colSums(keep[pts.idx,]) == 4)[1]

set.seed(56)
idx <- with(colData(bottomly.se), 
            c(sample(which(strain == "C" & batch == "4"), 1),
              sample(which(strain == "C" & batch == "6"), 1),
              sample(which(strain == "C" & batch == "7"), 1),
              sample(which(strain == "D" & batch == "4"), 1),
              sample(which(strain == "D" & batch == "6"), 1),
              sample(which(strain == "D" & batch == "7"), 1)))
ncts.sub <- ncts[,idx]

bigpar(2,2)
for (i in pts.idx) {
  plotCounts(bottomly.se, i, "strain", transform=FALSE, cex=2)
  points(rep(1:2,each=3), ncts.sub[i,], col="red", pch=4, cex=2, lwd=2)
}
