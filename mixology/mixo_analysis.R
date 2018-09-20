library(apeglm)
library(DESeq2)
library(edgeR)
library(limma)
library(ashr)
library(dplyr)
library(devtools)
library(REBayes)

source("common/funcs.R")

load("mixology/RSCE.mRNA.unstranded.rda")

RSCE.mRNA.unstranded$genes <- RSCE.mRNA.unstranded$genes[,c("GeneID", "Symbols")]

data <- RSCE.mRNA.unstranded
data <- calcNormFactors(data)
batch <- strsplit2(colnames(data), split="\\.")[,1]
mixes <- as.numeric(strsplit2(strsplit2(colnames(data), split="\\.")[,2], split="_")[,1])/100
design <- model.matrix(~data$samples$group-1)
v <- voom(data, design, plot=FALSE)
remove <- batch=="R2D"
fitmixgood <- fitmixture(v$E[,!remove], mixes[!remove])
c12 <- 2^(fitmixgood$M/2+fitmixgood$A)
c22 <- 2^(fitmixgood$A-fitmixgood$M/2)
rm(data, batch, mixes, design, v, remove, fitmixgood)

comparisons <- c("075vs025", "050vs025")
sample <- strsplit2(RSCE.mRNA.unstranded$samples$Sample, split="\\.")
colnames(sample) <- c("Batch", "Mixture")
sample <- as.data.frame(sample)

batch <- c("R1", "R2", "R3")

for (c in comparisons){
   cat("== Running comparison", c, "== \n")
   select <- c(strsplit2(c, split="vs"))
   sample.select <- which(sample$Mixture %in% select & sample$Batch %in% batch)
   data <- RSCE.mRNA.unstranded[, sample.select]   
   
   levels(data$samples$group) <- c("low", "high")
   coldata <- data$samples
   coldata$condition <- coldata$group

   conc <- as.numeric(select)
   p <- conc[1]/100
   q <- conc[2]/100
   #   Adjust predicted logFC based on sample concentrations
   logfc.adjusted <- log2(p*c12+(1-p)*c22)-log2(q*c12+(1-q)*c22)
   
   dds <- DESeqDataSetFromMatrix(data$counts, colData = coldata, design = ~condition)
   mcols(dds) <- DataFrame(mcols(dds), data$genes)

   res <- getLFCsSvals(test = dds, beta = logfc.adjusted)
  
   mle.dds <- DESeq(dds)
   mle.res <- results(mle.dds)
   mu <- mle.res$baseMean

   finalres <- list(lfcs = res$lfcs, sval = res$sval, baseMean = mu)
   save(finalres, file = paste0("mixology/rdata/", c, ".RData"))
 }
