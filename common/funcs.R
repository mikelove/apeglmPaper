##############################################
## This is the R script of functions used for
## the evaluation of apeglm (call from DESeq2) 
## lfcShrink() with yeast data. 
## For previous version of evaluation see
## lfcprior/eval/yeast/code/kdcode/funcs.R
## to run any other method, simply load the 
## .RData and use tid
##############################################

# run edgeR function
edger <- function(dds, pc = NULL){
   y <- DGEList(counts=counts(dds), group=dds$condition)
   y <- calcNormFactors(y)
   design <- model.matrix(~dds$condition)
   y <- estimateDisp(y,design)
  if (is.null(pc)){
    fite <- glmFit(y, design)
  } else {
	fite <- glmFit(y, design, prior.count = pc)
  }
  lrt <- glmLRT(fite)
  tte <- topTags(lrt,sort="none",n=Inf)
  return(tte)
}

# run ashr function
ashr.shrinklfc <- function(beta, betase){
   nm <- ash(beta, betase, mixcompdist = "normal", method="shrink")  #
  nm
}

# limma-voom function
limmav <- function(dds){
  design <- model.matrix(~dds$condition)
  dgel <- DGEList(counts(dds))
  dgel <- calcNormFactors(dgel)
  v <- voom(dgel,design,plot=FALSE)
  fitlma <- lmFit(v,design)
  fitebs <- eBayes(fitlma)
  fitebs 
}

# DESeq2 function
deseq2 <- function(dds) {
  dds <- nbinomWaldTest(dds, betaPrior=TRUE)
  results(dds)
}


getLFCsSvals <- function(test, beta){
  ## test is a dataset of DESeq format
  
  ## MLE estimates
  mle.dds <- DESeq(test)
  mle.res <- results(mle.dds)
  mle.res$log2FoldChange[mle.res$baseMean == 0] <- 0

  # apeglm with Adaptive Shrinkage
  aperes <- lfcShrink(dds = mle.dds, coef = 2, type = "apeglm", svalue = TRUE, returnList = FALSE, apeAdapt = TRUE)
  aperes$log2FoldChange[mle.res$baseMean == 0] <- 0

  # DESeq2 
  deseq2res <- lfcShrink(dds = mle.dds, coef=2)	
  deseq2res$log2FoldChange[mle.res$baseMean == 0] <- 0

  # edgeR 
  test.edgerres <- edger(test)
  
  # edgeR with prior count of 5
  test.edgerpcres <- edger(test, pc = 5)

  # ashr DESeq2 input
  ashrdeseq2.res <- lfcShrink(dds = mle.dds, coef = 2, type = "ashr", svalue = TRUE, returnList = FALSE)
							 
  # limma-voom
  
  test.limres <- limmav(test)
  limmares <- topTable(test.limres,sort="none",n=Inf)

  # ashr based on limma results
  
  ashr.res2 <- ashr.shrinklfc(topTable(test.limres,sort="none",n=Inf)$logFC,
                              sqrt(test.limres$s2.post)*test.limres$stdev.unscaled[,2])

  lfcs <- list()
  sval <- list()

  lfcs[["true"]] <- beta
  lfcs[["apeglm"]] <- aperes$log2FoldChange
  sval[["apeglm"]] <- aperes$svalue

  lfcs[["ashr.d"]] <- ashrdeseq2.res$log2FoldChange
  sval[["ashr.d"]] <- ashrdeseq2.res$svalue

  lfcs[["ashr.l"]] <- ashr.res2$result$PosteriorMean
  sval[["ashr.l"]] <- ashr::get_svalue(ashr.res2)
  
  lfcs[["DESeq2"]] <- deseq2res$log2FoldChange
  sval[["DESeq2"]] <- apeglm::svalue(0.5*deseq2res$pvalue)
  
  lfcs[["edgeRPC5"]] <- test.edgerpcres$table$logFC
  lfcs[["edgeR"]] <- test.edgerres$table$logFC
  lfcs[["limma"]] <- limmares$logFC
  
  return(list(lfcs=lfcs, sval=sval))

}
