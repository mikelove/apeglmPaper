eval.sim <- function(Y, condition, beta){
  dds <- DESeqDataSetFromMatrix(Y, DataFrame(condition), ~condition)
  res <- getLFCsSvals(test = dds, beta = beta)
  finalres <- list(lfcs = res$lfcs, sval = res$sval, Y = Y)
  finalres
}
