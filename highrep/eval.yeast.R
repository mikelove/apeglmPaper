eval.yeast <- function(i, fulldata, ind.one, ind.two, sampsize, MLE){
  tid <- c(sample(ind.one, sampsize, replace=F),sample(ind.two, sampsize, replace=F))
  test <- fulldata[,tid]
  res <- getLFCsSvals(test = test, beta = MLE)
  res
}
