generateNB <- function(meanDispPairs, is.pickrell, nper) {
  meanDispPairs <- meanDispPairs[meanDispPairs[,1] > 10,]
  n.per.group <- nper # number samples per group
  n <- n.per.group * 2
  m <- 10000 # number of rows
  condition <- factor(rep(letters[1:2],each=n.per.group))
  design <- model.matrix(~ condition)
  # only positive LFC
  if (is.pickrell) {
    beta <- abs(c(rnorm(.05*m,0,3),rnorm(.05*m,0,2),rnorm(.9*m,0,1)))
  } else {
    beta <- abs(c(rnorm(.05*m,0,1),rnorm(.05*m,0,.5),rnorm(.9*m,0,.25)))
  }
  idx <- sample(nrow(meanDispPairs), m, replace=TRUE)
  mu0 <- meanDispPairs[idx, 1]
  disp <- meanDispPairs[idx, 2]
  #plot(mu0, disp, log="xy")
  sf <- rep(1,n)
  betafull <- cbind(log2(mu0), beta)
  mu <- 2^(betafull %*% t(design))
  # now swap half of the betas to be negative
  swap <- sort(sample(m, m/2, replace=FALSE))
  mu[swap,] <- mu[swap,c((n.per.group+1):n,1:n.per.group)]
  beta[swap] <- -1 * beta[swap]
  muMat <- matrix(rep(mu, times=n) * rep(sf, each=m), ncol=n)
  Y <- matrix(rnbinom(m*n, mu=muMat, size=1/disp), ncol=n)
  list(Y=Y, beta=beta, condition=condition)
}
