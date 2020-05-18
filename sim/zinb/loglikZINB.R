logLikZINB <- function (y, x, beta, param, offset) {
    xbeta <- x %*% beta + offset
    mean.hat <- exp(xbeta)
    k <- 1/param[1]
    omega <- param[-1]
    out <- ZIM::dzinb(y, k=k, lambda=mean.hat, omega=omega, log=TRUE)
    out[!is.finite(out)] <- -1000
    out
}
