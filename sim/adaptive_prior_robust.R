set.seed(1001)
g <- gl(2, 5)
collected <- numeric(10000)
for (it in seq_along(collected)) {
    y <- rpois(10, lambda=1)
    fit <- glm(y~g, family=poisson)
    collected[it] <- summary(fit)$coefficients[,2][2]
    stopifnot(collected[it] < 1000)
}
hist(collected) # a number of strong outliers
hist(collected[collected < 100]) # still very dispersed

###

library(apeglm)
n <- 1000 # number of genes
nsamp <- 10
grid <- expand.grid(A=c(.25,.5,.75,1),lambda=c(1,5,10))
out <- lapply(seq_len(nrow(grid)), function(t) {
  A <- grid$A[t]
  res <- sapply(1:10, function(i) {
    set.seed(i)
    theta <- rnorm(n, 0, sqrt(A))
    mle <- matrix(ncol=2, nrow=n)
    for (j in 1:n) {
      lambdas <- rep(grid$lambda[t],nsamp) *
        exp(rep(c(0,theta[j]),each=nsamp/2))
      x <- factor(rep(1:2,each=nsamp/2))
      y <- rpois(nsamp, lambda=lambdas)
      mle[j,] <- summary(glm(y ~ x, poisson))$coefficients[2,1:2]
    }
    apeglm:::priorVar(mle)
  })
  data.frame(A=rep(grid$A[t],length(res)),
             lambda=rep(grid$lambda[t],length(res)),
             A.hat=res)
})
dat <- do.call(rbind, out)
dat$A <- factor(dat$A)
true.pts <- grid
true.pts$A.hat <- true.pts$A
true.pts$A <- factor(true.pts$A)
library(ggplot2)
library(cowplot)

pdf(file="unstable_se.pdf", width=10, height=4)
g <- ggplot(dat, aes(A, A.hat)) +
  geom_boxplot() + facet_wrap(~lambda, labeller=label_both) +
  geom_point(data=true.pts, color="blue", size=3) +
  theme(axis.text = element_text(size=14),
        axis.title = element_text(size=22, face="bold"))
plot_grid(g)
dev.off()

###

library(apeglm)
n <- 5000 # number of genes
nsamp <- 10
grid <- expand.grid(A=c(.25,.5,.75,1),zero.center=c(.8,.9,.95))
out <- lapply(seq_len(nrow(grid)), function(t) {
  A <- grid$A[t]
  res <- t(sapply(1:10, function(i) {
    set.seed(i)
    theta <- c(rnorm(n*grid$zero.center[t], 0, sqrt(A)),
               rnorm(n*(1-grid$zero.center[t]), 3*sqrt(A), sqrt(A)/2))
    #hist(theta, col="grey",breaks=30)
    D <- rexp(n, rate=.1)
    X <- rnorm(n, theta, sqrt(D))
    c(apeglm:::priorVar(cbind(X,sqrt(D))), var(theta))
  }))
  data.frame(A=rep(grid$A[t],nrow(res)),
             zero.center=rep(grid$zero.center[t],nrow(res)),
             A.hat=res[,1],
             var=res[,2])
})
dat <- do.call(rbind, out)
dat$A <- factor(dat$A)
true.pts <- grid
true.pts$A.hat <- true.pts$A
true.pts$A <- factor(true.pts$A)
library(dplyr)
var <- dat %>% group_by(A,zero.center) %>% summarize(A.hat=mean(var))
library(ggplot2)
cols <- c("total"="red","zero.center"="blue")

pdf(file="adaptive_bimodal.pdf", width=10, height=4)
ggplot(dat, aes(A, A.hat)) +
  geom_boxplot() + facet_wrap(~zero.center, labeller=label_both) +
  geom_point(data=true.pts, aes(col="zero.center"), size=3) + 
  geom_point(data=var, aes(col="total"), size=3) +
  scale_color_manual("Var",values=cols) +
  theme(axis.text = element_text(size=14),
        axis.title = element_text(size=22, face="bold"))
dev.off()
