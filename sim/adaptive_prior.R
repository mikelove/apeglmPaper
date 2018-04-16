sim <- function(s, m, n, dispersion) {
  n.per.group <- n/2
  condition <- factor(rep(letters[1:2],each=n.per.group))
  design <- model.matrix(~ condition)  
  beta.sd <- s
  beta.b <- rnorm(m,0,beta.sd)
  intercept <- runif(m,1,10)
  beta.mat <- cbind(intercept, beta.b)
  mu <- exp(t(design %*% t(beta.mat)))
  Y <- matrix(rnbinom(m*n, size=1/dispersion, mu=mu),ncol=n)
  list(Y=Y, beta.mat=beta.mat)
}

# prior.control for non-adaptive shrinkage
pc <- list(
  no.shrink = 1,
  prior.mean = 0,
  prior.scale = 1,
  prior.df = 1,
  prior.no.shrink.mean = 0,
  prior.no.shrink.scale = 15
)
pc0 <- pc

fsrCalc <- function(f, dat) {
  sig <- f$svalue < alpha
  if (sum(sig, na.rm=TRUE) == 0) {
    0
  } else {
    mean(sign(f$map[sig,2]) != sign(dat$beta.mat[sig,2]), na.rm=TRUE)
  }
}

alpha <- .01 # testing at FSR = 0.01
m <- 500 # number of rows
# the s-value inflation can be seen when there is only small biological variation
#dispersion <- .0001
dispersion <- .01
# a grid of standard deviations for the simulated betas: N(0,s^2)
s <- c(.01,.025,.05,.1,.25,.5,1)
# we will see how the scale of t prior links to scale of simulated betas
# `scale` = 1 is linking a t prior with scale=1 to a distribution of betas with SD=1
# `scale` > 1 is having a wider scale for the t prior
# `scale` < 1 is having a narrower scale for the t prior
scale <- c(.1, .25, .5, 1, 2, 4)
# a grid for the per-group sample size
# n=3 we notice not much error, n=100 big error
ns <- c(5,10,50)
grid <- expand.grid(ns=ns,s=s,scale=scale)
grid <- data.frame(lapply(grid, rep, each=100))

library(apeglm)

#library(pbapply)
# takes ~30 minutes on laptop with each=10
#res <- pbsapply(seq_len(nrow(grid)), function(i) {

library(BiocParallel)
BPPARAM <- MulticoreParam(8)
res <- bplapply(seq_len(nrow(grid)), function(i) {

  set.seed(i)
  n.per.group <- grid$ns[i]
  n <- n.per.group * 2
  offset <- matrix(0, nrow=m, ncol=n)
  condition <- factor(rep(letters[1:2],each=n.per.group))
  design <- model.matrix(~ condition)
  dat <- sim(grid$s[i],m,n,dispersion)
  Y <- dat$Y
  param <- rep(dispersion, nrow(Y))
  # scale the prior based on oracle knowledge of beta distribution
  pc$prior.scale <- grid$scale[i] * grid$s[i]
  # fit once with adaptive prior, once with constant scale = 1
  fit.adapt <- apeglm(Y=Y,x=design,log.lik=logLikNB,
                      param=param,coef=2,prior.control=pc,offset=offset)
  fit.const <- apeglm(Y=Y,x=design,log.lik=logLikNB,
                      param=param,coef=2,prior.control=pc0,offset=offset)
  # the "error inflation factor" = how much more MAE on betas
  # when we adapt prior to the scale of the beta distribution
  mae.adapt <- median(abs(dat$beta.mat[,2] - fit.adapt$map[,2]))
  mae.const <- median(abs(dat$beta.mat[,2] - fit.const$map[,2]))
  eif <- mae.adapt/mae.const
  # FSR for adaptive prior and for constant scale
  fsrs <- lapply(list(fit.adapt, fit.const), fsrCalc, dat)
  # return these three measures
  c(eif, fsrs[[1]], fsrs[[2]])
}, BPPARAM=BPPARAM)
#})

#save(grid, res, file="res.rda")

# re-assemble results into a data.frame with grid parameters
load("res.rda")

#res <- t(res)
res <- do.call(rbind, res)

colnames(res) <- c("eif","fsr","fsr.const")
gres <- data.frame(grid, res)
gres$scale <- factor(gres$scale)

fun.y <- median
fun.ymax <- function(x) quantile(x,.75)
fun.ymin <- function(x) quantile(x,.25)

library(ggplot2)
library(cowplot)
# the problem: the FSR with fixed prior scale=1
brks <- unique(gres$s)
alpha <- 0.01
g <- ggplot(gres, aes(s, fsr.const)) +
  stat_summary(fun.y=fun.y,geom="point",size=2,color="blue") +
  stat_summary(fun.y=fun.y,geom="line",color="blue") +
  #stat_summary(fun.ymax=fun.ymax,fun.ymin=fun.ymin,geom="errorbar",color="blue") +
  scale_x_log10() + facet_wrap(~ns) +
  geom_hline(yintercept=alpha, col="black") +
  xlab("true LFC SD") + ylab("false sign rate") +
  theme(axis.title=element_text(size=16, face="bold"))
plot_grid(g)

# the false sign rate when adapting, at various links btwn
# prior scale and the oracle scale of betas
g <- ggplot(gres, aes(s, fsr, col=scale, shape=scale, group=scale)) +
  stat_summary(fun.y=fun.y,geom="point",size=3) +
  stat_summary(fun.y=fun.y,geom="line") +
  #stat_summary(fun.ymax=fun.ymax,fun.ymin=fun.ymin,geom="errorbar") +
  scale_x_log10() + facet_wrap(~ns) +
  geom_hline(yintercept=alpha, col="black") +
  xlab("true LFC SD") + ylab("false sign rate") +
  theme(axis.title=element_text(size=16, face="bold")) + 
  coord_cartesian(ylim=c(0,0.06))
plot_grid(g)

g <- ggplot(gres, aes(s, fsr, col=scale, shape=scale, group=scale)) +
  stat_summary(fun.y=fun.y,geom="point",size=3) +
  stat_summary(fun.y=fun.y,geom="line") +
  #stat_summary(fun.ymax=fun.ymax,fun.ymin=fun.ymin,geom="errorbar") +
  scale_x_log10() + facet_grid(scale~ns) +
  geom_hline(yintercept=alpha, col="black") +
  xlab("true LFC SD") + ylab("false sign rate") +
  theme(axis.title=element_text(size=16, face="bold")) +
  coord_cartesian(ylim=c(0,0.03))
plot_grid(g)

# the 'error inflation factor': too much shrinkage
# when the distribution of betas is very small, and when
# we have many samples to estimate the betas with precision
g <- ggplot(gres, aes(s, eif, col=scale, shape=scale, group=scale)) +
  stat_summary(fun.y=fun.y,geom="point",size=3) +
  stat_summary(fun.y=fun.y,geom="line") +
  #stat_summary(fun.ymax=fun.ymax,fun.ymin=fun.ymin,geom="errorbar") +
  scale_x_log10() + facet_wrap(~ns) +
  geom_hline(yintercept=1, col="black") +
  xlab("true LFC SD") + ylab("error inflation factor") +
  theme(axis.title=element_text(size=16, face="bold"))
plot_grid(g)
