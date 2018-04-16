pt.size <- 2.5
line.size <- 1

fixdf <- function(df) {
  df$method <- factor(df$method, levels=
                 c("apeglm","ashr.d","ashr.l","DESeq2","edgeRPC5","edgeR","limma"))
  lvls <- levels(df$method)
  levels(df$method)[match(c("ashr.d","ashr.l","edgeRPC5"),lvls)] <-
    c("ashr-DESeq2","ashr-limma", "edgeR-PC5")
  df
}

plotLFCs <- function(lfcs, lim) {
  # bigpar(3, 3)
  par(mfrow=c(3,3), mar=c(2,2,2,1), cex = 0.5, cex.main = 1.6, cex.axis = 1.2, cex.lab = 1.2)
  methn <- names(lfcs)
  methn[match(c("ashr.d","ashr.l","edgeRPC5"),methn)] <- c("ashr-DESeq2","ashr-limma", "edgeR-PC5")
  for (i in 2:length(lfcs)) {
    plot(lfcs$true, lfcs[[i]], cex=1, pch=21,
         xlim=c(-lim,lim), ylim=c(-lim,lim),
         main=methn[i], xlab="", ylab="",
         col=rgb(0,0,0,.3), bg=rgb(0,0,0,.1))
    abline(0,1,col=rgb(1,0,0,0.5))
  }
}

plotOverLFC <- function(res, summary="mae", keeprule, error = TRUE) {
  if (summary == "mae") {
    sum.f <- function(x) mean(abs(x), na.rm=TRUE)
  } else if (summary == "rmse") {
    sum.f <- function(x) sqrt(mean(x^2, na.rm=TRUE))
  }
    brks <- c(0,.25,.5,.75, 1, 1.25, Inf)
	beta.cut <- cut(abs(res$lfcs$true[keeprule]), brks)
	err <- (as.data.frame(res$lfcs)[,-1] - res$lfcs$true)[keeprule, ]
    tab <- do.call(cbind, lapply(err, function(i) {
      tapply(i, beta.cut, sum.f)
    }))
    rownames(tab) <- c("0.25", "0.5", "0.75", "1", "1.25", "Inf")
    df <- melt(tab)
    names(df) <- c("breaks", "method", "error")
  df$breaks <- factor(df$breaks)
  levels(df$breaks) <- c("0.25", "0.5", "0.75", "1", "1.25", "Inf")
  df <- fixdf(df)
  fun.y <- mean
  fun.ymax <- function(x) mean(x) + sd(x)/sqrt(length(x))
  fun.ymin <- function(x) mean(x) - sd(x)/sqrt(length(x))
  g <- ggplot(df, aes(x=breaks, y=error, color=method, group=method)) +
    stat_summary(fun.y=fun.y,geom="point",size=pt.size) +
    stat_summary(fun.y=fun.y,geom="line",size=line.size) +
    xlab("bin by true LFC") + ylab("mean absolute error") 
  g <- g +
    stat_summary(data=subset(df, method=="apeglm"),
                 aes(x=breaks, y=error, color=method),
                 fun.y=fun.y, geom="point", size=pt.size) +
    stat_summary(data=subset(df, method=="apeglm"),
                 aes(x=breaks, y=error, color=method),
                 fun.y=fun.y, geom="line", size=line.size)
  if (error) {
    g <- g + stat_summary(fun.ymax=fun.ymax,fun.ymin=fun.ymin,geom="errorbar")
  }
  g
}

plotOverEST <- function(res, summary="mae", keeprule = NULL, error = TRUE) {
  if (summary == "mae") {
    sum.f <- function(x) mean(abs(x), na.rm=TRUE)
  } else if (summary == "rmse") {
    sum.f <- function(x) sqrt(mean(x^2, na.rm=TRUE))
  }
  
  	brks <- c(0,.25,.5,.75, 1, 1.25, Inf)
	err <- (as.data.frame(res$lfcs)[,-1] - res$lfcs$true)[keeprule, ]
	ests <- (as.data.frame(res$lfcs)[,-1])[keeprule, ] 
	bins <- apply(ests, 2, function(x) cut(abs(x), brks))
    tab <- matrix(NA, nrow = 6, ncol = 7)
    for (i in 1:ncol(err)){
		temp <- tapply(err[,i], bins[,i], sum.f)
		tab[1:length(temp),i] <- data.matrix(temp)
	}
    rownames(tab) <- c("0.25", "0.5", "0.75", "1", "1.25", "Inf")
    colnames(tab) <- names(res$lfcs)[2:8]
	df <- melt(tab)
    names(df) <- c("breaks", "method", "error")
  df$breaks <- factor(df$breaks)
  levels(df$breaks) <- c("0.25", "0.5", "0.75", "1", "1.25", "Inf")
  df <- fixdf(df)
  fun.y <- mean
  fun.ymax <- function(x) mean(x) + sd(x)/sqrt(length(x))
  fun.ymin <- function(x) mean(x) - sd(x)/sqrt(length(x))
  g <- ggplot(df, aes(x=breaks, y=error, color=method, group=method)) +
    stat_summary(fun.y=fun.y,geom="point",size=pt.size) +
    stat_summary(fun.y=fun.y,geom="line",size=line.size) +
    xlab("bin by estimated LFC") + ylab("mean absolute error") 
  g <- g +
    stat_summary(data=subset(df, method=="apeglm"),
                 aes(x=breaks, y=error, color=method),
                 fun.y=fun.y, geom="point", size=pt.size) +
    stat_summary(data=subset(df, method=="apeglm"),
                 aes(x=breaks, y=error, color=method),
                 fun.y=fun.y, geom="line", size=line.size)
  if (error) {
    g <- g + stat_summary(fun.ymax=fun.ymax,fun.ymin=fun.ymin,geom="errorbar")
  }
  g
}

catPlot <- function(res, keeprule = NULL, error = TRUE) {
  tops <- c(100,150,200,250,300,350,400)
  df <- do.call(rbind, lapply(res, function(x) {
    tab <- do.call(cbind, lapply(2:length(x$lfcs), function(i) {
      sapply(tops, function(top) {
	    if (is.null(keeprule)) {
		length(intersect(order(-abs(x$lfcs$true))[1:top],
                         order(-abs(x$lfcs[[i]]))[1:top]))/top
		} else {
        length(intersect(order(-abs(x$lfcs$true[keeprule]))[1:top],
                         order(-abs(x$lfcs[[i]][keeprule]))[1:top]))/top
		}
      })
    }))
    rownames(tab) <- tops
    colnames(tab) <- names(x$lfcs)[-1]
    df <- melt(tab)
    names(df) <- c("top", "method", "concordance")
    df
  }))
  df <- fixdf(df)
  fun.y <- mean
  fun.ymax <- function(x) mean(x) + sd(x)/sqrt(length(x))
  fun.ymin <- function(x) mean(x) - sd(x)/sqrt(length(x))
  g <- ggplot(df, aes(x=top, y=concordance, color=method, group=method)) +
    stat_summary(fun.y=fun.y,geom="point",size=pt.size) +
    stat_summary(fun.y=fun.y,geom="line",size=line.size) +
    xlab("top genes by LFC") + ylab("concordance")
  g <- g +
    stat_summary(data=subset(df, method=="apeglm"),
                 aes(x=top, y=concordance, color=method),
                 fun.y=fun.y, geom="point", size=pt.size) +
    stat_summary(data=subset(df, method=="apeglm"),
                 aes(x=top, y=concordance, color=method),
                 fun.y=fun.y, geom="line", size=line.size)
  if (error) {
    g <- g + stat_summary(fun.ymax=fun.ymax,fun.ymin=fun.ymin,geom="errorbar")
  }
  g
}


