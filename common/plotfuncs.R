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

plotOverMean <- function(res, summary="mae", mu.obs = NULL, keeprule = NULL, error = TRUE) {
  if (summary == "mae") {
    sum.f <- function(x) mean(abs(x), na.rm=TRUE)
  } else if (summary == "rmse") {
    sum.f <- function(x) sqrt(mean(x^2, na.rm=TRUE))
  }
  brks <- c(0, 20, 50, 100, 1000, Inf)
  df <- do.call(rbind, lapply(res, function(x) {
	if (is.null(mu.obs)){
		mu.obs <- rowMeans(x$Y)
	} 
    mu.cut <- cut(mu.obs, brks)
	if (is.null(keeprule)) {
		err <- (as.data.frame(x$lfcs)[,-1] - x$lfcs$true)
	} else {
		err <- (as.data.frame(x$lfcs)[,-1] - x$lfcs$true)[keeprule, ]
	}
    tab <- do.call(cbind, lapply(err, function(i) {
      tapply(i, mu.cut, sum.f)
    }))
    df <- melt(tab)
    names(df) <- c("mean", "method", "error")
    df$mean <- factor(df$mean, levels=levels(mu.cut))
    df
  }))
  df <- fixdf(df)
  fun.y <- mean
  fun.ymax <- function(x) mean(x) + sd(x)/sqrt(length(x))
  fun.ymin <- function(x) mean(x) - sd(x)/sqrt(length(x))
  g <- ggplot(df, aes(x=mean, y=error, color=method, group=method)) +
    stat_summary(fun.y=fun.y,geom="point",size=pt.size) +
    stat_summary(fun.y=fun.y,geom="line",size=line.size) +
    xlab("mean of counts") + ylab("mean absolute error")  +
	scale_x_discrete(labels=c("(0,20]" = brks[2], "(20,50]" = brks[3], "(50,100]"=brks[4], 
                               "(100,1e+03]"=brks[5], "(1e+03,Inf]"=brks[6])) 
  g <- g +
    stat_summary(data=subset(df, method=="apeglm"),
                 aes(x=mean, y=error, color=method),
                 fun.y=fun.y, geom="point", size=pt.size) +
    stat_summary(data=subset(df, method=="apeglm"),
                 aes(x=mean, y=error, color=method),
                 fun.y=fun.y, geom="line", size=line.size)
  if (error) {
    g <- g + stat_summary(fun.ymax=fun.ymax, fun.ymin=fun.ymin, geom="errorbar")
  }
  g
}

plotOverLFC <- function(res, summary="mae", keeprule = NULL, error = TRUE) {
  if (summary == "mae") {
    sum.f <- function(x) mean(abs(x), na.rm=TRUE)
  } else if (summary == "rmse") {
    sum.f <- function(x) sqrt(mean(x^2, na.rm=TRUE))
  }
  ps <- c(.75, .9, .95, .975, .99, .995)
  df <- do.call(rbind, lapply(res, function(x) {
    qs <- round(quantile(abs(x$lfcs$true), ps), 2)  
    qmax <- round(max(abs(x$lfcs$true)),2) + .01
    brks <- c(0, qs, qmax)
	if (is.null(keeprule)) {
		beta.cut <- cut(abs(x$lfcs$true), brks)
		err <- (as.data.frame(x$lfcs)[,-1] - x$lfcs$true)
	} else {
		beta.cut <- cut(abs(x$lfcs$true[keeprule]), brks)
		err <- (as.data.frame(x$lfcs)[,-1] - x$lfcs$true)[keeprule, ]
	}
    
    tab <- do.call(cbind, lapply(err, function(i) {
      tapply(i, beta.cut, sum.f)
    }))
    rownames(tab) <- c(ps,1)
    df <- melt(tab)
    names(df) <- c("quantile", "method", "error")
    df
  }))
  df <- fixdf(df)
  df$quantile <- factor(df$quantile)
  levels(df$quantile) <- c(ps,1)
  fun.y <- mean
  fun.ymax <- function(x) mean(x) + sd(x)/sqrt(length(x))
  fun.ymin <- function(x) mean(x) - sd(x)/sqrt(length(x))
  g <- ggplot(df, aes(x=quantile, y=error, color=method, group=method)) +
    stat_summary(fun.y=fun.y,geom="point", size=pt.size) +
    stat_summary(fun.y=fun.y,geom="line", size=line.size) +
    xlab("quantile of true LFC") + ylab("mean absolute error") 
  g <- g +
    stat_summary(data=subset(df, method=="apeglm"),
                 aes(x=quantile, y=error, color=method),
                 fun.y=fun.y, geom="point", size=pt.size) +
    stat_summary(data=subset(df, method=="apeglm"),
                 aes(x=quantile, y=error, color=method),
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
    stat_summary(fun.y=fun.y,geom="point", size=pt.size) +
    stat_summary(fun.y=fun.y,geom="line", size=line.size) +
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

ArrowPlot <- function(lfcs, mu.obs = NULL){
  # bigpar(3, 3)
  par(mfrow=c(2,3), mar=c(2,2,2,1), cex = 0.5, cex.main = 1.6, cex.axis = 1.2, cex.lab = 1.2)
    nms <- names(lfcs[[1]])
	nms[match(c("ashr.d","ashr.l","edgeRPC5"),nms)] <-
    c("ashr-DESeq2","ashr-limma", "edgeR-PC5")
	methds <- c("apeglm","ashr-DESeq2","ashr-limma","DESeq2","edgeR-PC5")
	ids <- match(methds, nms)
	if (is.null(mu.obs)){
		mu.obs <- rowMeans(x$Y)
	} 
	for (i in ids){
		methn <- nms[i]
		lfcmat <- sapply(lfcs, function(x) x[[i]])
		lfcave <- rowMeans(lfcmat)
		plot(mu.obs, lfcave, log="x", xlim=c(1, 2.5e5), ylim=c(-6,6),
			type="n", xlab = "", 
			ylab = "",
			cex.lab=1, main = methn)
		long <- abs(lfcave - lfcs[[1]]$true) > 0.5
		points(mu.obs[!long], lfcave[!long], cex=.5, pch=20, col=rgb(0,0,0,.2))
		arrows(mu.obs[long], lfcs[[1]]$true[long], 
			mu.obs[long], lfcave[long], length=.1, lwd=2,
			col = "blue")
	}
}
