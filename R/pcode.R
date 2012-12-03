## 10/18/2012. Make a package instead of a flat R file.

## wrapper for linear (homogeneous or inhomogeneous) ODE
## solver. const=T/F: whether the equation system is homogeneous or
## inhomogeneous.

Xfun <- function(pars, Ts, xinit, const=TRUE){
  K <- length(xinit)
  if (const) {                          #inhomogeneous
    xderiv <- function(t, x, pars){
      ## pars is the vector form of both matrix A and constant vector
      ## b.
      CC <- matrix(pars, ncol=K+1)
      Ahat <- CC[,-1]; bvec <- CC[,1]
      list(t(Ahat %*% x)+bvec)
    }} else {                             #homogeneous
      xderiv <- function(t, x, pars){
        Ahat <- matrix(pars, ncol=K)
        list(t(Ahat %*% x))
      }
    }
  out <- ode(y=xinit, parms=pars, times=Ts, func=xderiv)
  as.data.frame(out)
}

## wrapper for low-dim (intrinsic) ODE estimation. Input: Ts; xhats
## and xhats.init estimated from pcafun(); C.init estimated from
## C.init.est(). const=T/F: whether the equation system is homogeneous
## or inhomogeneous.
lowdim.est <- function(Ts, xhats, xhats.curves, method=c("FME", "pda"), lambda=0.01, const=TRUE){
  K <- ncol(xhats); pcnames <- paste("PC",1:K,sep="")
  ## must make sure both xhats and xinit has the same, non-null names
  colnames(xhats) <- pcnames
  mybasis <- xhats.curves[["basis"]]
  mypar <- fdPar(mybasis, 2, lambda=lambda)

  Xcost <- function(pars, xinit){
    out <- Xfun(pars, Ts, xinit=xinit, const=const)
    cost <- modCost(model=out, obs=cbind(time=Ts, xhats))
    return(cost)
  }

  if (const) {                           #inhomogeneous
    ## z is the smoothed xhats plus a constant 1 (for intercept terms)
    z <- cbind(Const=rep(1,length(Ts)),eval.fd(Ts,xhats.curves))
    ## estimates of the derivatives
    z.deriv <- cbind(Const=rep(0,length(Ts)), eval.fd(Ts,deriv(xhats.curves)))
    ## then use simple linear regression to obtain an approximate C0.
    ## remember CC is the t(CC) because regression matrix and ODE matrix
    ## differ by one transpose.
    CC0 <- t(solve(t(z) %*% z) %*% t(z) %*% z.deriv)
    ## pars0[1:K] (the first column) are the constant terms. Others form
    ## the Ahat
    pars0 <- as.vector(CC0[-1,])
    ## must make sure both xhats and xinit has the same, non-null names
    x0 <- z[1,-1]; names(x0) <- colnames(xhats)
    myfit <- modFit(f=Xcost, p=pars0, xinit=x0)
    CC <- matrix(coef(myfit), ncol=K+1)
    Ahat <- CC[,-1]; bvec <- CC[,1]
    ## don't need the "time" column anymore
    xhats.fit <- as.matrix(Xfun(pars=coef(myfit), Ts, xinit=x0, const=const)[,-1])
  } else {                              #homogeneous
    z <- eval.fd(Ts,xhats.curves)
    z.deriv <- eval.fd(Ts,deriv(xhats.curves))
    CC0 <- t(solve(t(z) %*% z) %*% t(z) %*% z.deriv)
    pars0 <- as.vector(CC0)
    ## must make sure both xhats and xinit has the same, non-null names
    x0 <- z[1,]; names(x0) <- colnames(xhats)
    myfit <- modFit(f=Xcost, p=pars0, xinit=x0)

    Ahat <- matrix(coef(myfit), ncol=K); bvec <- rep(0,K)
    xhats.fit <- as.matrix(Xfun(pars=coef(myfit), Ts, xinit=x0, const=const)[,-1])
  }
  dimnames(Ahat) <- list(pcnames, pcnames)
  names(bvec) <- pcnames
  rownames(xhats.fit) <- Ts
  ## a spline representation of xhats
  xhats.fit.curves <- smooth.basis(Ts, xhats.fit, mypar)[["fd"]]

  return(list(Ahat=Ahat, bvec=bvec,
              deviance=deviance(myfit),
              iterations=myfit[["iterations"]],
              xhats.fit=xhats.fit,
              xhats.fit.curves=xhats.fit.curves,Ts=Ts))
}

## The main function
PCODE <- function(y, Ts, K, lambda=0.01, pca.method=c("fpca", "pca", "spca"), lowdim.method=c("FME","pda"), center=FALSE, spca.para=2^seq(K)/2, const=TRUE){
  pca.method <- match.arg(pca.method); lowdim.method <- match.arg(lowdim.method)
  ## test if y is a list of subjects or just one subject
  if (is.list(y)) {                     #many subjects
    m <- ncol(y[[1]])
    pca.results <- group.pcafun(y, Ts=Ts, K=K, method=pca.method, center=center, spca.para=spca.para)
  } else {                              #just one subject
    m <- ncol(y)
    pca.results <- pcafun(y, Ts, K=K, lambda=lambda, method=pca.method, center=center, spca.para=spca.para)
  }
  ## now use square root of varprop to weight xhats, xhats.curves
  lambda.root <- sqrt(pca.results[["varprop"]])
  wxhats <- pca.results[["xhats"]] %*% diag(lambda.root)
  xhats.curves <- pca.results[["xhats.curves"]]
  mybasis <- xhats.curves[["basis"]]
  ## below starts the fitting of the intrinsic system.  Note that we
  ## use sqrt(varprop) to weight xhats and xhats.curves.
  wxhats.curves <- fd(coef(xhats.curves) %*% diag(lambda.root), mybasis)
  intrinsic.system <- lowdim.est(Ts, xhats=wxhats, xhats.curves=wxhats.curves, method=lowdim.method, lambda=lambda, const=const)
  wAhat <- intrinsic.system[["Ahat"]]; wbvec <- intrinsic.system[["bvec"]]
  Ahat <- diag(1/lambda.root) %*% wAhat %*% diag(lambda.root)
  bvec <- diag(1/lambda.root) %*% wbvec
  xhats.fit <- intrinsic.system[["xhats.fit"]] %*% diag(1/lambda.root)
  pcnames <- paste("PC",1:K,sep=""); colnames(xhats.fit) <- pcnames
  ## Done fitting intrinsic system. Now translate these weighted
  ## curves back.
  wxhats.fit.curves <- intrinsic.system[["xhats.fit.curves"]]
  xcoefs <- coef(wxhats.fit.curves) %*% diag(1/lambda.root)
  colnames(xcoefs) <- pcnames
  xhats.fit.curves <- fd(xcoefs,mybasis)
  muvec <- matrix(rep(pca.results[["centers"]], m), nrow=length(Ts))
  y.fit <- xhats.fit %*% t(pca.results[["Bhat"]]) + muvec
  mucoef <- matrix(rep(coef(pca.results[["meancur"]]), m), ncol=m)
  y.fit.curves <- fd(coef(xhats.fit.curves) %*% t(pca.results[["Bhat"]]) + mucoef, mybasis)

  if (is.list(y)) {                     #many subjects
    rss <- mean(sapply(y, function(yn) sum((yn-y.fit)^2)))/m
  } else {                              #just one subject
    rss <- sum((y-y.fit)^2/m)
  }
  return (list(Ts=Ts, xhats.fit=xhats.fit,
               xhats.fit.curves=xhats.fit.curves,
               y.fit=y.fit, y.fit.curves=y.fit.curves,
               rss=rss, varprop=pca.results[["varprop"]],
               Ahat=Ahat, bvec=bvec,
               Bhat=pca.results[["Bhat"]], Binv=pca.results[["Binv"]],
               meancur=pca.results[["meancur"]],
               centers=pca.results[["centers"]],
               lambda=lambda,
               pca.results=pca.results,
               weighted.intrinsic.system=intrinsic.system))
}

## between-subject fitting.  Note that y0.new must be a column vector.  As of ver 0.02, this prediction function does not work well with const != 0 case.
predict.pcode1 <- function(pcode.fit, y0.new, Ts.new=NULL){
  Ts <- pcode.fit[["Ts"]]; centers <- pcode.fit[["centers"]]
  if (is.null(Ts.new)){
    Ts.new <- Ts
  }
  xhat0 <- as.vector(pcode.fit[["Binv"]] %*% (y0.new - centers[1]))
  CC <- cbind(pcode.fit[["bvec"]], pcode.fit[["Ahat"]])
  pars1 <- as.vector(CC)
  xhats.fit <- Xfun(pars=pars1, Ts.new, xinit=xhat0, const=TRUE)[,-1]
  y.fit <- sweep(as.matrix(xhats.fit) %*% t(pcode.fit[["Bhat"]]), 1, -centers)
  return (y.fit)
}

## wrapper for between subject cross-validation. As of ver 0.02, this function does not work with const=TRUE case.
CV.group <- function(Ylist, Ts, K, center=FALSE, const=FALSE, ...){
  N <- length(Ylist)
  trainsys <- foreach(n=1:length(Ylist)) %dopar%{
    Ylist.n <- Ylist[-n]
    Ts=Ts; K=K
    PCODE(Ylist.n, Ts=Ts, K=K, center=center, const=const, ...)
  }
  y.fit.list <- lapply(1:N, function(n) predict.pcode1(trainsys[[n]], Ylist[[n]][1,]))
  rss <- sapply(1:N, function(n) sum((Ylist[[n]]-y.fit.list[[n]])^2))/ncol(Ylist[[1]])
  ## ## To view RSS from anthother angle, the explained L^2 norm
  ## squared. Scratch this. Not an orthogonal decomposition so no
  ## var-decomposition.  totalL2 <- sum(sapply(Ylist, function(y)
  ## sum(y^2))) L2.fit <- sum(sapply(y.fit.list, function(y)
  ## sum(y^2))) L2prop <- L2.fit/totalL2 produce the CV fitted curves
  ## as well
  mybasis <- trainsys[[1]][["y.fit.curves"]][["basis"]]
  lambda <- trainsys[[1]][["lambda"]]
  mypar <- fdPar(mybasis, 2, lambda=lambda)
  y.fit.curves.list <- lapply(y.fit.list, function(y) smooth.basis(Ts, y, mypar)[["fd"]])
  ## Assign proper subject names for convenience
  names(y.fit.list) <- names(Ylist); names(y.fit.curves.list) <- names(Ylist)
  names(rss) <- names(Ylist)
  return(list(trainsys=trainsys, y.fit.list=y.fit.list, y.fit.curves.list=y.fit.curves.list, rss=rss, lambda=lambda, Ts=Ts))
}


## ## wrapper for plotting
plot.pcode <- function(y, pcode.result, true.y.curves=NULL, genes=NULL, plot.orig.curves=TRUE, ...){
  ngenes <- ncol(y); Ts <- pcode.result[["Ts"]]
  ## By default, plot the first 12 genes. Otherwise, plot genes listed
  ## in the genes.
  if (is.null(genes)) genes=seq(min(ngenes,12))

  y2 <- y[,genes]; y2.fit.curves <- pcode.result[["y.fit.curves"]][genes]
  if (is.null(colnames(y2))) {
    ynames <- paste("y", genes, sep="")
  } else {
    ynames <- colnames(y2)
  }
  ## Smooth the original data and plot the smoothed curves as well.
  if (plot.orig.curves){
    lambda <- pcode.result[["lambda"]]
    mybasis <- y2.fit.curves[["basis"]]
    mypar <- fdPar(mybasis, 2, lambda=lambda)
    y2.fit.orig <- smooth.basis(Ts, y2, mypar)[["fd"]]
  }  

  if (!is.null(true.y.curves)) true.y2.curves <- true.y.curves[genes]

  ## Now plot these curves
  for (i in 1:length(genes)){
    plot(Ts, y2[,i], ylab=ynames[i], xlab="Time", ...)
    lines(y2.fit.curves[i])
    if (plot.orig.curves) lines(y2.fit.orig[i], col="grey")
    if (!is.null(true.y.curves)) lines(y2.fit.orig[i], col="blue")
  }
}

## another wrapper for CV results.
plot.cv <- function(Ylist, cv.results, true.y.curves=NULL, genes=NULL, subjects=NULL, plot.orig.curves=TRUE, ...){
  N <- length(Ylist); ngenes <- ncol(Ylist[[1]]); Ts <- cv.results[["Ts"]]
  ## If not specified, use the first three subjects.
  if (is.null(subjects)) subjects <- seq(min(N,3))
  ## By default, plot the first 4 genes per each subjects. Otherwise,
  ## plot genes listed in the genes.
  if (is.null(genes)) genes=seq(min(ngenes,4))

  for (ss in subjects){
    y2 <- Ylist[[ss]][,genes]
    y2.fit.curves <- cv.results[["y.fit.curves.list"]][[ss]][genes]
    if (is.null(colnames(y2))) {
      ynames <- paste("y", genes, sep="")
    } else {
      ynames <- colnames(y2)
    }
    ## Smooth the original data and plot the smoothed curves as well.
    if (plot.orig.curves){
      lambda <- cv.results[["lambda"]]
      mybasis <- y2.fit.curves[["basis"]]
      mypar <- fdPar(mybasis, 2, lambda=lambda)
      y2.fit.orig <- smooth.basis(Ts, y2, mypar)[["fd"]]
    }
    
    if (!is.null(true.y.curves)) true.y2.curves <- true.y.curves[[ss]][genes]

    ## Now plot these curves
    for (i in 1:length(genes)){
      plot(Ts, y2[,i], ylab=ynames[i], xlab="Time", main=ss, ...)
      lines(y2.fit.curves[i])
      if (plot.orig.curves) lines(y2.fit.orig[i], col="grey")
      if (!is.null(true.y.curves)) lines(y2.fit.orig[i], col="blue")
    }
  }
}
