## 10/18/2012. Make a package instead of a flat R file.

## wrapper for linear (homogeneous or inhomogeneous) ODE
## solver. const=T/F: whether the equation system is homogeneous or
## inhomogeneous.
ode.fit <- function(Ts, xinit, A, b=NULL){
  if (is.null(b)) b <- rep(0,ncol(A))
  pars <- list(A=A,b=b)
  xderiv <- function(Time, x, pars){
    with(as.list(c(x, pars)), {
      return(list(t(A %*% x)+b))
    })
  }
  out <- ode(y=xinit, parms=pars, times=Ts, func=xderiv)
  as.data.frame(out)
}

## The two-stage backend
.est.2stage <- function(Xt){
  X <- eval.fd(Ts, Xt); X.deriv <- eval.fd(Ts,deriv(Xt))
  Ahat <- t(X.deriv) %*% X %*% solve(t(X) %*% X)
  return(Ahat)
}

## The unweighted PDA backend.
.est.pda <- function(Xt) {
  Sigma.xx <- inprod(Xt, Xt)
  Sigma.xderivx <- inprod(deriv(Xt), Xt)
  Ahat <- Sigma.xderivx %*% solve(Sigma.xx)
  return(Ahat)
}

## This is a full least-square method
.est.FME <- function(Ts, Xt, xhats, const) {
  Xcost <- function(pars, Ts, xinit, const){
    if (const) {
      Ab <- matrix(pars, ncol=length(xinit)+1)
      A <- Ab[-1,-1]; b <- Ab[-1,1]
      out <- ode.fit(Ts, xinit, A, b)
    } else {
      A <- matrix(pars, ncol=length(xinit))
      out <- ode.fit(Ts, xinit, A)
    }
    cost <- modCost(model=out, obs=cbind(time=Ts, xhats))
    return(cost)
  }

  ## use 2-stage method to estimate the initial values
  Ab0 <- .est.2stage(Xt)

  X <- eval.fd(Ts, Xt)
  if (const) {                           #inhomogeneous
    x0 <- X[1,-1]
  } else {                              #homogeneous
    x0 <- X[1,]
  }
  ## Make sure x0 and xhats have the same (non.null names)
  pcnames <- paste("PC",1:length(x0),sep="")
  names(x0) <- pcnames; colnames(xhats) <- pcnames
  
  myfit <- modFit(f=Xcost, p=as.vector(Ab0), Ts=Ts, xinit=x0, const=const)
  Ahat <- matrix(coef(myfit), ncol=ncol(X))
  return(Ahat)
}

## wrapper for low-dim (intrinsic) ODE estimation. Input: Ts; xhats
## and xhats.init estimated from pcafun(); C.init estimated from
## C.init.est(). const=T/F: whether the equation system is homogeneous
## or inhomogeneous.
lowdim.est <- function(Ts, xhats, xhats.curves, method=c("pda", "two.stage", "FME"), lambda=0.01, const=TRUE){
  K <- ncol(xhats); pcnames <- paste("PC",1:K,sep="")
  ## must make sure both xhats and xinit has the same, non-null names
  colnames(xhats) <- pcnames
  mybasis <- xhats.curves[["basis"]]
  mypar <- fdPar(mybasis, 2, lambda=lambda)

  ## now deal with the inhomogeneous case
  if (const) {
    Xt <- fd(cbind("Const"=1, coef(xhats.curves)), mybasis)
  } else {
    Xt <- xhats.curves
  }

  ## Different backends
  method <- match.arg(method)
  if (method=="pda"){
    Ab <- .est.pda(Xt)
  } else if (method=="two.stage"){
    Ab <- .est.pda(Xt)
  } else if (method=="FME") {
    Ab <- .est.FME(Ts, Xt, xhats, const=const)
  } else {
    stop("Only the following low dimensional ODE estimation methods are implemented: pda, two.stage, FME.")
  }

  ## Disentangle A and b.
  if (const) {
    Ahat <- Ab[-1, -1]; bvec <- as.vector(Ab[-1, 1])
  } else {
    Ahat <- Ab; bvec <- rep(0,K)
  }
  dimnames(Ahat) <- list(pcnames, pcnames)
  names(bvec) <- pcnames

  ## Now produce some additional information based on A and b.
  ## don't need the "time" column anymore
  x0 <- t(eval.fd(Ts[1], xhats.curves))
  xhats.fit <- as.matrix(ode.fit(Ts, x0, Ahat, bvec)[,-1])
  rownames(xhats.fit) <- Ts
  ## a spline representation of xhats
  xhats.fit.curves <- smooth.basis(Ts, xhats.fit, mypar)[["fd"]]

  return(list(Ahat=Ahat, bvec=bvec,
              xhats.fit=xhats.fit,
              xhats.fit.curves=xhats.fit.curves,Ts=Ts))
}

## The main function
PCODE <- function(y, Ts, K, lambda=0.01, pca.method=c("fpca", "pca", "spca"), lowdim.method=c("pda", "two.stage", "FME"), center=FALSE, spca.para=2^seq(K)/2, const=TRUE){
  pca.method <- match.arg(pca.method); lowdim.method <- match.arg(lowdim.method)
  ## test if y is a list of subjects or just one subject
  if (is.list(y)) {                     #many subjects
    m <- ncol(y[[1]])
    pca.results <- group.pcafun(y, Ts=Ts, K=K, method=pca.method, center=center, spca.para=spca.para)
  } else {                              #just one subject
    m <- ncol(y)
    pca.results <- pcafun(y, Ts, K=K, lambda=lambda, method=pca.method, center=center, spca.para=spca.para)
  }

  xhats=pca.results[["xhats"]]; xhats.curves <- pca.results[["xhats.curves"]]
  intrinsic.system <- lowdim.est(Ts, xhats=xhats, xhats.curves=xhats.curves,
                                 method=lowdim.method, lambda=lambda, const=const)
  xhats.fit <- intrinsic.system[["xhats.fit"]]
  xhats.fit.curves <- intrinsic.system[["xhats.fit.curves"]]

  pcnames <- paste("PC",1:K,sep=""); colnames(xhats.fit) <- pcnames
  muvec <- matrix(rep(pca.results[["centers"]], m), nrow=length(Ts))
  y.fit <- xhats.fit %*% t(pca.results[["Bhat"]]) + muvec
  mucoef <- matrix(rep(coef(pca.results[["meancur"]]), m), ncol=m)
  mybasis <- xhats.curves[["basis"]]
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
               Ahat=intrinsic.system[["Ahat"]],
               bvec=intrinsic.system[["bvec"]],
               Bhat=pca.results[["Bhat"]], Binv=pca.results[["Binv"]],
               meancur=pca.results[["meancur"]],
               centers=pca.results[["centers"]],
               lambda=lambda,
               pca.results=pca.results,
               intrinsic.system=intrinsic.system))
}

## The main function
PCODE.weighted <- function(y, Ts, K, lambda=0.01, pca.method=c("fpca", "pca", "spca"), lowdim.method=c("pda", "two.stage", "FME"), center=FALSE, spca.para=2^seq(K)/2, const=TRUE){
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
  pcnames <- paste("PC",1:K,sep="")
  dimnames(Ahat) <- list(pcnames, pcnames)
  bvec <- as.vector(diag(1/lambda.root) %*% wbvec); names(bvec) <- pcnames
  xhats.fit <- intrinsic.system[["xhats.fit"]] %*% diag(1/lambda.root)
  colnames(xhats.fit) <- pcnames
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
  Ts <- pcode.fit[["Ts"]]; Ahat <- pcode.fit[["Ahat"]]
  bvec <- pcode.fit[["bvec"]]
  ## No matter what Ts.new is, center0 is the first center estimated
  ## from the original system to match y0.new.
  center0 <- pcode.fit[["centers"]][1]
  xhat0 <- as.vector(pcode.fit[["Binv"]] %*% (y0.new - center0))
  if (is.null(Ts.new)){                 #just use the original Ts
    Ts.new <- Ts; centers <- pcode.fit[["centers"]]
    xhats.fit <- as.matrix(ode.fit(Ts.new, xhat0, Ahat, bvec)[,-1])
  } else {                              #estimate the centers for new Ts
    centers <- eval.fd(Ts.new, pcode.fit[["meancur"]])
    ## Due to the design of ode(), we need to append Ts1 to Ts.new
    ## before fitting the curve, and remove it from the results.
    xhats.fit <- as.matrix(ode.fit(c(Ts[1], Ts.new), xhat0, Ahat, bvec)[-1,-1])
  }
  y.fit <- sweep(as.matrix(xhats.fit) %*% t(pcode.fit[["Bhat"]]), 1, -centers)
  return (y.fit)
}

## individual cross-validation
CV <- function(Y, Ts, K, center=FALSE, const=FALSE, ...){
  J <- length(Ts)
  ## Due to the fact that it is hard to extrapolate, we will drop the
  ## first time point for now.
  trainsys <- foreach(j=2:(J-1)) %dopar%{
    ## Ts=Ts; K=K
    PCODE(Y[-j,], Ts=Ts[-j], K=K, center=center, const=const, ...)
  }
  y.fits <- t(sapply(2:(J-1), function(j) predict.pcode1(trainsys[[j-1]], Y[1,], Ts.new=Ts[j])))
  rss <- sum((Y[2:(J-1),] - y.fits)^2)/ncol(Y)
  return(rss)
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

  y2 <- as.matrix(y[,genes]); y2.fit.curves <- pcode.result[["y.fit.curves"]][genes]
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

## another wrapper for group CV results.
plot.cv <- function(Ylist, cv.results, true.y.curves=NULL, genes=NULL, subjects=NULL, plot.orig.curves=TRUE, ...){
  N <- length(Ylist); ngenes <- ncol(Ylist[[1]]); Ts <- cv.results[["Ts"]]
  ## If not specified, use the first three subjects.
  if (is.null(subjects)) subjects <- seq(min(N,3))
  ## By default, plot the first 4 genes per each subjects. Otherwise,
  ## plot genes listed in the genes.
  if (is.null(genes)) genes=seq(min(ngenes,4))

  for (ss in subjects){
    y2 <- as.matrix(Ylist[[ss]][,genes])
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
