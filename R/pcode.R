## 9/13/2012.  The sign of y.fit is still random.

## useful functions for high-dim ODE project
projected.fnorm <- function(V, W){
  ## projection Frobenius norm.
  Qv <- qr.Q(qr(V)); Qw <- qr.Q(qr(W))
  diffmat <- Qv %*% t(Qv) - Qw %*% t(Qw)
  return(norm(diffmat, "F")/sqrt(2))
}

grassmann.mean <- function(Ws){
  ## A simple way of computing centroid of points on Gr(m,n)
  N <- length(Ws)                       #sample size
  K <- ncol(Ws[[1]])                    #dim of eigen-space
  Qws <- sapply(Ws, function(W) qr.Q(qr(W)))

  mean.mat2 <- apply(sapply(Qws, function(Qw) W%*%t(Qw), simplify="array"), c(1,2), mean)
  eigen(apply(mean.mat2, c(1,2), mean))[["vectors"]][,1:K]
}



## A wrapper for conducting various PCA and produces the same
## output. As a "bonus", it produces the initial values of X through
## smoothing.
pcafun <- function(y,Ts,K, lambda=0.01, method=c("fpca", "pca", "spca"),
                   center=TRUE, spca.para=2^seq(K)/2){
  method <- match.arg(method)
  mybasis <- create.bspline.basis(range(Ts), length(Ts)+4-2, 4, Ts)
  mypar <- fdPar(mybasis, 2, lambda=lambda)  #under-smooth

  ycurves <- smooth.basis(Ts, y, mypar)[["fd"]]
  pcnames <- paste("PC",1:K,sep="")
  if (method=="fpca"){
    rr <- pca.fd(ycurves, centerfns=center, nharm=K)
    ## as for the centers, fPCA is an odd ball.
    if (center) {
      centers <- as.vector(eval.fd(rr[["meanfd"]],Ts))
    } else {
      centers <- rep(0,nrow(y))
    }; names(centers) <- rownames(y)
    Bhat <- rr[["scores"]]; dimnames(Bhat) <- list(colnames(y),pcnames)
    Bhat.qr <- qr(Bhat)
    ## rotation <- qr.Q(Bhat.qr); dimnames(rotation) <- list(colnames(y),pcnames)
    ## sdev <- abs(diag(qr.R(Bhat.qr))); names(sdev) <- pcnames
    ## xhats.curves are just the harmonics
    xhats.curves <- rr[["harmonics"]]
    ## for fPCA, we evaluate eigen-functions at Ts directly. They are
    ## near, but not exactly, orthogonal.
    xhats <- eval.fd(xhats.curves, Ts)
  } else if (method=="pca"){
    rr <- prcomp(y, center=center)
    if (center) {
      centers <- rr[["center"]]
    } else {
      centers <- rep(0,ncol(y))
    }
    rotation <- rr[["rotation"]][,1:K]
    sdev <- rr[["sdev"]][1:K]           #square root of eigenvalues
    names(sdev) <- pcnames
    Bhat <- rotation %*% diag(sdev)     #y = Bhat %*% xhats
    colnames(Bhat) <- pcnames
    ## make sure xhats, the eigenvectors, are orthonormal.
    xhats <- rr[["x"]][,1:K] %*% diag(1/sdev)
    colnames(xhats) <- pcnames
    ## The xhat curves are represented by smoothed splines of xhats
    xhats.curves <- smooth.basis(Ts, xhats, mypar)[["fd"]]
  } else if (method=="spca"){
    if (center){
      centers <- colMeans(y)
      yy <- t(scale(y, scale=FALSE)) %*% scale(y, scale=FALSE)
    } else {
      centers <- rep(0,ncol(y))
      yy <- t(y)%*%y
    }
    rr <- spca(yy, K=K, type="Gram", sparse="penalty", para=spca.para)
    rotation <- rr[["loadings"]]
    ## eigenvectors are only approximately orthonormal
    sdev <- sqrt(rr[["pev"]]*rr[["var.all"]])
    names(sdev) <- pcnames
    xhats <- y %*% rotation %*% diag(1/sdev)
    colnames(xhats) <- pcnames
    Bhat <- rotation %*% diag(sdev)
    colnames(Bhat) <- pcnames
    xhats.curves <- smooth.basis(Ts, xhats, mypar)[["fd"]]
  } else {
    stop("Only the following PCA methods are implemented: fpca, pca, spca.")
  }
  ## finally, estimate Binv, the generalized inverse of Bhat. This
  ## matrix is used in between-subject fitting
  ## 09/08/2012. The following shortcut is no good: it has random signs.
  ## Binv <- diag(1/sdev) %*% t(rotation)
  Binv <- solve(t(Bhat) %*% Bhat) %*% t(Bhat)

  return(list(centers=centers, Bhat=Bhat, Binv=Binv, Ts=Ts, xhats=xhats, xhats.curves=xhats.curves, parameters=list(K=K, lambda=lambda, method=method,center=center, spca.para=spca.para)))
}

## This function takes a linear subspace (results from a PCA) and
## apply it to another set of data
reprojection <- function(Y2, PCA1){
  method <- PCA1[["parameters"]][["method"]]
  lambda <- PCA1[["parameters"]][["lambda"]]
  Ts <- PCA1[["Ts"]]; K <- PCA1[["parameters"]][["K"]]
  Xt <- PCA1[["xhats.curves"]]
  pcnames <- paste("PC",1:K,sep="")
  if (method=="fpca"){
    bs1 <- Xt[["basis"]]
    mypar <- fdPar(bs1, 2, lambda=lambda)  #under-smooth
    ycurves <- smooth.basis(Ts, Y2, mypar)[["fd"]]
    Bhat <- inprod(ycurves, Xt)
  } else if (method=="pca"){
    stop("not implemented now.")
  } else if (method=="spca"){
    stop("not implemented now.")
  } else {
    stop("Only the following PCA methods are implemented: fpca, pca, spca.")
  }
  rownames(Bhat) <- colnames(Y2); colnames(Bhat) <- pcnames
  return(Bhat)
}


## group.pcafun estimates a centroid for a group of data (xhats).
group.pcafun <- function(Ylist, Ts, K, ...){
  n <- length(Ylist); genenames <- colnames(Ylist[[1]])
  pcnames <- paste("PC",1:K,sep="")
  pcalist <- lapply(Ylist, pcafun, Ts=Ts, K=K, method=method, ...)
  if (method=="fpca"){
    ## uses eigen-functions (curves) as the tool to compute the centroid
    Xts <- lapply(pcalist, function(pp) coef(pp[["xhats.curves"]]))
    Qxs2 <- lapply(Xts, function(X) {Qx=qr.Q(qr(X)); Qx%*%t(Qx)})
    Phat <- Reduce("+", Qxs2)/n
    Ehat <- eigen(Phat)[["vectors"]][,1:K]

    bs1 <- pcalist[[1]][["xhats.curves"]][["basis"]]
    Beta <- fd(diag(bs1[["nbasis"]]), bs1)
    SigmaBeta <- inprod(Beta, Beta); ee <- eigen(SigmaBeta)
    Tbeta <- ee[["vectors"]]; Lambda.negroot <- diag(1/sqrt(ee[["values"]]))
    xmean.curve <- fd(Tbeta %*% Lambda.negroot %*% Ehat, bs1)
    xmean <- eval.fd(xmean.curve, Ts)
  } else {
    stop("not implemented yet.")
  }
  ## wrap it up in a new pca object for use with reprojection
  pca.mean <- pcalist[[1]]            #the template
  pca.mean[["xhats"]] <-xmean
  pca.mean[["xhats.curves"]] <- xmean.curve
  ## recompute Bhat for each data
  Bhats2 <- lapply(Ylist, reprojection, PCA1=pca.mean)
  Bmean <- Reduce("+", Bhats2)/n
  Bmean.inv <- solve(t(Bmean) %*% Bmean) %*% t(Bmean)
  Bmean.qr <- qr(Bmean)
  ## update other terms in pca.mean
  pca.mean[["Bhat"]] <- Bmean; pca.mean[["Binv"]] <- Bmean.inv
  pca.mean[["centers"]] <- Reduce("+", lapply(pcalist,function(pp) pp[["centers"]]))/n
  return(pca.mean)
}

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
lowdim.est <- function(Ts, xhats, xhats.curves, method=c("FME", "pda"), const=TRUE){
  K <- ncol(xhats)
  Xcost <- function(pars, xinit){
    out <- Xfun(pars, Ts, xinit=xinit, const=const)
    cost <- modCost(model=out, obs=cbind(time=Ts, xhats))
    return(cost)
  }

  if (const) {                           #inhomogeneous
    ## z is the smoothed xhats plus a constant 1 (for intercept terms)
    z <- cbind(Const=rep(1,length(Ts)),eval.fd(xhats.curves,Ts))
    ## estimates of the derivatives
    z.deriv <- cbind(Const=rep(0,length(Ts)), eval.fd(deriv(xhats.curves), Ts))
    ## then use simple linear regression to obtain an approximate C0.
    ## remember CC is the t(CC) because regression matrix and ODE matrix
    ## differ by one transpose.
    CC0 <- t(solve(t(z) %*% z) %*% t(z) %*% z.deriv)
    ## pars0[1:K] (the first column) are the constant terms. Others form
    ## the Ahat
    pars0 <- as.vector(CC0[-1,])
    myfit <- modFit(f=Xcost, p=pars0, xinit=z[1,-1])
    CC <- matrix(coef(myfit), ncol=K+1)
    Ahat <- CC[,-1]; bvec <- CC[,1]
    ## don't need the "time" column anymore
    xhats.fit <- Xfun(pars=coef(myfit), Ts, xinit=z[1,-1], const=const)[,-1]
  } else {                              #homogeneous
    z <- eval.fd(xhats.curves,Ts)
    z.deriv <- eval.fd(deriv(xhats.curves), Ts)
    CC0 <- t(solve(t(z) %*% z) %*% t(z) %*% z.deriv)
    pars0 <- as.vector(CC0)
    myfit <- modFit(f=Xcost, p=pars0, xinit=z[1,])
    Ahat <- matrix(coef(myfit), ncol=K); bvec <- rep(0,K)
    xhats.fit <- Xfun(pars=coef(myfit), Ts, xinit=z[1,], const=const)[,-1]
  }
  pcnames <- paste("PC",1:K,sep="")
  dimnames(Ahat) <- list(pcnames, pcnames)
  names(bvec) <- pcnames
  rownames(xhats.fit) <- Ts
  return(list(Ahat=Ahat, bvec=bvec, deviance=deviance(myfit), iterations=myfit[["iterations"]],xhats.fit=xhats.fit, Times=Ts))
}

## The main function
PCODE <- function(y, Ts, K, lambda=0.01, pca.method=c("fpca", "pca", "spca"), lowdim.method=c("FME","pdf"), center=TRUE, spca.para=2^seq(K)/2, const=TRUE){
  pca.method <- match.arg(pca.method)
  lowdim.method <- match.arg(lowdim.method)
  pca.results <- pcafun(y, Ts, K=K, lambda=lambda, method=pca.method,
                        center=center, spca.para=spca.para)
  intrinsic.system <- lowdim.est(Ts, xhats=pca.results[["xhats"]], xhats.curves=pca.results[["xhats.curves"]], method=lowdim.method, const=const)
  xhats.fit <- intrinsic.system[["xhats.fit"]]

  y.fit0 <- as.matrix(xhats.fit) %*% t(pca.results[["Bhat"]])

  ## "centering" refers to different meanings in fPCA and other PCA
  ## methods.
  if (center) {
    if (pca.method=="fpca"){
      ## row centers
      rcenters <- pca.results[["centers"]]
      y.fit <- y.fit0 + matrix(rep(rcenters,ncol(y)), nrow=nrow(y))
    } else {
      ## column centers
      ccenters <- pca.results[["centers"]]
      y.fit <- y.fit0 + matrix(rep(ccenters,nrow(y)), nrow=nrow(y), byrow=T)
    }
  } else {                              #centers are assumed to be zeros
    y.fit <- y.fit0
  }

  return (list(Times=Ts, xhats.fit=xhats.fit, y.fit=y.fit, residuals=y-y.fit,
               Ahat=intrinsic.system[["Ahat"]],
               bvec=intrinsic.system[["bvec"]],
               Bhat=pca.results[["Bhat"]], Binv=pca.results[["Binv"]],
               pca.results=pca.results,
               intrinsic.system=intrinsic.system))
}

## between-subject fitting.  Note that y0.new must be a column vector.
predict.pcode1 <- function(pcode.fit, y0.new, Ts.new="same"){
  ## This function only works with un-centered version as of 09/08/2012.
  Ts <- pcode.fit[["Times"]];
  if (Ts.new=="same"){
    Ts.new <- Ts
  }
  xhat0 <- as.vector(pcode.fit[["Binv"]] %*% y0.new)
  xhats.fit <- Xfun(pars=as.vector(pcode.fit[["Ahat"]]), Ts.new, xinit=xhat0, const=FALSE)[,-1]
  y.fit <- as.matrix(xhats.fit) %*% t(pcode.fit[["Bhat"]])
  return (y.fit)
}

## between subjects fitting, group version.
## predict.pcdode <- ...


## ## wrapper for plotting
## plot.pcode <- function(pcode.result){

## }


## karcher.mean <- function(Ws) {
##   ## the karcher mean/centroid on Gr(m,n)
##   N <- length(Ws)                       #sample size
##   objfun <- function(V){
##     sum(sapply(1:N, function(i) projected.fnorm(V, Ws[[i]])^2))
##   }
##   ...
## }


## align.grassmann <- function(V, W){

## }


## geodesic.norm <- function(V, W){
##   ## The geodesic distance on SOn
## }
