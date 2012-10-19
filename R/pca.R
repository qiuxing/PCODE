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

