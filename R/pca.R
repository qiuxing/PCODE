## A wrapper for conducting various PCA and produces the same
## output. As a "bonus", it produces the initial values of X through
## smoothing.
pcafun <- function(y,Ts,K, lambda=0.01, method=c("fpca", "pca", "spca"),
                   center=FALSE, spca.para=2^seq(K)/2){
  method <- match.arg(method)
  mybasis <- create.bspline.basis(range(Ts), length(Ts)+4-2, 4, Ts)
  mypar <- fdPar(mybasis, 2, lambda=lambda)  #under-smooth

  ycurves <- smooth.basis(Ts, y, mypar)[["fd"]]
  pcnames <- paste("PC",1:K,sep="")
  if (method=="fpca"){
    rr <- pca.fd(ycurves, centerfns=center, nharm=K)
    ## as for the centers, fPCA is an odd ball.
    if (center) {
      meancur <- rr[["meanfd"]]
      ## centers <- as.vector(eval.fd(rr[["meanfd"]],Ts))
      centers <- as.vector(Ts, eval.fd(rr[["meanfd"]]))
    } else {
      centers <- rep(0.0,nrow(y))
      meancur <- fd(rep(0.0, mybasis$nbasis), mybasis)
    }; names(centers) <- rownames(y)
    Bhat <- rr[["scores"]]; dimnames(Bhat) <- list(colnames(y),pcnames)
    Bhat.qr <- qr(Bhat)
    ## rotation <- qr.Q(Bhat.qr); dimnames(rotation) <- list(colnames(y),pcnames)
    ## sdev <- abs(diag(qr.R(Bhat.qr))); names(sdev) <- pcnames
    ## xhats.curves are just the harmonics
    xhats.curves <- rr[["harmonics"]]
    ## for fPCA, we evaluate eigen-functions at Ts directly. They are
    ## near, but not exactly, orthogonal.
    ## xhats <- eval.fd(xhats.curves, Ts)
    xhats <- eval.fd(Ts, xhats.curves)
  } else if (method=="pca"){
    rr <- prcomp(t(y), center=center)
    if (center) {
      centers <- rr[["center"]]
    } else {
      centers <- rep(0.0,nrow(y))
    }
    meancur <- smooth.basis(Ts, centers, mypar)[["fd"]]
    xhats <- rr[["rotation"]][,1:K]  #The first K eigenvectors
    Bhat <- rr[["x"]][,1:K]
    ## The xhat curves are represented by smoothed splines of xhats
    xhats.curves <- smooth.basis(Ts, xhats, mypar)[["fd"]]
  } else if (method=="spca"){
    stop("Currently not available.")
  } else {
    stop("Only the following PCA methods are implemented: fpca, pca, spca.")
  }
  ## finally, estimate Binv, the generalized inverse of Bhat. This
  ## matrix is used in between-subject fitting
  ## 09/08/2012. The following shortcut is no good: it has random signs.
  ## Binv <- diag(1/sdev) %*% t(rotation)
  Binv <- solve(t(Bhat) %*% Bhat) %*% t(Bhat)

  return(list(centers=centers, meancur=meancur, Bhat=Bhat, Binv=Binv, xhats=xhats, xhats.curves=xhats.curves, parameters=list(Ts=Ts, K=K, lambda=lambda, method=method,center=center, spca.para=spca.para)))
}

## group.pcafun estimates a centroid for a group of data (xhats).
group.pcafun <- function(Ylist, Ts, K, method=c("fpca", "pca", "spca"), ...){
  method <- match.arg(method)
  pcalist <- lapply(Ylist, pcafun, Ts=Ts, K=K, method=method, ...)
  ## The Graff mean of these PCA results
  pca.mean <- graff.mean(pcalist)
  ## recompute Bhat for each data
  Bhats2 <- lapply(Ylist, reprojection, PCA1=pca.mean, Bhat.only=TRUE)
  Bmean <- Reduce("+", Bhats2)/length(Ylist)
  Bmean.inv <- solve(t(Bmean) %*% Bmean) %*% t(Bmean)
  ## append these two terms in pca.mean. Note that affine projection
  ## does not alter meancur/centers
  pca.mean[["Bhat"]] <- Bmean; pca.mean[["Binv"]] <- Bmean.inv
  return(pca.mean)
}

