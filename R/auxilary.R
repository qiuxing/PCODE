## useful functions for high-dim ODE project
projected.fnorm <- function(V, W){
  ## projection Frobenius norm.
  Qv <- qr.Q(qr(V)); Qw <- qr.Q(qr(W))
  diffmat <- Qv %*% t(Qv) - Qw %*% t(Qw)
  return(norm(diffmat, "F")/sqrt(2))
}

## the procrustes mean on Graff(K,m). Needed by group.pcafun().
graff.mean <- function(pcalist){
  ## given a list of PCA or Graff objects (typically returned by
  ## pcafun()), returns a new PCA mean.
  n <- length(pcalist); parameters <- pcalist[[1]][["parameters"]]
  K <- parameters[["K"]]; method=parameters[["method"]]
  lambda <- parameters[["lambda"]]
  mybasis <- pcalist[[1]][["xhats.curves"]][["basis"]]
  mypar <- fdPar(mybasis, 2, lambda=lambda)  #under-smooth

  centers <- lapply(pcalist, function(pp) pp[["centers"]])
  Xs <- lapply(pcalist, function(pp) pp[["xhats"]])
  meancurs.list <- lapply(pcalist, function(pp) pp[["meancur"]])
  meancurs.coefs <- sapply(meancurs.list, coef)
  ## use the vector form for inprod
  meancurs <- fd(meancurs.coefs, mybasis)
  ## Xts: the coefficients of xhats.curves
  Xts <- lapply(pcalist, function(pp) coef(pp[["xhats.curves"]]))
  if (method=="pca"){
    ## Xs returned by pcafun() with methods=(pca|spca) should be
    ## orthonormal; we re-orthonormalize them just to be on the right
    ## track.
    Qxs <- lapply(Xs, function(X) qr.Q(qr(X)))
    Qxs2 <- lapply(Qxs, function(Qx) {Qx%*%t(Qx)})
    mutilde <- lapply(1:n, function(i) centers[[i]]-Qxs2[[i]] %*% centers[[i]])
    Phat <- Reduce("+", Qxs2)/n
    mean.centers <- Reduce("+", mutilde)/n
    mean.xhats <- eigen(Phat)[["vectors"]][,1:K]
    mean.meancur <- smooth.basis(Ts, mean.centers, mypar)[["fd"]]
    mean.xhats.curves <- smooth.basis(Ts, mean.xhats, mypar)[["fd"]]
  } else if (method=="fpca") {
    Beta <- fd(diag(mybasis[["nbasis"]]), mybasis)
    SigmaBeta <- inprod(Beta, Beta); ee <- eigen(SigmaBeta)
    Tbeta <- ee[["vectors"]]
    Lambda.root <- diag(sqrt(ee[["values"]]))
    Lambda.negroot <- diag(1/sqrt(ee[["values"]]))

    ## Qx are the Q_{\tilde{x}} in the manuscript
    Qxs <- lapply(Xts, function(X) qr.Q(qr(Lambda.root %*% t(Tbeta) %*% X)))
    Qxs2 <- lapply(Qxs, function(Qx) {Qx%*%t(Qx)})
                 
    Phat <- Reduce("+", Qxs2)/n
    Ehat <- eigen(Phat)[["vectors"]][,1:K]

    mean.xhats.curves <- fd(Tbeta %*% Lambda.negroot %*% Ehat, mybasis)
    meancur.proj.coefs <- coef(mean.xhats.curves) %*% inprod(mean.xhats.curves, meancurs)
    mean.meancur <- mean(fd(meancurs.coefs - meancur.proj.coefs, mybasis))
    mean.xhats <- eval.fd(mean.xhats.curves, Ts)
    mean.centers <- eval.fd(mean.meancur, Ts)
  } else if  (method=="spca") {
    stop("Currently not available.")
  } else {
    stop("Only the following PCA methods are implemented: fpca, pca, spca.")
  }
  return(list(centers=mean.centers, meancur=mean.meancur,
              xhats=mean.xhats, xhats.curves=mean.xhats.curves,
              parameters=parameters))
}


## This function takes a linear subspace (results from a PCA) and
## apply it to another set of data
reprojection <- function(Y2, PCA1){
  params <- PCA1[["parameters"]]; method <- params[["method"]]
  Ts <- params[["Ts"]]; K <- params[["K"]]
  pcnames <- paste("PC",1:K,sep="")
  if (method=="fpca"){
    Xt <- PCA1[["xhats.curves"]];  bs1 <- Xt[["basis"]]
    mypar <- fdPar(bs1, 2, lambda=params[["lambda"]])
    ycurves <- smooth.basis(Ts, Y2, mypar)[["fd"]]
    ## Since Xt may only be approx orthonormal, I use the following
    ## safer formula to compute Bhat
    Bhat <- t(solve(inprod(Xt, Xt)) %*% inprod(Xt, ycurves))
  } else if (method=="pca"){
    X <- PCA1[["xhats"]]
    Bhat <- t(solve(t(X) %*% X) %*% t(X) %*% Y2)
  } else if (method=="spca"){
    X <- PCA1[["xhats"]]
    Bhat <- t(solve(t(X) %*% X) %*% t(X) %*% Y2)
  } else {
    stop("Only the following PCA methods are implemented: fpca, pca, spca.")
  }
  rownames(Bhat) <- colnames(Y2); colnames(Bhat) <- pcnames
  return(Bhat)
}

## karcher.mean <- function(Ws) {
##   ## the karcher mean/centroid on Gr(m,n)
##   N <- length(Ws)                       #sample size
##   objfun <- function(V){
##     sum(sapply(1:N, function(i) projected.fnorm(V, Ws[[i]])^2))
##   }
##   ...
## }


## geodesic.norm <- function(V, W){
##   ## The geodesic distance on SOn
## }
