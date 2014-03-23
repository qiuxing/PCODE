## A Ridge-like matrix inverse for symmetric matrices
ridge.inv <- function(SymMat, lambda.prop=.1^4){
    ee <- eigen(SymMat); LL <- ee[["values"]]; TT <- ee[["vectors"]]
    lambda <- sum(LL) * lambda.prop
    inv.mat <- TT %*% diag(1/(pmax(LL, lambda))) %*% t(TT)
    return(inv.mat)
}

## A wrapper to select smoothing parameter automatically
autolambda <- function(Ts, Y, loglams=seq(-2,4,.5), subset=TRUE, nsubset=50){
  K <- length(loglams); lambdas <- 10^loglams
  ngenes <- ncol(Y)
  if (subset && ngenes>nsubset) {
    Ysub <- Y[,sample(ngenes,nsubset)]
  } else {
    Ysub <- Y
  }
  mybasis <- create.bspline.basis(range(Ts), length(Ts)+4-2, 4, Ts)
  gcvs <- rep(0,K); dfs <- rep(0,K)
  for (k in 1:K){
    par.k <- fdPar(mybasis, 2, lambdas[k])
    curves.k <- smooth.basis(Ts, Ysub, par.k)
    gcvs[k] <- mean(curves.k[["gcv"]])
    dfs[k] <- curves.k[["df"]]
  }
  kstar <- which.min(gcvs)
  return(list(gcvs=gcvs, dfs=dfs, lambdas=lambdas,
              kstar=kstar, df.star=dfs[kstar],
              lambda.star=lambdas[kstar]))
}

## useful functions for high-dim ODE project
projected.fnorm <- function(V, W){
  ## projection Frobenius norm.
  if (is.fd(V) && is.fd(W)){             #fPCA objects or curves
    mybasis <- V[["basis"]]
    Beta <- fd(diag(mybasis[["nbasis"]]), mybasis)
    SigmaBeta <- inprod(Beta, Beta); ee <- eigen(SigmaBeta)
    Tbeta <- ee[["vectors"]]
    Lambda.root <- diag(sqrt(ee[["values"]]))
    ## Qv, Qw are the representation of V, W under an orthonormal basis
    Qv <- qr.Q(qr(Lambda.root %*% t(Tbeta) %*% coef(V)))
    Qw <- qr.Q(qr(Lambda.root %*% t(Tbeta) %*% coef(W)))
  } else if (is.matrix(V) && is.matrix(W)){ #pca objects
    Qv <- qr.Q(qr(V)); Qw <- qr.Q(qr(W))
  } else {
    stop("Only the following PCA methods are implemented: fpca, pca, spca.")
  }
  diffmat <- as.matrix(Qv %*% t(Qv) - Qw %*% t(Qw))
  return(norm(diffmat, "F")/sqrt(2))
}

## Given two PCA results, produces the difference and the translation
## matrix between the two.
pdist2 <- function(pca1, pca2){
  parameters <- pca1[["parameters"]]; method=parameters[["method"]]
  if (method=="pca"){
    centers1 <- pca1[["centers"]]; centers2 <- pca2[["centers"]]
    X1 <- pca1[["xhats"]]; X2 <- pca2[["xhats"]]
    Q1 <- qr.Q(qr(X1)); mu1 <- centers1 - Q1 %*% t(Q1) %*% centers1
    Q2 <- qr.Q(qr(X2)); mu2 <- centers2 - Q2 %*% t(Q2) %*% centers2
    mu.dist2 <- sum((mu1-mu2)^2)
  } else if (method=="fpca"){
    lambda <- parameters[["lambda"]]
    mybasis <- pca1[["xhats.curves"]][["basis"]]
    mypar <- fdPar(mybasis, 2, lambda=lambda)  #under-smooth
    Xt1 <- pca1[["xhats.curves"]]; Xt2 <- pca2[["xhats.curves"]]
    meancur1 <- pca1[["meancur"]]; meancur2 <- pca2[["meancur"]]

    Beta <- fd(diag(mybasis[["nbasis"]]), mybasis)
    SigmaBeta <- inprod(Beta, Beta); ee <- eigen(SigmaBeta)
    Tbeta <- ee[["vectors"]]
    Lambda.root <- diag(sqrt(ee[["values"]]))
    Q1 <- qr.Q(qr(Lambda.root %*% t(Tbeta) %*% coef(Xt1)))
    Q2 <- qr.Q(qr(Lambda.root %*% t(Tbeta) %*% coef(Xt2)))

    mu1 <- fd(coef(meancur1) -coef(Xt1)%*%inprod(Xt1, meancur1), mybasis)
    mu2 <- fd(coef(meancur2) -coef(Xt2)%*%inprod(Xt2, meancur2), mybasis)
    mu.dist2 <- as.real(inprod(mu1, mu2))
  } else if  (method=="spca") {
    stop("Currently not available.")
  } else {
    stop("Only the following PCA methods are implemented: fpca, pca, spca.")
  }
  return(norm(Q1%*%t(Q1) - Q2%*%t(Q2),"F")^2 + mu.dist2)
}

## The distance of the resulting ODE of Y
fdist2 <- function(pcode1, pcode2){
  Tlength <- diff(range(pcode1[["Ts"]]))
  A1 <- pcode1[["Ahat"]]; B1 <- pcode1[["Bhat"]]; Binv1 <- pcode1[["Binv"]]
  A2 <- pcode2[["Ahat"]]; B2 <- pcode2[["Bhat"]]; Binv2 <- pcode2[["Binv"]]
  C1 <- solve(A1) %*% (expm(Tlength * A1) - diag(nrow(A1)))
  C2 <- solve(A2) %*% (expm(Tlength * A2) - diag(nrow(A2)))
  qr1 <- qr(B1); Q1 <- qr.Q(qr1); R1 <- qr.R(qr1)
  qr2 <- qr(B2); Q2 <- qr.Q(qr2); R2 <- qr.R(qr2)
  Q12 <- t(Q1) %*% Q2
  diffmat <- as.matrix(R1%*%C1%*% solve(R1) - Q12%*% R2%*% C2 %*% solve(R2) %*% t(Q12))
  return(norm(diffmat, "F")^2)
}

## the procrustes mean on Graff(K,m). Needed by group.pcafun().
graff.mean <- function(pcalist){
  ## given a list of PCA or Graff objects (typically returned by
  ## pcafun()), returns a new PCA mean.
  n <- length(pcalist); parameters <- pcalist[[1]][["parameters"]]
  Ts <- parameters[["Ts"]];  K <- parameters[["K"]]
  method=parameters[["method"]]; lambda <- parameters[["lambda"]]
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
    mean.xhats <- eval.fd(Ts, mean.xhats.curves)
    mean.centers <- as.vector(eval.fd(Ts, mean.meancur))
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
## apply it to another set of data.  Its main purpose is to produce
## the Bhat matrix, but it can be used to produce a full list of other
## objects as well.  Note that the affine projection does not change
## meancur/centers.
reprojection <- function(Y2, PCA1, Bhat.only=FALSE, lambda.prop=.1^4){
  params <- PCA1[["parameters"]]; method <- params[["method"]]
  Ts <- params[["Ts"]]; K <- params[["K"]]
  pcnames <- paste("PC",1:K,sep="")
  if (method=="fpca"){
    Xt <- PCA1[["xhats.curves"]];  bs1 <- Xt[["basis"]]
    mybasis <- Xt[["basis"]]
    meancur <- PCA1[["meancur"]]
    mypar <- fdPar(bs1, 2, lambda=params[["lambda"]])
    ycurves <- smooth.basis(Ts, Y2, mypar)[["fd"]]
    ycurves.centered <- fd(sweep(coef(ycurves),1,as.vector(coef(meancur))),bs1)
    ymat <- coef(ycurves.centered)
    ## Since Xt may only be approx orthonormal, I use the following
    ## safer formula to compute Bhat
    Bhat <- t(ridge.inv(inprod(Xt, Xt), lambda.prop=lambda.prop) %*% inprod(Xt, ycurves.centered))
    ngenes <- ncol(coef(ycurves))
    ## Calculating the total variance
    Beta <- fd(diag(mybasis[["nbasis"]]), mybasis)
    SigmaBeta <- inprod(Beta, Beta); ee <- eigen(SigmaBeta)
    Tbeta <- ee[["vectors"]]; Lambda.root <- diag(sqrt(ee[["values"]]))
    totalvar <- sum((t(ymat) %*% Tbeta %*% Lambda.root)^2)
    ## Now the variance of the projections.
    Lambdas <- sqrt(eigen(inprod(Xt, Xt))[["values"]])
    varlist <- sapply(1:K, function(k) sum((Bhat[,k]*Lambdas[k])^2))
  } else if (method=="pca"){
    X <- PCA1[["xhats"]]
    Y2.centered <- sweep(Y2, 1, PCA1[["centers"]])
    Bhat <- t(ridge.inv(t(X) %*% X, lambda.prop=lambda.prop) %*% t(X) %*% Y2.centered)
    totalvar <- sum(Y2.centered^2)
    varlist <- sapply(1:K, function(k) sum((Bhat[,k] %*% t(X[,k]))^2))
  } else if (method=="spca"){
    ## identical to the method used for pca. In the future may have a
    ## different implementation.
    X <- PCA1[["xhats"]]
    Y2.centered <- sweep(Y2, 1, PCA1[["centers"]])
    Bhat <- t(ridge.inv(t(X) %*% X, lambda.prop=lambda.prop) %*% t(X) %*% Y2.centered)
    totalvar <- sum(Y2.centered^2)
    varlist <- sapply(1:K, function(k) sum((Bhat[,k] %*% t(X[,k]))^2))
  } else {
    stop("Only the following PCA methods are implemented: fpca, pca, spca.")
  }
  rownames(Bhat) <- colnames(Y2); colnames(Bhat) <- pcnames

  if (Bhat.only){
    return(Bhat)
  } else {
    PCA1[["Bhat"]] <- Bhat; PCA1[["varprop"]] <- varlist/totalvar
    return(PCA1)
  }
}

## This function takes B and A, produces a sparse representation of C=BAB-
sparseC <- function(A, B, prior=NULL){
  m <- nrow(B); K <- ncol(B)
  ## dirty trick: find the best combination within 2*K candidate genes
  if (is.null(prior)) {
    np <- min(2*K,m); prior <- 1:(np)
  }
  B0 <- scale(B[prior,], center=FALSE)
  Sigma0 <- B0 %*% t(B0)
  combmat <- combn(np,K)
  dets <- sapply(1:ncol(combmat), function(k) det(Sigma0[combmat[,k],combmat[,k]]))
  bestgenes <- combmat[,which.max(dets)]
  ord <- order(c(bestgenes, setdiff(1:m, bestgenes)))
  ## break the matrix into B1 and B2
  B1 <- B[bestgenes,]; B1inv <- solve(B1); B2 <- B[-bestgenes,]
  Cmat <- cbind(rbind(B1 %*% A %*% B1inv, B2 %*% A %*% B1inv),matrix(0,nrow=m,ncol=(m-K)))
  Cmat <- Cmat[ord,ord]
  attr(Cmat, "bestgenes") <- bestgenes
  return(Cmat)
}

sortequiv <- function(rmat){
    ## This function takes a nx2 matrix of equivalence, returns its
    ## disconnected components.
    nodes <- unique(as.vector(rmat))
    if (nrow(rmat)==0){                 #empty matrix
        return (NULL)
    } else {
        equivs <- lapply(1:nrow(rmat), function(i) rmat[i,])
        if (length(equivs) >1){
            for (i in 1:(length(equivs)-1)){
                for (j in (i+1):length(equivs)){
                    if (any(equivs[[j]] %in% equivs[[i]])) {
                        equivs[[i]] <- unique(c(equivs[[i]], equivs[[j]]))
                        equivs[[j]] <- "null"
                    } 
                }
            }
        } 
        ## now remove all the "null" elements
        equivs <- equivs[equivs != "null"]
    }
    return (equivs)
}


