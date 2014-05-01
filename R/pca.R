## A wrapper for conducting various PCA and produces the same
## output. As a "bonus", it produces the initial values of X through
## smoothing.
pcafun <- function(y,Ts,K, lambda=0.01, method=c("fpca", "pca"),
                   center=FALSE){
    method <- match.arg(method)
    mybasis <- create.bspline.basis(range(Ts), length(Ts)+4-2, 4, Ts)
    mypar <- fdPar(mybasis, 2, lambda=lambda)  #under-smooth
    ycurves <- smooth.basis(Ts, y, mypar)[["fd"]]
    pcnames <- paste("PC",1:K,sep="")
    if (method=="fpca"){
        rr <- pca.fd(ycurves, centerfns=center, nharm=K)
        varprop <- rr[["varprop"]]; names(varprop) <- pcnames
        ## as for the centers, fPCA is an odd ball.
        if (center) {
            meancur <- rr[["meanfd"]]
            ## centers <- as.vector(eval.fd(rr[["meanfd"]],Ts))
            centers <- as.vector(eval.fd(Ts,rr[["meanfd"]]))
        } else {
            centers <- rep(0.0,nrow(y))
            meancur <- fd(rep(0.0, mybasis$nbasis), mybasis)
        }; names(centers) <- rownames(y)
        ## 3/20/2014. Make Bhat an (near-)orthogonal matrix
        Bhat0 <- rr[["scores"]]
        scale.coefs <- sqrt(nrow(Bhat0) * rr[["values"]][1:K])
        Bhat <- Bhat0 %*% diag(1/scale.coefs)
        dimnames(Bhat) <- list(colnames(y),pcnames)
        Bhat.qr <- qr(Bhat)
        ## rotation <- qr.Q(Bhat.qr); dimnames(rotation) <- list(colnames(y),pcnames)
        ## sdev <- abs(diag(qr.R(Bhat.qr))); names(sdev) <- pcnames
        ## xhats.curves are just the harmonics (scaled by the eigenvals)
        xhats.curves0 <- rr[["harmonics"]]
        xhats.curves <- fd(coef(xhats.curves0) %*% diag(scale.coefs), mybasis)
        ## for fPCA, we evaluate eigen-functions at Ts directly. They are
        ## near, but not exactly, orthogonal.
        ## xhats <- eval.fd(xhats.curves, Ts)
        xhats <- eval.fd(Ts, xhats.curves)
    } else if (method=="pca"){
        rr <- prcomp(t(y), center=center)
        varlist <- rr[["sdev"]]^2; varprop <- varlist[1:K]/sum(varlist)
        names(varprop) <- pcnames
        if (center) {
            centers <- rr[["center"]]
        } else {
            centers <- rep(0.0,nrow(y))
        }
        meancur <- smooth.basis(Ts, centers, mypar)[["fd"]]
        ## 3/20/2014. Scale the Bhat matrix
        Bhat0 <- rr[["x"]][,1:K]
        scale.coefs <- rr[["sdev"]][1:K] * sqrt(nrow(Bhat0)-1)
        Bhat <- Bhat0 %*% diag(1/scale.coefs)
        xhats0 <- rr[["rotation"]][,1:K]  #The first K eigenvectors
        xhats <- xhats0 %*% diag(scale.coefs)
        ## The xhat curves are represented by smoothed splines of xhats
        xhats.curves <- smooth.basis(Ts, xhats, mypar)[["fd"]]
    } else {
        stop("Only the following PCA methods are implemented: fpca, pca.")
    }
    ## finally, estimate Binv, the generalized inverse of Bhat. This
    ## matrix is used in between-subject fitting
    ## 09/08/2012. The following shortcut is no good: it has random signs.
    ## Binv <- diag(1/sdev) %*% t(rotation)
    Binv <- solve(t(Bhat) %*% Bhat) %*% t(Bhat)

    return(list(varprop=varprop, centers=centers, meancur=meancur, Bhat=Bhat, Binv=Binv, xhats=xhats, xhats.curves=xhats.curves, parameters=list(Ts=Ts, K=K, lambda=lambda, method=method,center=center)))
}

## group.pcafun estimates a centroid for a group of data (xhats).
group.pcafun <- function(Ylist, Ts, K, method=c("fpca", "pca"), ...){
    method <- match.arg(method); N <- length(Ylist)
    pcalist <- lapply(Ylist, pcafun, Ts=Ts, K=K, method=method, ...)
    ## The Graff mean of these PCA results
    pca.mean <- graff.mean(pcalist)
    ## recompute Bhat for each data
    reproj <- lapply(Ylist, reprojection, PCA1=pca.mean)
    Bhats2 <- lapply(reproj, function(rr) rr[["Bhat"]])
    Bmean <- Reduce("+", Bhats2)/N
    Bmean.inv <- solve(t(Bmean) %*% Bmean) %*% t(Bmean)
    ## calculate the varprop
    varprops <- lapply(reproj, function(rr) rr[["varprop"]])
    varprop <- Reduce("+", varprops)/N
    ## append these two terms in pca.mean. Note that affine projection
    ## does not alter meancur/centers

    pca.mean[["Bhat"]] <- Bmean; pca.mean[["Binv"]] <- Bmean.inv
    pca.mean[["varprop"]] <- varprop
    return(pca.mean)
}

## Cross-validation functions
predict.pca <- function(pca.fit, Ts.new=NULL){
    Bhat <- pca.fit[["Bhat"]]
    if (is.null(Ts.new)){                 #just use the original Ts
        xhats.fit <- pca.fit[["xhats"]]; centers <- pca.fit[["centers"]]
    } else {                              #estimate the centers for new Ts
        centers <- eval.fd(Ts.new, pca.fit[["meancur"]])
        xhats.fit <- eval.fd(Ts.new, pca.fit[["xhats.curves"]])
    }
    y.fit <- sweep(as.matrix(xhats.fit) %*% t(Bhat), 1, -centers)
    return (y.fit)
}


PCACV <- function(Y, Ts, Ks, folds=5, lambda=.1^5, pca.method="fpca", center=FALSE, const=FALSE,  ...){
    J <- length(Ts)
    ## Generate the index of testing time points first.  Due to the
    ## fact that it is hard to extrapolate, we will drop the first and
    ## last time points for now.
    idx <- lapply(2:min(folds+1,J-1), function(j) seq(j,J-1,folds))
    rss <- rep(0, length(Ks)); names(rss) <- paste("K", Ks, sep="")
    for (k in Ks){
        rss[paste("K",k,sep="")] <- foreach(j=1:length(idx), .combine="+") %dopar%{
            ## Ts=Ts; K=K
            Y.j <- Y[-idx[[j]],]; Ts.j <- Ts[-idx[[j]]]
            trainsys.j <- pcafun(Y.j, Ts=Ts.j, K=k, lambda=lambda, method=pca.method, center=center)
            y.fits <- predict.pca(trainsys.j, Ts.new=Ts[idx[[j]]])
            sum((Y[idx[[j]], ] - y.fits)^2)
        }
    }
    return(rss)
}


