## 6/24/2013.  TODO: per-trajectory CV

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

## wrapper for low-dim (intrinsic) ODE estimation. Input: Ts; xhats
## and xhats.init estimated from pcafun(); C.init estimated from
## C.init.est(). const=T/F: whether the equation system is homogeneous
## or inhomogeneous.
lowdim.est <- function(Ts, xhats, xhats.curves, weights=NULL, est.method=c("pda", "two.stage", "pda0", "two.stage0"), est.pen=.1^4, stab.method=c("eigen-bound", "eigen-bound2", "random", "zero", "none"), refine.method=c("pelos", "none"), lambda=.1^4, const=FALSE, verbose=FALSE, ...){
    K <- ncol(xhats); pcnames <- paste("PC",1:K,sep="")
    ## must make sure both xhats and xinit has the same, non-null names
    colnames(xhats) <- pcnames
    mybasis <- xhats.curves[["basis"]]
    mypar <- fdPar(mybasis, 2, lambda=lambda)
    ## normalize weights
    if (!is.null(weights)) weights <- weights/sum(weights)

    ## now deal with the inhomogeneous case
    if (const) {
        Xt <- fd(cbind("Const"=1, coef(xhats.curves)), mybasis)
        xhats <- cbind(1, xhats)
        ## pad weights with 1/(K+1) for now
        if (!is.null(weights)) weights <- c(1/(K+1), K/(K+1)*weights)
    } else {
        Xt <- xhats.curves
    }
    ## Different backends
    est.method <- match.arg(est.method)
    stab.method <- match.arg(stab.method)
    refine.method <- match.arg(refine.method)

    Ab1 <- .est(Ts, Xt, method=est.method, est.pen=est.pen)
    Ab2 <- .stabilize(Ab1, Ts, method=stab.method)
    RF <- .refine(Ts, Xt, xhats, Ab2, method=refine.method, weights=weights, verbose=verbose)
    Ab <- RF[["Ahat"]]; x0 <- RF[["x0"]]

    ## Disentangle A and b.
    if (const) {
        Ahat <- Ab[-1, -1]; bvec <- as.vector(Ab[-1, 1]); x0 <- x0[-1]
    } else {
        Ahat <- Ab; bvec <- rep(0,K)
    }
    dimnames(Ahat) <- list(pcnames, pcnames)
    names(bvec) <- pcnames

    ## Now produce some additional information based on A and b.
    ## don't need the "time" column anymore
    xhats.fit <- as.matrix(ode.fit(Ts, x0, Ahat, bvec)[,-1])
    rownames(xhats.fit) <- Ts
    ## a spline representation of xhats
    xhats.fit.curves <- smooth.basis(Ts, xhats.fit, mypar)[["fd"]]

    return(list(Ahat=Ahat, bvec=bvec,
                xhats.fit=xhats.fit,
                xhats.fit.curves=xhats.fit.curves,Ts=Ts))
}

## The main function
PCODE <- function(y, Ts, K, lambda=.1^5, pca.method=c("fpca", "pca"), weight=c("varprop", "none"), prescreen.prop=.7, lambda1=1e-5, lambda2="auto", thresh=1e-11, cutoff=.05, est.method=c("pda", "two.stage", "pda0", "two.stage0"), est.pen=.1^4, stab.method=c("eigen-bound", "eigen-bound2", "random", "zero", "none"), refine.method=c("none", "pelos"), compensate=TRUE, gamma.approx=TRUE, center=FALSE, const=FALSE, verbose=FALSE, ...){
    ## The ... arguments are used by .backfit().
    pca.method <- match.arg(pca.method)
    weight <- match.arg(weight)
    est.method <- match.arg(est.method)
    stab.method <- match.arg(stab.method)
    refine.method <- match.arg(refine.method)
    ## test if y is a list of subjects or just one subject
    if (is.list(y)) {                     #many subjects
        m <- ncol(y[[1]])
        pca.results <- group.pcafun(y, Ts=Ts, K=K, lambda=lambda, method=pca.method, center=center)
    } else {                              #just one subject
        m <- ncol(y)
        pca.results <- pcafun(y, Ts, K=K, lambda=lambda, method=pca.method, center=center)
    }

    xhats=pca.results[["xhats"]]
    Bhat <- pca.results[["Bhat"]]; xhats.curves <- pca.results[["xhats.curves"]]
    if (weight=="varprop") {
        weights <- pca.results[["varprop"]]
    } else if (weight=="none"){
        weights <- NULL
    } else {
        stop("Supported weighting methods are varprop or none.")
    }

    intrinsic.system <- lowdim.est(Ts, xhats=xhats, xhats.curves=xhats.curves, weights=weights,
                                   est.method=est.method, est.pen=est.pen, stab.method=stab.method,
                                   refine.method=refine.method, lambda=lambda, const=const, verbose=verbose)
    xhats.fit <- intrinsic.system[["xhats.fit"]]
    xhats.fit.curves <- intrinsic.system[["xhats.fit.curves"]]
    Ahat <- intrinsic.system[["Ahat"]]

    pcnames <- paste("PC",1:K,sep=""); colnames(xhats.fit) <- pcnames
    muvec <- matrix(rep(pca.results[["centers"]], m), nrow=length(Ts))
    y.fit <- xhats.fit %*% t(Bhat) + muvec
    mucoef <- matrix(rep(coef(pca.results[["meancur"]]), m), ncol=m)
    mybasis <- xhats.curves[["basis"]]
    y.fit.curves <- fd(coef(xhats.fit.curves) %*% t(Bhat) + mucoef, mybasis)
    ## backfit the original eqn system
    bfit <- backfit(Bhat, Ahat, xhats.fit[1,], xhats, Ts, lambda1=lambda1, lambda2=lambda2, compensate=compensate, gamma.approx=gamma.approx, thresh=thresh, cutoff=cutoff, prescreen.prop=prescreen.prop, ...)
    Theta <- bfit[["Theta"]]
    ## Compute RSS
    if (is.list(y)) {                     #many subjects
        rss <- mean(sapply(y, function(yn) sum((yn-y.fit)^2)))/m
    } else {                              #just one subject
        rss <- sum((y-y.fit)^2/m)
    }
    return (list(Ts=Ts, xhats.fit=xhats.fit,
                 xhats.fit.curves=xhats.fit.curves,
                 y.fit=y.fit, y.fit.curves=y.fit.curves,
                 rss=rss, varprop=pca.results[["varprop"]],
                 Ahat=Ahat,
                 bvec=intrinsic.system[["bvec"]],
                 Bhat=Bhat, Binv=pca.results[["Binv"]],
                 bfit=bfit,
                 Theta=Theta,
                 meancur=pca.results[["meancur"]],
                 centers=pca.results[["centers"]],
                 lambda=lambda,
                 pca.results=pca.results,
                 intrinsic.system=intrinsic.system))
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
CV <- function(Y, Ts, Ks, folds=10, center=FALSE, const=FALSE, refine.method="none", ...){
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
            trainsys.j <- PCODE(Y.j, Ts=Ts.j, K=k, center=center, const=const, refine.method=refine.method, ...)
            y.fits <- predict.pcode1(trainsys.j, Y[1,], Ts.new=Ts[idx[[j]]])
            sum((Y[idx[[j]], ] - y.fits)^2)
        }
    }
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

## Estimating K.
## Kest <- function(pcaobj, Ts, Ks, pca.method=c("fpca", "pca", "spca"), method=c("varprop", "cv", "pca.cv"), var.prop=.1^3, ...){
##     pca.results <- pcafun(y, Ts, K=K, lambda=lambda, method=pca.method, center=center, spca.para=spca.para)
##     if (method=="varprop"){
##         ## just cutoff at a given varprop
        
##     } else if (method=="cv"){
##     } else if (method=="pca.cv"){
##     } else {
##         stop(paste("Method",method,"is not implemented in function Kest()!"))
##     }
    

##     return(Kstar)
## }


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
        if (!is.null(true.y.curves)) lines(true.y2.curves[i], col="blue")
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
            if (!is.null(true.y.curves)) lines(true.y2.curves[i], col="blue")
        }
    }
}
