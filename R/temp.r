 ## ## A garbage can to collect old functions that *might* be useful.

## ## The main function
## PCODE.weighted <- function(y, Ts, K, lambda=0.01, pca.method=c("fpca", "pca", "spca"), lowdim.method=c("pelos", "pda", "two.stage", "FME"), center=FALSE, spca.para=2^seq(K)/2, const=TRUE){
##   pca.method <- match.arg(pca.method); lowdim.method <- match.arg(lowdim.method)
##   ## test if y is a list of subjects or just one subject
##   if (is.list(y)) {                     #many subjects
##     m <- ncol(y[[1]])
##     pca.results <- group.pcafun(y, Ts=Ts, K=K, method=pca.method, center=center, spca.para=spca.para)
##   } else {                              #just one subject
##     m <- ncol(y)
##     pca.results <- pcafun(y, Ts, K=K, lambda=lambda, method=pca.method, center=center, spca.para=spca.para)
##   }
##   ## now use square root of varprop to weight xhats, xhats.curves
##   lambda.root <- sqrt(pca.results[["varprop"]])
##   wxhats <- pca.results[["xhats"]] %*% diag(lambda.root)
##   xhats.curves <- pca.results[["xhats.curves"]]
##   mybasis <- xhats.curves[["basis"]]
##   ## below starts the fitting of the intrinsic system.  Note that we
##   ## use sqrt(varprop) to weight xhats and xhats.curves.
##   wxhats.curves <- fd(coef(xhats.curves) %*% diag(lambda.root), mybasis)
##   intrinsic.system <- lowdim.est(Ts, xhats=wxhats, xhats.curves=wxhats.curves, method=lowdim.method, lambda=lambda, const=const)
##   wAhat <- intrinsic.system[["Ahat"]]; wbvec <- intrinsic.system[["bvec"]]
##   Ahat <- diag(1/lambda.root) %*% wAhat %*% diag(lambda.root)
##   pcnames <- paste("PC",1:K,sep="")
##   dimnames(Ahat) <- list(pcnames, pcnames)
##   bvec <- as.vector(diag(1/lambda.root) %*% wbvec); names(bvec) <- pcnames
##   xhats.fit <- intrinsic.system[["xhats.fit"]] %*% diag(1/lambda.root)
##   colnames(xhats.fit) <- pcnames
##   ## Done fitting intrinsic system. Now translate these weighted
##   ## curves back.
##   wxhats.fit.curves <- intrinsic.system[["xhats.fit.curves"]]
##   xcoefs <- coef(wxhats.fit.curves) %*% diag(1/lambda.root)
##   colnames(xcoefs) <- pcnames
##   xhats.fit.curves <- fd(xcoefs,mybasis)
##   muvec <- matrix(rep(pca.results[["centers"]], m), nrow=length(Ts))
##   y.fit <- xhats.fit %*% t(pca.results[["Bhat"]]) + muvec
##   mucoef <- matrix(rep(coef(pca.results[["meancur"]]), m), ncol=m)
##   y.fit.curves <- fd(coef(xhats.fit.curves) %*% t(pca.results[["Bhat"]]) + mucoef, mybasis)

##   if (is.list(y)) {                     #many subjects
##     rss <- mean(sapply(y, function(yn) sum((yn-y.fit)^2)))/m
##   } else {                              #just one subject
##     rss <- sum((y-y.fit)^2/m)
##   }
##   return (list(Ts=Ts, xhats.fit=xhats.fit,
##                xhats.fit.curves=xhats.fit.curves,
##                y.fit=y.fit, y.fit.curves=y.fit.curves,
##                rss=rss, varprop=pca.results[["varprop"]],
##                Ahat=Ahat, bvec=bvec,
##                Bhat=pca.results[["Bhat"]], Binv=pca.results[["Binv"]],
##                meancur=pca.results[["meancur"]],
##                centers=pca.results[["centers"]],
##                lambda=lambda,
##                pca.results=pca.results,
##                weighted.intrinsic.system=intrinsic.system))
## }


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


## ## This is a full least-square method
## .est.FME <- function(Ts, Xt, xhats, const) {
##   Xcost <- function(pars, Ts, xinit, const){
##     if (const) {
##       Ab <- matrix(pars, ncol=length(xinit)+1)
##       A <- Ab[-1,-1]; b <- Ab[-1,1]
##       out <- ode.fit(Ts, xinit, A, b)
##     } else {
##       A <- matrix(pars, ncol=length(xinit))
##       out <- ode.fit(Ts, xinit, A)
##     }
##     cost <- modCost(model=out, obs=cbind(time=Ts, xhats))
##     return(cost)
##   }

##   ## use 2-stage method to estimate the initial values
##   Ab0 <- .est.2stage(Ts,Xt)

##   X <- eval.fd(Ts, Xt)
##   if (const) {                           #inhomogeneous
##     x0 <- X[1,-1]
##   } else {                              #homogeneous
##     x0 <- X[1,]
##   }
##   ## Make sure x0 and xhats have the same (non.null names)
##   pcnames <- paste("PC",1:length(x0),sep="")
##   names(x0) <- pcnames; colnames(xhats) <- pcnames
  
##   myfit <- modFit(f=Xcost, p=as.vector(Ab0), Ts=Ts, xinit=x0, const=const)
##   Ahat <- matrix(coef(myfit), ncol=ncol(X))
##   return(Ahat)
## }

