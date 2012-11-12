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
  K <- ncol(xhats)
  mybasis <- xhats.curves[["basis"]]
  mypar <- fdPar(mybasis, 2, lambda=lambda)

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
    xhats.fit <- as.matrix(Xfun(pars=coef(myfit), Ts, xinit=z[1,-1], const=const)[,-1])
  } else {                              #homogeneous
    z <- eval.fd(xhats.curves,Ts)
    z.deriv <- eval.fd(deriv(xhats.curves), Ts)
    CC0 <- t(solve(t(z) %*% z) %*% t(z) %*% z.deriv)
    pars0 <- as.vector(CC0)
    myfit <- modFit(f=Xcost, p=pars0, xinit=z[1,])
    Ahat <- matrix(coef(myfit), ncol=K); bvec <- rep(0,K)
    xhats.fit <- as.matrix(Xfun(pars=coef(myfit), Ts, xinit=z[1,], const=const)[,-1])
  }
  pcnames <- paste("PC",1:K,sep="")
  dimnames(Ahat) <- list(pcnames, pcnames)
  names(bvec) <- pcnames
  rownames(xhats.fit) <- Ts
  ## a spline representation of xhats
  xhats.fit.curves <- smooth.basis(Ts, xhats.fit, mypar)[["fd"]]

  return(list(Ahat=Ahat, bvec=bvec, deviance=deviance(myfit), iterations=myfit[["iterations"]],xhats.fit=xhats.fit, xhats.fit.curves=xhats.fit.curves,Times=Ts))
}

## The main function
PCODE <- function(y, Ts, K, lambda=0.01, pca.method=c("fpca", "pca", "spca"), lowdim.method=c("FME","pda"), center=FALSE, spca.para=2^seq(K)/2, const=TRUE){
  pca.method <- match.arg(pca.method); lowdim.method <- match.arg(lowdim.method)
  ## test if y is a list of subjects or just one subject
  if (is.list(y)) {                     #many subjects
    ngenes <- ncol(y[[1]])
    pca.results <- group.pcafun(y, Ts=Ts, K=K, method=pca.method, center=center, spca.para=spca.para)
  } else {                              #just one subject
    ngenes <- ncol(y)
    pca.results <- pcafun(y, Ts, K=K, lambda=lambda, method=pca.method, center=center, spca.para=spca.para)
  }

  intrinsic.system <- lowdim.est(Ts, xhats=pca.results[["xhats"]], xhats.curves=pca.results[["xhats.curves"]], method=lowdim.method, lambda=lambda, const=const)
  xhats.fit <- intrinsic.system[["xhats.fit"]]
  xhats.fit.curves <- intrinsic.system[["xhats.fit.curves"]]
  mybasis <- xhats.fit.curves[["basis"]]

  muvec <- matrix(rep(pca.results[["centers"]], ngenes), nrow=length(Ts))
  y.fit <- xhats.fit %*% t(pca.results[["Bhat"]]) + muvec
  mucoef <- matrix(rep(coef(pca.results[["meancur"]]), ngenes), ncol=ngenes)
  y.fit.curves <- fd(coef(xhats.fit.curves) %*% t(pca.results[["Bhat"]]) + mucoef, mybasis)

  if (is.list(y)) {                     #many subjects
    rss <- mean(sapply(y, function(yn) sum((yn-y.fit)^2)))
  } else {                              #just one subject
    rss <- sum((y-y.fit)^2)
  }
  return (list(Times=Ts, xhats.fit=xhats.fit,
               xhats.fit.curves=xhats.fit.curves,
               y.fit=y.fit, y.fit.curves=y.fit.curves,
               rss=rss,
               Ahat=intrinsic.system[["Ahat"]],
               bvec=intrinsic.system[["bvec"]],
               Bhat=pca.results[["Bhat"]], Binv=pca.results[["Binv"]],
               pca.results=pca.results,
               intrinsic.system=intrinsic.system))
}

## between-subject fitting.  Note that y0.new must be a column vector.  As of ver 0.02, this prediction function does not work well with const != 0 case.
predict.pcode1 <- function(pcode.fit, y0.new, Ts.new="same"){
  ## This function only works with un-centered version as of 09/08/2012.
  Ts <- pcode.fit[["Times"]];
  if (is.character(Ts.new) && Ts.new=="same"){
    Ts.new <- Ts
  }
  xhat0 <- as.vector(pcode.fit[["Binv"]] %*% y0.new)
  xhats.fit <- Xfun(pars=as.vector(pcode.fit[["Ahat"]]), Ts.new, xinit=xhat0, const=FALSE)[,-1]
  y.fit <- as.matrix(xhats.fit) %*% t(pcode.fit[["Bhat"]])
  return (y.fit)
}

## wrapper for between subject cross-validation. As of ver 0.02, this function does not work with const=TRUE case.
CV.group <- function(Ylist, Ts, K, const=FALSE, ...){
  N <- length(Ylist)                    #number of subjects
  y.fit.list <- list(); rss <- rep(0,N)
  for (n in 1:N){
    Ylist.n <- Ylist[-n]; Yn <- Ylist[[n]]
    meansys <- PCODE(Ylist.n, Ts=Ts, K=K, ...)
    y.fit.list[[n]] <- predict.pcode1(meansys, Yn[1,])
    rss[n] <- sum((Yn-y.fit.list[[n]])^2)
  }
  return(list(y.fit.list=y.fit.list, rss=rss))
}


## ## wrapper for plotting
## plot.pcode <- function(pcode.result){

## }
