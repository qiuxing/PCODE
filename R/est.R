## The two-stage backends

.est <- function(Ts, Xt, method="pda"){
    if (method=="two.stage"){
        X <- eval.fd(Ts, Xt); X.deriv <- eval.fd(Ts,deriv(Xt))
        Ahat <- t(X.deriv) %*% X %*% solve(t(X) %*% X)
    } else if (method=="pda"){
        Sigma.xx <- inprod(Xt, Xt)
        Sigma.xderivx <- inprod(deriv(Xt), Xt)
        Ahat <- Sigma.xderivx %*% solve(Sigma.xx)
    } else {
        stop(paste("Method",method,"is not implemented in function .est()!"))
    }
    return(Ahat)
}

## The stabilizing procedure
.stabilize <- function(Ahat, Ts, method="eigen-bound", tol=1e-5){
    if (method=="eigen-bound" | method=="eigen-bound2"){
        ## First, compute the eigen-bounds from the Ts
        Tlength=max(Ts)-min(Ts); dT <- min(abs(diff(Ts)))
        real.upper <- -1e-4 / Tlength   #slightly negative for the real eigenvalues
        real.lower <- -3/Tlength        #e^-3 is approx. 0.05
        imaginary.bound=pi/dT           #at least 1/2 period within one dT

        Tmat <- eigen(Ahat)[["vectors"]]
        lambdas <- eigen(Ahat)[["values"]]
        N <- length(lambdas)
        ## deal with near-tied eigenvalues. Add diag to avoid 0
        distmat <- Mod(outer(lambdas, lambdas, "-"))+tol*diag(N)
        if (min(distmat)<tol){
            ties <- which(distmat<tol & upper.tri(distmat), arr.ind = TRUE)
            equivsets <- sortequiv(ties)
            ## re-orthogonalize eigen vectors
            for (i in length(equivsets)){
                lambdas[equivsets[[i]]] <- mean(lambdas[equivsets[[i]]])
                Tmat[, equivsets[[i]]] <- qr.Q(qr(Tmat[, equivsets[[i]]]))
            }
        }
        lambdas2 <- complex(real=pmax(pmin(Re(lambdas), real.upper), real.lower),
                            imaginary=pmax(pmin(Im(lambdas), imaginary.bound), -1*imaginary.bound))
        Sigma <- diag(lambdas2)
        ## Note that I can't use more efficient Conj(t(Tmat)) here
        ## because when a) Ahat is not of full rank; b) there are
        ## cerntain nontrivial invariant subspaces, Tmat is NOT
        ## unitary.
        Ahat2a <- Re(Tmat %*% Sigma %*% solve(Tmat))
    }
    
    if (method=="eigen-bound"){
        Ahat2 <- Ahat2a
    } else if (method=="eigen-bound2") {
        Ahat2 <- 0.1*Ahat2a             #This is a compromise between eigen-bound and zero
    } else if (method=="random"){
        sigma=.1/Tlength
        N <- nrow(Ahat)
        Ahat2 <- sigma*matrix(rnorm(N^2), nrow=N, ncol=N)
    } else if (method=="zero"){
        N <- nrow(Ahat)
        Ahat2 <- matrix(0, nrow=N, ncol=N)
    } else if (method=="divide"){
        lambdas <- eigen(Ahat)[["values"]]
        Ahat2 <- bound/max(Re(lambdas)) * Ahat
    } else if (method=="none") {
        Ahat2 <- Ahat
    } else {
        stop(paste("Method",method,"is not implemented in function .stabilize()!"))
    }
    return(Ahat2)
}

## refinement methods. A0 is typically the result of .stabilize(.est(...))
.refine <- function(Ts, Xt, xhats, A0, method="pelos", weights=NULL, verbose=FALSE, ...) {
    if (method=="pelos"){
        ## Leqin's method.  
        J <- nrow(xhats); K <- nrow(A0)
        ## A0 <- cbind(x=as.vector(A0), 1:length(as.vector(A0)))
        ## initial values
        x0 <- as.vector(eval.fd(Ts[1], Xt))
        ## transform xhats into the sparse format
        if (is.null(weights)) {
            weights <- rep(1/K,K*J)
        } else {
            weights <- rep(weights/sum(weights), each=J)
        }
        ## OBS <- cbind(value=as.vector(xhats), var=rep(1:K, each=J),Time=rep(Ts,K), weight=weights)
        OBS <- rbind(Ts, t(xhats))
        if (verbose) {
            myfit <- pelos(OBS, structure=matrix(TRUE, K, K),
                           coefficient=A0, initial_value=x0,
                           num_print=8, ...)
        } else {
            capture.output(myfit <- pelos(OBS, structure=matrix(TRUE, K, K), coefficient=A0, initial_value=x0, num_print=0, ...))
        }
        ## Ahat <- matrix(myfit[["coefficient"]], nrow=K)
        Ahat <- myfit[["coefficient"]]
        x0 <- myfit[["initial_value"]]
    } else if (method=="none"){
        Ahat <- A0
        x0 <- t(eval.fd(Ts[1], Xt))
    } else {
        stop(paste("Method",method,"is not implemented in function .refine()!"))
    }
    return(list("Ahat"=Ahat, "x0"=x0))
}

.backfit <- function(Bhat, Astar, X0, Ts, method=c("linearization", "lasso", "gen.inv"), prior=NULL, ...){
    method <- match.arg(method)
    if (method=="lasso"){
        m <- nrow(Bhat)
        Y <- t(Astar) %*% t(Bhat)
        Theta <- matrix(0, nrow=m, ncol=m)
        for (i in 1:m){
            larscoefs <- coef(lars(x=t(Bhat), y=Y[,i], intercept=FALSE, use.Gram=FALSE, ...))
            Theta[i,] <- larscoefs[nrow(larscoefs),]
        }
    } else if (method=="linearization"){
        Theta <- Linearization(Astar, Bhat, X0, Ts, ...)
    } else if (method=="gen.inv"){
        Binv <- solve(t(Bhat) %*% Bhat) %*% t(Bhat)
        Theta <- Bhat %*% Astar %*% Binv
    } else {
        stop("Supported weighting methods are lasso and gen.inv.")
    }
    return(Theta)
}


