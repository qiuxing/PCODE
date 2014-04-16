.SigmaSeries <- function(A, X0, Ts, coefs=c("integral", "sum"), level=5){
    coefs <- match.arg(coefs)
    p <- nrow(A)
    ## compute the building blocks first
    An <- lapply(0:level, function(n) A %^% n)
    AnX <- lapply(0:level, function(n) as.vector(An[[n+1]] %*% X0))
    AA <- array(0, c(1+level,1+level,p,p))
    AXXA <- array(0, c(1+level,1+level,p,p))
    for (n1 in 0:level){
        for (n2 in 0:level){
            AA[n1+1, n2+1, ,] <- t(An[[n1+1]]) %*% An[[n2+1]]
            AXXA[n1+1, n2+1, ,] <- AnX[[n1+1]] %*% t(AnX[[n2+1]])
        }
    }
    ## put everything together
    Sigma <- matrix(0, p^2, p^2)
    for (m1 in 0:level){
        for (n1 in 0:level){
            for (m2 in 0:level){
                for (n2 in 0:level){
                    if (coefs == "integral"){
                        cc <- 1/(m1+n1+m2+n2+3)
                    } else if (coefs == "sum"){
                        cc <- mean(Ts^(m1+n1+m2+n2+2))
                    } else {
                        stop("coefs must be: 1. integral; 2. sum.")
                    }
                    sm <- AXXA[n1+1, n2+1, ,] %x% AA[m1+1, m2+1, ,]
                    Sigma <- Sigma + cc*sm/(factorial(m1+n1+1) *factorial(m2+n2+1))
                }}}}
    return(Sigma)
}

.vf <- function(A, X0, Ts, coefs=c("integral", "sum"), level=5){
    coefs <- match.arg(coefs)
    p <- nrow(A)
    combs <- expand.grid(0:level, 0:level)
    V <- matrix(0, p, p)
    XX <- X0 %*% t(X0)
    for (i2 in 1:nrow(combs)) {
        cc <- unlist(combs[i2,])
        n1 <- cc[1]; n2 <- cc[2]
        if (coefs == "integral"){
            cc <- 1/(n1+n2+3)
        } else if (coefs == "sum"){
            cc <- mean(Ts^(n1+n2+2))
        } else {
            stop("coefs must be: 1. integral; 2. sum.")
        }
        ss <- (A %^% n1) %*% XX %*% t(A %^% n2)
        V <- V + cc*ss/(factorial(n1+1) * factorial(n2+1))
    }
    return(V)
}

gammafun <- function(A, B, X0, Y, Ts){
    ## First, compute the residual vector
    p <- length(X0); J <- length(Ts)
    ## we will drop Ts[0] in calculation of integral
    Xhat <- t(sapply(Ts[-1], function(tt) expm(tt*A) %*% X0))
    tLX <- array(0, c(p*(J-1), p, p))
    for (k in 1:p){
        for (l in 1:p){
            Ekl <- matrix(0, p, p); Ekl[k,l] <- 1
            L.kl <- sapply(Ts[-1], function(tt) tt*expmFrechet(tt*A, Ekl, expm=F)[["Lexpm"]] %*% X0)
            tLX[, k, l] <-  as.vector(t(L.kl))
        }}
    Zt <- as.vector(t(Y[-1,] - Xhat)) + sapply(1:(p*(J-1)), function(tt) sum(tLX[tt, ,] * A))
    dZt <- Zt*diff(Ts)
    gamma <- as.vector(B %*% apply(tLX, c(2,3), function(vec) sum(vec * dZt)) %*% t(B))
    return(gamma)
}

linearize <- function(A, B, X0, Y, Ts, compensate=TRUE, gamma.approx=TRUE, coefs=c("integral", "sum"), level=5, eigen.prop.thresh=0.1^5, lambda2.min=0.1^5){
    ## This function transforms the original ODE fitting problem into
    ## a LASSO regression problem by computational shortcuts
    coefs <- match.arg(coefs)           #for Sigma/V
    p <- length(X0); P <- nrow(B)
    ## First, let's linearize the small system
    Sigma <- .SigmaSeries(A, X0, Ts, coefs=coefs, level=level)
    ## gamma <- Sigma %*% as.vector(A)
    V <- .vf(A, X0, Ts, coefs=coefs, level=level)
    ## The simplified XX and YY
    ee <- eigen(Sigma); eeV <- eigen(V)
    LL0 <- ee[["values"]]; LLV0 <- eeV[["values"]]
    ## determine how many eigenvalues we need
    eigen.thresh <- LL0[1]*eigen.prop.thresh
    LL <- LL0[LL0>eigen.thresh]
    pk <- length(LL); TA <- ee[["vectors"]][, 1:pk]
    LLV <- LLV0[LLV0>eigen.thresh]
    TV <- eeV[["vectors"]][, 1:length(LLV)]
    ## The big system
    BB <- B %x% B
    if (P-p >0){
        Bperp <- eigen(diag(P) - B %*% t(B))[["vectors"]][, 1:(P-p)]
        TW <- (B %*% TV) %x% Bperp
        TB <- cbind(BB %*% TA, TW)
        LLB <- c(LL, rep(LLV, each=P-p))
        if (pk == p^2 & length(LLV) == p){  #no cutoff
            lambda2 <- lambda2.min
        } else {                            #some values are cut
            lambda2 <- max(lambda2.min, setdiff(c(LL0, LLV0), c(LL, LLV)))
        }
    } else {
        TB <- BB %*% TA
        LLB <- LL
        if (pk == p^2){  #no cutoff
            lambda2 <- lambda2.min
        } else {                            #some values are cut
            lambda2 <- max(lambda2.min, setdiff(LL0, LL))
        }
    }
    if (compensate) {
        ## The compensated LL
        LLB2 <- LLB - lambda2
    } else {
        LLB2 <- LLB
    }
    XX <- diag(sqrt(LLB2)) %*% t(TB)
    ## Compute YY based on shortcut or definition
    if (gamma.approx) {
        YY <- drop((diag(sqrt(LLB2)) + lambda2 * diag(1/sqrt(LLB2))) %*% t(TB) %*% as.vector(B %*% A %*% t(B)))
    } else {
        gamma <- gammafun(A, B, X0, Y, Ts)
        YY <- drop(diag(1/sqrt(LLB2)) %*% t(TB) %*% (gamma - lambda2*as.vector(B %*% A %*% t(B))))
    }
    return(list(XX=XX, YY=YY, V=V, lambda2=lambda2))
}

.prescreen <- function(A, X, method=c("L1", "L2", "marginal", "none"), prop=.5){
    ## A wrapper for various pre-screening methods
    ## First, normalize X so that their columns are comparable
    X.std <- X %*% diag(as.vector(A))
    method <- match.arg(method)
    if (method=="none"){
        return (list(X=X, zeroindx=integer(0)))
    } else if (method=="L1"){
        Svec <- colSums(abs(X.std))
    } else if (method=="L2"){
        Svec <- colSums(X.std^2)
    } else if (method=="marginal"){
        Svec <- colSums(X.std)^2
    } else {
        stop("Only these pre-screening methods are implemented: L2 (default), L1, and none.")
    }
    zeroindx <- which(Svec<=quantile(Svec, prop))
    if (length(zeroindx)>0){
        X2 <- X[, -zeroindx]
    } else {
        X2 <- X
    }
    return (list(X=X2, zeroindx=zeroindx))
}

lin2cv.enet <- function(linobj, A, lambda2="auto", K=5, s.seq=seq(.5, 1, .01), prescreen="L1", prescreen.prop=.5, ...){
    X <- linobj[["XX"]]; Y <- linobj[["YY"]]
    p <- sqrt(ncol(X))
    ## Pre-screening
    PS <- .prescreen(A, X, method=prescreen, prop=prescreen.prop)
    X2 <- PS[["X"]]; zeroindx <- PS[["zeroindx"]]
    if (lambda2=="auto") lambda2 <- linobj[["lambda2"]]
    ENCV <- cv.enet(X2, Y, K=K, lambda=lambda2, s=s.seq, mode="fraction", intercept=FALSE, normalize=FALSE, eps=.1^6, ...)
    best.s <- ENCV[["s"]][which.min(ENCV[["cv"]])]
    return(list(encv.obj=ENCV, best.s=best.s, zeroindx=zeroindx))
}

lin2enet <- function(linobj, A, lambda2="auto", prescreen="L1", prescreen.prop=.5, ...){
    X <- linobj[["XX"]]; Y <- linobj[["YY"]]
    p <- sqrt(ncol(X))
    ## Pre-screening
    PS <- .prescreen(A, X, method=prescreen, prop=prescreen.prop)
    X2 <- PS[["X"]]; zeroindx <- PS[["zeroindx"]]
    if (lambda2=="auto") lambda2 <- linobj[["lambda2"]]
    ## limit number of LARS-EN steps to 3*length(Y) instead of 50*length(Y).
    ## N.steps <- 3*length(Y)
    EN <- enet(X2, Y, lambda=lambda2, intercept=FALSE, normalize=FALSE, eps=.1^5, ...)
    return(list(enet.obj=EN, zeroindx=zeroindx))
}

lin2glmnet <- function(linobj, A, lambda2="auto", prescreen="L1", prescreen.prop=.5, ...){
    X <- linobj[["XX"]]; Y <- linobj[["YY"]]
    p <- sqrt(ncol(X))
    ## Pre-screening
    PS <- .prescreen(A, X, method=prescreen, prop=prescreen.prop)
    X2 <- PS[["X"]]; zeroindx <- PS[["zeroindx"]]
    if (lambda2=="auto") lambda2 <- linobj[["lambda2"]]
    ## Let us use an arbitrary alpha for now.  It seems that we can't
    ## control lambda2 directly in glmnet
    GN <- glmnet(X2, Y, alpha=.99, intercept=FALSE, standardize=FALSE, thresh=.1^5, ...)
    return(list(glmnet.obj=EN, zeroindx=zeroindx))
}


glmnet2A <- function(lin.glmnet.obj, cutoff=.1, s=200, ...){
    GN <- lin.glmnet.obj[["glmnet.obj"]]; zeroindx=lin.enet.obj[["zeroindx"]]
    Avec2 <- coef(GN, s=s)
    p <- sqrt(length(Avec2) + length(zeroindx))
    ## Add pre-determined zeros back. Needs to be written in more efficient code
    if (length(zeroindx)>0){
        Avec <- rep(0, p^2)
        Avec[-zeroindx] <- Avec2
    } else {
        Avec <- Avec2
    }
    ## cut really small values
    Avec.cut <- replace(Avec, abs(Avec)<cutoff, 0)
    Ahat <- matrix(Avec.cut, p, p)
    return(Ahat)
}

enet2A <- function(lin.enet.obj, cutoff=.1, s=.9, ...){
    EN <- lin.enet.obj[["enet.obj"]]; zeroindx=lin.enet.obj[["zeroindx"]]
    Avec2 <- predict(EN, s=s, type="coef", mode="fraction")[["coefficients"]]
    p <- sqrt(length(Avec2) + length(zeroindx))
    ## Add pre-determined zeros back. Needs to be written in more efficient code
    if (length(zeroindx)>0){
        Avec <- rep(0, p^2)
        Avec[-zeroindx] <- Avec2
    } else {
        Avec <- Avec2
    }
    ## cut really small values
    Avec.cut <- replace(Avec, abs(Avec)<cutoff, 0)
    Ahat <- matrix(Avec.cut, p, p)
    return(Ahat)
}

## The big wrapper function
Linearization <- function(Astar, Bhat, X0, Y, Ts, compensate=TRUE, gamma.approx=TRUE, coefs=c("integral", "sum"), level=5, eigen.prop.thresh=0.1^5, lambda2.min=0.1^5, prescreen="L1", prescreen.prop=.5, cutoff=0.1, L1pen=0.9, ...){
    coefs <- match.arg(coefs)
    Binv <- solve(t(Bhat) %*% Bhat) %*% t(Bhat)
    ThetaHat <- Bhat %*% Astar %*% Binv
    ## 3-step computation
    LIN <- linearize(Astar, Bhat, X0, Y, Ts, compensate=compensate, gamma.approx=gamma.approx, coefs=coefs, level=level, eigen.prop.thresh=eigen.prop.thresh, lambda2.min=lambda2.min)
    LinEnet <- lin2enet(LIN, ThetaHat, lambda2=LIN[["lambda2"]], prescreen=prescreen, prescreen.prop=prescreen.prop, ...)
    Ahat2 <- enet2A(LinEnet, cutoff=cutoff, s=L1pen)
    return(Ahat2)
}
