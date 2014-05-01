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
    ## min of max Svec/1e-3 or quantile, to make it more robust
    Scut <- max(quantile(Svec, prop), max(Svec)*1e-3)
    zeroindx <- which(Svec<=Scut)
    if (length(zeroindx)>0){
        X2 <- X[, -zeroindx]
    } else {
        X2 <- X
    }
    return (list(X=X2, zeroindx=zeroindx))
}

lin2glmnet <- function(linobj, Theta0, lambda1=0.5, lambda2="auto", thresh=1e-11, prescreen="L1", prescreen.prop=.7, cutoff=.05, ...){
    X <- linobj[["XX"]]; Y <- linobj[["YY"]]
    p <- sqrt(ncol(X))
    ## Pre-screening
    PS <- .prescreen(Theta0, X, method=prescreen, prop=prescreen.prop)
    X2 <- PS[["X"]]; zeroindx <- PS[["zeroindx"]]
    if (lambda2=="auto") lambda2 <- linobj[["lambda2"]]
    ## glmnet needs a little computation
    N <- length(Y); sdY <- sd(Y) * sqrt((N-1)/N)
    X.std <- X2/sdY; Y.std <- Y/sdY
    lambda1.std <- max(lambda1/(sdY)^2, 1e-7)
    lambda2.std <- max(lambda2/(sdY)^2, 1e-7)
    s <- (lambda1.std + lambda2.std)/N; alpha <- lambda1.std/(N*s)
    LL <- c(0, 10^seq(-2,2, .01) * s, max(1, 100*s))
    gnet <- glmnet(X.std, Y.std, alpha=alpha, standardize=FALSE, intercept=FALSE, lambda=LL, thresh = thresh, ...)
    ## the first item is the intercept
    Thetavec2 <- as.numeric(coef(gnet, s=s)[-1, ])
    ## Add pre-determined zeros back. Needs to be written in more efficient code
    if (length(zeroindx)>0){
        Thetavec <- rep(0, p^2)
        Thetavec[-zeroindx] <- Thetavec2
    } else {
        Thetavec <- Thetavec2
    }
    ## cut really small values
    Thetavec.cut <- replace(Thetavec, abs(Thetavec)<cutoff, 0)
    Theta <- matrix(Thetavec.cut, p, p)
    return(list(Theta=Theta, glmnet.obj=gnet, zeroindx=zeroindx, s=s, alpha=alpha, cutoff=cutoff))
}

## The big wrapper function
backfit <- function(Bhat, Astar, X0, Y, Ts, lambda1=.5, lambda2="auto", thresh=1e-11, cutoff=0.1, compensate=TRUE, gamma.approx=TRUE, coefs=c("integral", "sum"), level=5, eigen.prop.thresh=0.1^5, lambda2.min=0.1^5, prescreen="L1", prescreen.prop=.5,  ...){
    coefs <- match.arg(coefs)
    Binv <- solve(t(Bhat) %*% Bhat) %*% t(Bhat)
    ThetaHat <- Bhat %*% Astar %*% Binv
    ## linearization
    LIN <- linearize(Astar, Bhat, X0, Y, Ts, compensate=compensate, gamma.approx=gamma.approx, coefs=coefs, level=level, eigen.prop.thresh=eigen.prop.thresh, lambda2.min=lambda2.min)
    ## glmnet
    bfit <- lin2glmnet(LIN, ThetaHat, lambda1=lambda1, lambda2=lambda2, prescreen=prescreen, prescreen.prop=prescreen.prop, thresh=thresh, ...)
    return(bfit)
}
