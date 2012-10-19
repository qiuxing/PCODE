## useful functions for high-dim ODE project
projected.fnorm <- function(V, W){
  ## projection Frobenius norm.
  Qv <- qr.Q(qr(V)); Qw <- qr.Q(qr(W))
  diffmat <- Qv %*% t(Qv) - Qw %*% t(Qw)
  return(norm(diffmat, "F")/sqrt(2))
}

## the procrustes mean on Graff(K,m) for pca/spca
graff.mean.pca <- function(pcalist){
  
}




grassmann.mean <- function(Ws){
  ## A simple way of computing centroid of points on Gr(m,n)
  N <- length(Ws)                       #sample size
  K <- ncol(Ws[[1]])                    #dim of eigen-space
  Qws <- sapply(Ws, function(W) qr.Q(qr(W)))

  mean.mat2 <- apply(sapply(Qws, function(Qw) W%*%t(Qw), simplify="array"), c(1,2), mean)
  eigen(apply(mean.mat2, c(1,2), mean))[["vectors"]][,1:K]
}



## karcher.mean <- function(Ws) {
##   ## the karcher mean/centroid on Gr(m,n)
##   N <- length(Ws)                       #sample size
##   objfun <- function(V){
##     sum(sapply(1:N, function(i) projected.fnorm(V, Ws[[i]])^2))
##   }
##   ...
## }


## align.grassmann <- function(V, W){

## }


## geodesic.norm <- function(V, W){
##   ## The geodesic distance on SOn
## }
