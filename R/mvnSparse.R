## matrices.R --   Part of the sparseHessianFD package for the R programming language.
##
## Copyright (C) 2012 Michael Braun


rmvn.sparse <- function(n, mu, L, prec=FALSE) {

  ## n = number of draws to collect.  Each column is a draw
  ## params is a list:
  ##   mu - mean of the MVN
  ##   L - lower-triangular Cholesky of the precision matrix (if prec==TRUE)
  ##        or of the covariance matrix (if prec==FALSE)
  
  k <- length(mu)
  x <- matrix(rnorm(n*k),k,n)  
  if (prec) {
    y <- as(solve(Matrix:::t(L),x),"matrix") ## L = chol(prec)
  } else {
    y <- as(L %*% x,"matrix") ## L = chol(cov)
  }
  
  res <- y+matrix(mu,k,n,byrow=FALSE)

  return(t(res))

}

dmvn.sparse <- function(x, mu, L, prec=FALSE) {

 
  ## INPUT:
  ##    x : point at which density is evaluated
  ##    params:  a list with the following elements:
  ##        mu : the mean of the MVN

  ##    L - lower-triangular Cholesky of the precision matrix (if prec==TRUE)
  ##        or of the covariance matrix (if prec==FALSE)

  ## OUTPUT:
  ##    log density of MVN evaluated at x
  
  k <- NCOL(x)
  M <- NROW(x)
  
  xmu <- t(x-matrix(mu,nrow=M, ncol=k,byrow=TRUE))

  ## something weird is going on with namespaces
  if (prec) {
    y <- Matrix:::t(L) %*% xmu
    log.dens <- -k*log(2*pi)/2 + sum(log(Matrix:::diag(L))) - Matrix:::colSums(y*y)/2
  } else {
    y <- solve(L, xmu)
    log.dens <- -k*log(2*pi)/2 - sum(log(Matrix:::diag(L))) - Matrix:::colSums(y*y)/2    
  }
    
  return(as.numeric(log.dens))
}

