## mvn.R -- this file is part of sparseHessianFD, a contributed package
## for the R statistical programming platform.
##
## Copyright (C) 2012 Michael Braun

## this file demonstrates the use of the dmvn.sparse and rmvn.sparse
## functions for sampling from, and computing the log density of,
## a multivariate normal distribution, given a sparse lower triangular
## precision (inverse covariance) matrix.


rm(list=ls())
gc()

library(Matrix)
library(plyr)
library(sparseHessianFD)
library(mvtnorm)

set.seed(123)

N <- 50
k <- 2
n.draws <- 201
nu <- k+5

## simulate mean
mu <- rnorm(N*k)

## simulate a block-diagonal precision matrix
W <- rWishart(N, nu, diag(k))
L <- alply(W,3,function(x) return(x))
prec <- bdiag(L)

## sparse Cholesky
chol.prec <- t(chol(prec))

## recover covariance matrix
Sig.sparse <- crossprod(solve(chol.prec))
Sig.dense <- as(Sig.sparse, "matrix")

s1 <- system.time(x.sparse <- rmvn.sparse(n.draws, mu, chol.prec, prec=TRUE))
s2 <- system.time(x.dense <- rmvnorm(n.draws, mu, Sig.dense))

print(s1,units="secs")
print(s2,units="secs")

c1 <- diag(crossprod(x.sparse))/(N*k)
c2 <- diag(crossprod(x.dense))/(N*k)

y <- as(x.sparse,"matrix")

t1 <- system.time(d.sparse <- dmvn.sparse(y, mu, chol.prec, prec=TRUE))
t2 <- system.time(d.dense <- dmvnorm(y, mu, Sig.dense,log=TRUE))

print(t1,units="secs")
print(t2,units="secs")

print(all.equal(d.sparse, d.dense))

## check prec vs cov

chol.Sig <- t(chol(Sig.sparse))
x.cov <- rmvn.sparse(n.draws, mu, chol.Sig, prec=FALSE)
d.cov <- dmvn.sparse(y, mu, chol.Sig, prec=FALSE)
print(all.equal(d.sparse, d.dense))
c3 <- diag(crossprod(x.cov))/(N*k)
