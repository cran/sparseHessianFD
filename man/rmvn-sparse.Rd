\name{mvn.sparse}
\alias{dmvn.sparse}
\alias{rmvn.sparse}
\title{Sampling and log density functions for multivariate normal distribution, given mean and the Cholesky decomposition of the precision or covariance matrix}
\description{Efficient sampling and density calculation from a multivariate normal, when the Cholesky decomposition of the precision matrix is already available.}
\usage{
dmvn.sparse(x, mu, L, prec=FALSE)
rmvn.sparse(n, mu, L, prec=FALSE)
}
\arguments{
	\item{x}{numeric matrix, where each row is an mvn sample.}
	\item{n}{number of draws}
	\item{mu}{mean (numeric vector)}
	\item{L}{Either the lower triangle of the Cholesky decomposition of the covariance or precision (inverse covariance) matrix, stored in a sparse format as defined in the Matrix package. Depends on value of prec.}
	\item{prec}{If TRUE, L is the lower triangle of the Cholesky of the precision matrix.  If FALSE, L is the lower triangle of the Cholesky decomposition of the covariance matrix.}
}
\details{
This function uses sparse matrix operations to sample from, or compute the log density of, a multivariate normal distribution  The user must compute the Cholesky decomposition first.  This must be a LOWER triangle; note that this is not the default behavior of the base \code{chol} function, which returns an upper triangle. For stability and efficiency, we recommend passing L as an object of class \code{"triangularMatrix"}, as defined in the \code{Matrix} package.
}
\value{
For \code{rmvn.sparse}, a matrix with n rows, with each draw in each row. The number of columns is equal to the length of mu. 

For \code{dmvn.sparse}, a vector with a length equal to the number of rows of \code{x}, containing the MVN log densities evaluated at each row of \code{x}.
}
