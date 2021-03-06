% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/matrices.R
\name{Coord.to.Pointers}
\alias{Coord.to.Pointers}
\title{Convert a matrix defined by row and column indices to one
defined by a row- or column-oriented compression scheme.}
\usage{
Coord.to.Pointers(rows, cols, dims, triangle = TRUE, lower = TRUE,
  symmetric = FALSE, order = c("column", "row", "triplet"),
  index1 = TRUE)
}
\arguments{
\item{rows, cols}{row and column indices of non-zero elements}

\item{dims}{2-element vector for number of rows and columns.}

\item{triangle}{Is input intended to be a triangular (TRUE) or full
(FALSE) matrix. See details for how this argument is interpreted
for different values of \code{order}.}

\item{lower}{If \code{triangular} is true, this argument identifies
the input matrix as lower- or upper-triangular.  This argument is
ignored if \code{triangle} is FALSE.}

\item{symmetric}{If TRUE, and matrix is triangular, then the matrix
will be treated as symmetric, with only the triangular elements
provided.  If matrix is neither triangular nor symmetric, then
symmetric=TRUE will probably trigger an error.}

\item{order}{Determines the indexing/compression scheme for the
output matrix.  Use "triplet" to get row and column indices.
Defaults to the same class as M.}

\item{index1}{TRUE if using 1-based indexing.  FALSE for 0-based indexing.}
}
\value{
A list.  See Matrix.to.Pointers (no values are included in
return list).
}
\description{
Returns indices and pointers that define a sparse
Hessian in compressed format.  Inputs are the row and column indices.
}
\details{
\code{triangle} and \code{order} have the following interpretation:
\describe{
\item{triangle=TRUE}{Input \code{rows} and {cols} represent lower
or upper triangle of a matrix. If \code{order="symmetric"}, then
the output list will be for a full, symmetric matrix. Otherwise,
the output list will be for only the lower or upper triangle.  Any
elements outside of the specified triangle will trigger an error.}
\item{triangle=FALSE}{Input \code{rows} and {cols} represent a full
matrix. If that matrix is not symmetric, then
\code{order=="symmetric"} will trigger an error.}
If \code{symmetric=FALSE} and \code{order='triplet'}, the output
list should contain the same row and column indices as the input
list.}
}
\seealso{
Matrix.to.Pointers
}
