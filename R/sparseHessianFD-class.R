## This file is part of the sparseHessianFD package
## Copyright (C) 2015-2017 Michael Braun
##
##' @title sparseHessianFD
##' @name sparseHessianFD
##' @description A reference class for computing sparse Hessians
##' @docType class
##' @field fn1 A closure for calling fn(x, ...).
##' @field gr1 A closure for calling gr(x, ...).
##' @field iRow,jCol Numeric vectors with row and column indices of
##' the non-zero elements in the lower triangle (including diagonal) of
##' the Hessian.
##' @field delta The perturbation amount for finite differencing of
##' the gradient to compute the Hessian. Defaults to 1e-07.
##' @field index1 TRUE if rows and cols use 1-based (R format)
##' indexing (FALSE for 0-based (C format) indexing.
##' @field complex TRUE if Hessian will be computed using the complex
##' step method, and FALSE (default) if using finite differences.
##' @field D raw finite differences (internal use only)
##' @field nvars Number of variables (length of x)
##' @field nnz Number of non-zero elements in the lower triangle of the Hessian.
##' @field ready TRUE if object has been initialized, and Hessian has
##' been partitioned.
##' @field idx,pntr Column indices and row pointers for non-zero
##' elements in lower triangle of the permuted Hessian.  Row-oriented
##' compressed storage.
##' @field colors A vector representation of the partitioning of the columns.
##' There are nvars elements, one for each column of the permuted
##' Hessian.  The value corresponds to the "color" for that column.
##' @field perm,invperm Permutation vector and its inverse
##' @details The sparseHessianFD function calls the initializer for the
##' sparseHessianFD class, and returns a sparseHessianFD object.
##' \preformatted{sparseHessianFD(x, fn, gr, rows, cols, delta, index1, complex, ...)}
##' The function, gradient and sparsity pattern are
##' declared as part of the initialization.
##'
##' Once initialized, the $hessian method will evaluate the Hessian at x.
##' \preformatted{
##' obj <- sparseHessian(x, fn, gr, rows, cols, ...)
##' obj$hessian(x)
##' }
##' For convenience, the class provides wrapper methods to the \code{fn} and \code{gr} functions that were specified in the initializer.
##' \preformatted{
##' obj$fn(P) ## wrapper to objective function
##' obj$gr(P) ## wrapper to gradient
##' obj$fngr(P) ## list of obj function and gradient
##' obj$fngrhs(P) ## list of obj function, gradient and Hessian.
##' }
##' \subsection{Arguments to initializer}{
##' \describe{
##' \item{x}{an vector at which the function, gradient and Hessian are initialized and tested.}
##' \item{fn, gr}{R functions that return the function value and
##' gradient, evaluated at x.}
##' \item{rows, cols}{Numeric vectors with row and column indices of
##' the non-zero elements in the lower triangle (including diagonal) of
##' the Hessian.}
##' \item{delta}{The perturbation amount for finite difference (or
##' complex step) of the gradient to compute the Hessian.  Defaults to 1e-07.}
##' \item{index1}{TRUE if rows and cols use 1-based (R format)
##' indexing (FALSE for 0-based (C format) indexing.}
##' \item{complex}{TRUE if Hessian will be computed using the complex
##' step method, and FALSE (default) if using finite differences. If
##' TRUE, both fn and gr must accept complex arguments and return complex values.}
##' \item{...}{other arguments to be passed to fn and gr.}
##' }
##' }
##' Other methods are described below.  Do not access any of the
##'     fields directly.  The internal structure is subject to change in future versions.
##' @examples
##' ## Log posterior density of hierarchical binary choice model. See vignette.
##' set.seed(123)
##' data("binary_small")
##' N <- length(binary[["Y"]])
##' k <- NROW(binary[["X"]])
##' T <- binary[["T"]]
##' P <- rnorm((N+1)*k)
##' priors <- list(inv.Sigma = rWishart(1,k+5,diag(k))[,,1],
##'                inv.Omega = diag(k))
##' true.hess <- binary.hess(P, binary, priors)
##' pattern <- Matrix.to.Coord(Matrix::tril(true.hess))
##' str(pattern)
##' obj <- sparseHessianFD(P, fn=binary.f, gr=binary.grad,
##'        rows=pattern[["rows"]], cols=pattern[["cols"]],
##'                       data=binary, priors=priors)
##' hs <- obj$hessian(P)
##' all.equal(hs, true.hess)
##'
##' f <- obj$fn(P) ## obj function
##' df <- obj$gr(P) ## gradient
##' fdf <- obj$fngr(P) ## list of obj function and gradient
##' fdfhs <- obj$fngrhs(P) ## list of obj function, gradient and Hessian.
##' @export sparseHessianFD
sparseHessianFD <-
    setRefClass("sparseHessianFD",
                fields=list(
                    fn1 = "function",
                    gr1 = "function",
                    iRow = "numeric",
                    jCol = "numeric",
                    delta = "numeric",
                    index1 = "logical",
                    complex = "logical",
                    D = "matrix",
                    nvars = "integer",
                    nnz = "integer",
                    ready = "logical",
                    idx = "integer",
                    pntr = "integer",
                    colors = "integer",
                    perm = "integer",
                    invperm = "integer"),
                methods = list(
                    initialize = function(x, fn, gr, rows, cols,
                                          delta=1e-7,
                                          index1 = TRUE, complex=FALSE, ...) {
                        "Initialize object with functions to compute the objective function and gradient (fn and gr), row and column indices of non-zero elements (rows and cols), an initial variable vector x at which fn and gr can be evaluated, a finite differencing parameter delta, flags for 0 or 1-based indexing (index1), whether the complex step method will be used, and other arguments (...) to be passed to fn and gr."


                        validate(fn, gr, rows, cols, x, delta, index1, complex, ...)
                        ww <- which(cols <= rows)

                        initFields(fn1 = function(y) fn(y, ...),
                                   gr1 = function(y) gr(y, ...),
                                   iRow = as.integer(rows[ww]),
                                   jCol = as.integer(cols[ww]),
                                   delta = delta,
                                   index1 = index1,
                                   complex = complex,
                                   nvars = length(x),
                                   ready = FALSE)
                        if (any(cols > rows)) {
                            warning("Some elements are in upper triangle, and will be ignored.")
                        }
                        nnz <<- length(iRow)

                        tmp <- sparseMatrix(i=iRow, j=jCol,
                                            index1=index1, symmetric=TRUE,
 ## fix for Matrix 1.3?     ##              giveCsparse=FALSE)
                                            repr = 'T')
                        perm <<- order(Matrix::rowSums(tmp), decreasing=TRUE)
                        invperm <<- invPerm(perm)

                        L <- tril(tmp[perm,perm])
                        ptr <- Matrix.to.Pointers(L, order="row", index1=index1)
                        idx <<- ptr[[1]]
                        pntr <<- ptr[[2]]

                        colors <<- coloring(L)
                        ncolors <- max(colors)+1

                        D <<- matrix(0, nvars, ncolors)
                        D[cbind(perm, colors+1)] <<- delta

                        ready <<- TRUE
                    },

                    partition = function() {
                        "Return the partitioning used to compute finite differences"
                        stopifnot(ready)
                        colors
                    },

                    validate = function(fn, gr, rows, cols, x,
                                        delta, index1, complex, ...) {

                        stopifnot(is.numeric(x),
                                  is.function(fn),
                                  is.function(gr),
                                  !is.null(rows),
                                  !is.null(cols),
                                  length(rows)==length(cols),
                                  all(is.finite(rows)),
                                  all(is.finite(cols)))

                        gradient <- gr(x,...)
                        val <- fn(x,...)
                        k <- length(x)

                        I1 <- as.integer(index1)
                        check.index1.row <- min(rows) >= I1 &
                          (max(rows) <= (k + I1 - 1))
                        check.index1.col <- min(cols) >= I1 &
                          (max(cols) <= (k + I1 - 1))

                        stopifnot(all(is.finite(gradient)),
                                  length(gradient)==k,
                                  check.index1.row,
                                  check.index1.col,
                                  is.finite(val)
                                  )
                        if(complex) {
                            h0 <- rep(0, k)
                            h0[1] <- delta * 1i
                            fc <- try(fn(x + h0, ...))
                            if(inherits(fc, "try-error"))
                                stop("function does not accept complex argument as required by complex step method")
                            if(!is.complex(fc))
                                stop("function does not return a complex value as required by complex step method")
                            gc <- try(gr(x + h0, ...))
                            if(inherits(gc, "try-error"))
                                stop("gradient does not accept complex argument as required by complex step method")
                            if(!all(is.complex(gc)))
                                stop("gradient does not return a complex value as required by complex step method")
                        }

                        dups <- which(duplicated(cbind(rows, cols)))
                        if (length(dups)>0) {
                            cat("Duplicates in sparsity pattern:\n",dups,"\n")
                            stop("sparseHessianFD does not allow duplicates")
                        }

                    },

                    fd = function(d, x, grad.x) {
                        gr1(x+d) - grad.x
                    },

                    fd_complex = function(d, x) {
                        Im(gr1(x + 1i * d))
                    },

                    hessian = function(x) {
                        "Return sparse Hessian, evaluated at x, as a dgCMatrix object."

                        stopifnot(ready)

                        if(complex) {
                            usingMethods(fd_complex)
                            Y2 <- apply(D, 2, fd_complex, x = x)
                        } else {
                            usingMethods(fd)
                            grad.x <- gr1(x)
                            Y2 <- apply(D, 2, fd, x = x, grad.x = grad.x)
                        }
                        Y <- Y2[perm,]
                        res <- subst(Y, colors,
                                     idx-index1, pntr-index1, delta, nvars, nnz)
                        return(res[invperm,invperm])
                    },

                    fn = function(x) {
                        "Return function value, evaluated at x: fn(x, ...)"
                        fn1(x)
                    },

                    gr = function(x) {
                        "Return gradient, evaluated at x:  gr(x,...)"
                        gr1(x)
                    },

                    fngr = function(x) {
                        "Return list of function value and gradient, evaluated at x"
                        list(fn = fn(x), gr = gr(x))
                    },

                    fngrhs = function(x) {
                        "Return list of function value, gradient, and Hessian, evaluated at x"
                        list(fn = fn(x),
                             gr = gr(x),
                             hessian=hessian(x))
                    },

                    pointers = function(out.index1=index1) {
                        "Return list with indices (idx) and pointers (pntr) for sparsity pattern of the compressed sparse Hessian.  Since the Hessian is symmetric, the indices and pointers for row-oriented and column-oriented storage patterns are the same."
                        stopifnot(ready)
                        list(idx = idx + out.index1 - index1,
                             pntr = pntr + out.index1 - index1)
                    },

                    get_nnz= function() {
                        "Return number of non-zero elements in lower triangle of Hessian"
                        stopifnot(ready)
                        nnz
                    },
                    get_nvars = function() {
                        "Return dimension (number of rows or columns) of Hessian"
                        stopifnot(ready)
                        nvars
                    },
                    get_perm = function() {
                        "Return integer vector of permutation used for computing Hessian"
                        stopifnot(ready)
                        perm
                    },
                    get_invperm = function() {
                        "Return integer vector of inverse of permutation used for computing Hessian"
                        stopifnot(ready)
                        invperm
                    },
                    get_pattern = function() {
                        "Return pattern matrix of lower triangle of Hessian"
                        stopifnot(ready)
                        sparseMatrix(i=iRow, j=jCol, index1=index1, symmetric=FALSE)
                    },
                    get_perm_pattern = function() {
                        "Return pattern matrix of lower triangle of *permuted* Hessian"
                        stopifnot(ready)
                        sparseMatrix(j=idx, p=pntr-index1, index1=index1, symmetric=FALSE)
                    }
                    )
                )
