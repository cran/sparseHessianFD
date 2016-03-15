# NEWS file for sparseHessianFD package 

## VERSION 0.3.0

*  An even more major rewrite of the package.  All ACM code was
   removed, and replaced with original R/C++ implementations.

*  The sparseHessianFD class is now implemented as an R reference
  class, and not as an Rcpp module.  The `sparseHessianFD.new`
  function is deprecated.  Instead, use `sparseHessianFD` to
  initialize an object.  Initialization once again takes place in a single step.

*  The 'direct' computation method has been removed.  All computation
   uses the 'indirect' triangular substitution method.  The 'direct'
   argument in the initializer for the sparseHessianFD class is now
   deprecated, and remains solely for compatibility with older
   versions of the package.

*  There is a new vignette with a lot more detail about what the
   package does and how it works.

*  New matrix helper functions  `Matrix.to.Pointers` and
   `Coord.to.Pointers`.  The `Coord.to.Pattern.Matrix` function is deprecated.
   Use the `spMatrix` or `sparseMatrix` functions in the *Matrix*
   package instead.

*  With the removal of ACM-copyrighted code, this package is now
   licensed under the MPL 2.0.

## VERSION 0.2.0 (Jan. 28, 2015)

*  Essentially a complete re-write of the package.

*  New vignette, using a binary choice model as an example. Functions for this model are in binary.R.  Access sample simluated data with data(binary).

*  Documentation written using roxygen2

*  Added unit tests using testthat

*  Core class has been renamed sparseHessianFD. Construction and initialization are now two separate steps.

*  New function sparseHessianFD.new is a wrapper around the construction and initialization steps.

*  Deprecated functions in Version 0.2.0
    +  `new.sparse.hessian.object`.  Use `sparseHessianFD.new` instead.
    +  `Coord.to.Sym.Pattern.Matrix`. Use `Coord.to.Pattern.Matrix` instead, with option `symmetric=TRUE`.
    +  `Sym.CSC.to.Matrix`.  Use the `spMatrix` function in the *Matrix* package instead.
   


## VERSION 0.1.1 (Nov. 5, 2013)

*  Removed functions for sampling from, and computing the density of, a multivariate normal distribution.  These functions are now available in the *sparseMVN* package.


## VERSION 0.1.0 (Nov. 5, 2012)

*  Initial upload to CRAN.
