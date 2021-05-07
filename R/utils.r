#' Create block diagonal matrix
#'
#' @param A a numeric matrix forming each block.
#' @param k an integer value indicating the number of blocks.
#' @return Return a block diagonal matrix from the matrix A.
#'
#' @examples
#' Bdiag(matrix(1:4,2,2), 3)
#' @export
Bdiag <- function(A, k) {
  if(k < 1 | abs(k -round(k)) > 1e-6)
    stop("k must be a positive integer")
  m <- nrow(A)
  n <- ncol(A)
  B <- matrix(0, nrow=k*m, ncol=k*n)
  for(i in 0L:(k-1)) {
    B[i*m + seq(m), i*n + seq(n)] <- A
  }
  return(B)
}



#' Return inverse square root of square positive-definite matrix.
#'
#' @return inverse square root of square positive-definite matrix.
#' @param A a numeric matrix
#' @examples
#' require(Rsfar)
#' X <- Bdiag(matrix(1:4,2,2), 3)
#' invsqrt(t(X) %*% X)
#' @export
invsqrt <- function(A) {
  if (!is.matrix(A)) {
    stop("A is not a matrix")
  } else if (nrow(A) != ncol(A)) {
    stop("A is not an sqare matrix")
  }
  Q <- eigen(A)
  lambda <- Q$values
  if (!all(lambda > 0)) {
    stop("A is not a positive definite matrix")
  }
  U <- Q$vectors
  r <- length(lambda)
  M <- diag(1 / sqrt(lambda), nrow = r, ncol = r)
  return(U %*% M %*% t(U))
}
