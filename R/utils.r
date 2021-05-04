#' Create block diagonal matrix
#'
#' @param A a numeric matrix forming each block.
#' @param k an integer value indicating the number of blocks.
#' @examples
#' Bdiag(matrix(1:4,2,2), 3)
#' @export
Bdiag <- function(A, k) {
  m <- nrow(A)
  n <- ncol(A)
  O <- matrix(0L, nrow = m, ncol = n)
  S1 <- A
  if (k > 1) {
    for (j in seq(2L, k))
      S1 <- cbind(S1, O)
    for (i in seq(2L, k)) {
      S2 <- O
      for (j in seq(2L, k)) {
        if (i == j) {
          S2 <- cbind(S2, A)
        } else {
          S2 <- cbind(S2, O)
        }
      }
      S1 <- rbind(S1, S2)
    }
  }
  else if (k < 1) {
    stop("k must be a positive integer")
  }
  return(S1)
}

#' Return inverse square root of square positive-definite matrix.
#'
#' @param A a numeric matrix
#' @examples
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
