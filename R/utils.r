#' Create Block Diagonal matrix
#'
#' Get an kernel function k(s,t) corresponds to an integral operator and a
#' basis system, return the corresponding matrix of the operator with respect
#' to given basis system.
#' @param A a numeric vectors with values in [0,1].
#' @param k a  bivariate
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
  return(S1)
}

#' Create Block Diagonal matrix
#'
#' Get an kernel function k(s,t) corresponds to an integral operator and a
#' basis system, return the corresponding matrix of the operator with respect
#' to given basis system.
#' @param A a numeric vectors with values in [0,1].
#' @export
invsqrt <- function(A){
  if(!is.matrix(A)) stop("A is not a matrix") else
    if(nrow(A)!=ncol(A)) stop("A is not an sqare matrix") else
      if(!all(eigen(A)$values>0)) stop("A is not a positive definite matrix")
  n <- nrow(A)
  Q <- eigen(A)
  U <- Q$vectors
  lambda <- Q$values
  r <- length(lambda)
  M <- diag(1/sqrt(lambda), nrow = r, ncol = r)
  return(U %*% M %*% t(U))
}
