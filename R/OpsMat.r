#' The corresponding matrix of an integral operator
#'
#' Get an kernel function k(s,t) corresponds to an integral operator and a
#' basis system, return the corresponding matrix of the operator with respect
#' to given basis system.
#' @param u a numeric vectors with values in [0,1].
#' @param ker a  bivariate numeric corresponds to an integral operator.
#' @param basis a functional basis object defining the basis.
#' @return A numeric square matrix of the rank of basis system.
#' @export
OpsMat <- function(u, ker, basis) {
  n <- length(u)
  K_mat <- outer(u, u, FUN = ker)
  K_t <- smooth.basis(u, K_mat, basis)$fd # the kernel function convert
  A <- inprod(K_t, basis) # An n*d matrix.
  K <- smooth.basis(u, A, basis)$fd # d fd object.
  B <- inprod(K, basis) # An d*d matrix.
  G <- inprod(basis, basis) # Gram matrix
  return(solve(G) %*% B)
}
