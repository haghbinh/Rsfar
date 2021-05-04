# The corresponding matrix of an integral operator
#
# Find a kernel function k(s,t) corresponding to an integral operator and a
# basis system, and return the corresponding matrix of the operator with respect
# to the given basis system.
# @param u a numeric vector with values in [0,1].
# @param ker a function taking two inputs corrresponding to an integral operator.
# @param basis a functional basis object defining the basis.
# @return A numeric square matrix of the rank of basis system.
#' @importFrom fda smooth.basis inprod
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
