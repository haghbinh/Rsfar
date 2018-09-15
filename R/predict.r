#' Prediction of an SFAR model
#'
#' The h step ahead prediction for an SFAR(1) model.
#' @param object a fit sfar object.
#' @param h number of steps ahead at which to predict.
#' @param ... other.
#' @export
predict.sfar <- function(object, h, ...) {
  S <- object$seasonal
  Phi <- object$Phi
  basis <- object$kernel$sbasis
  N <- ncol(object$X$coefs)
  a <- h %/% S
  c0 <- h %% S
  Phi_a <- Phi
  for (i in 1:(a-1)) Phi_a <- Phi_a %*% Phi
  c_1 <- Phi_a %*% as.matrix(object$X$coefs[,(N-S+c0)])
  return(fda::fd(c_1,basis))
}
