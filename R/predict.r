#' Prediction of an SFAR model
#'
#' Compute h-step-ahead prediction for an SFAR(1) model. Only the h-step
#' predicted function is returned, not the predictions for 1,2,...,h.
#' @param object an `sfar` object containing a fitted SFAR(1) model.
#' @param h number of steps ahead to predict.
#' @param ... Other parameters, not currently used.
#' @examples
#' # Generate Brownian motion noise
#' N <- 300 # the length of the series
#' n <- 200 # the sample rate that each function will be sampled
#' u <- seq(0, 1, length.out = n) # argvalues of the functions
#' d <- 45 # the number of bases
#' basis <- create.fourier.basis(c(0, 1), d) # the basis system
#' sigma <- 0.05 # the std of noise norm
#' Z0 <- matrix(rnorm(N * n, 0, sigma), nrow = n, nc = N)
#' Z0[, 1] <- 0
#' Z_mat <- apply(Z0, 2, cumsum) # N standard Brownian motion
#' Z <- smooth.basis(u, Z_mat, basis)$fd
#'
#' # Simulate random SFAR(1) data
#' kr <- function(x, y) {
#'   (2 - (2 * x - 1)^2 - (2 * y - 1)^2) / 2
#' }
#' s <- 5 # the period number
#' X <- rsfar(kr, s, Z)
#' plot(X)
#'
#' # SFAR(1) model parameter estimation:
#' Model1 <- sfar(X, seasonal = s, kn = 1)
#'
#' # Forecasting 3 steps ahead
#' fc <- predict(Model1, h = 3)
#' plot(fc)
#' @importFrom fda fd
#' @export
predict.sfar <- function(object, h, ...) {
  S <- object$seasonal
  Phi <- object$Phi
  basis <- object$kernel$sbasis
  N <- ncol(object$X$coefs)
  a <- h %/% S
  c0 <- h %% S
  Phi_a <- Phi
  for (i in seq(a - 1L)) {
    Phi_a <- Phi_a %*% Phi
  }
  c_1 <- Phi_a %*% as.matrix(object$X$coefs[, (N - S + c0)])
  return(fda::fd(c_1, basis))
}
