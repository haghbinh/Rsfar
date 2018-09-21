#' Simulation of a Seasonal Functional Autoregressive SFAR(1) process.
#'
#' Simulation of a SFAR(1) process on a Hilbert space of L2[0,1].
#'
#' @return A sample of functional time series from a SFAR(1) model of the class `fd`.
#'
#' @param phi a kernel function corresponding to the seasonal autoregressive operator.
#' @param seasonal a positive integer variable specifying the seasonal period.
#' @param Z the functional noise object of the class 'fd'.
#'
#' @examples
#' # Set up Brownian motion noise process
#' N <- 300 # the length of the series
#' n <- 200 # the sample rate that each function will be sampled
#' u <- seq(0, 1, length.out = n) # argvalues of the functions
#' d <- 15 # the number of basis functions
#' basis <- create.fourier.basis(c(0, 1), d) # the basis system
#' sigma <- 0.05 # the stdev of noise norm
#' Z0 <- matrix(rnorm(N * n, 0, sigma), nr = n, nc = N)
#' Z0[, 1] <- 0
#' Z_mat <- apply(Z0, 2, cumsum) # N standard Brownian motion
#' Z <- smooth.basis(u, Z_mat, basis)$fd
#'
#' # Compute the standardized constant of a kernel function with respect to a given HS norm.
#' gamma0 <- function(norm, kr) {
#'   f <- function(x) {
#'     g <- function(y) {
#'       kr(x, y)^2
#'     }
#'     return(integrate(g, 0, 1)$value)
#'   }
#'   f <- Vectorize(f)
#'   A <- integrate(f, 0, 1)$value
#'   return(norm / A)
#' }
#' # Definition of parabolic integral kernel:
#' norm <- 0.99
#' kr <- function(x, y) {
#'   2 - (2 * x - 1)^2 - (2 * y - 1)^2
#' }
#' c0 <- gamma0(norm, kr)
#' phi <- function(x, y) {
#'   c0 * kr(x, y)
#' }
#'
#' # Simulating a path from an SFAR(1) process
#' s <- 5 # the period number
#' X <- rsfar(phi, s, Z)
#' plot(X)
#' @export
rsfar <- function(phi, seasonal, Z) {
  u <- seq(Z$basis$rangeval[1], Z$basis$rangeval[2], by = 0.01)
  basis <- Z$basis
  k_mat <- OpsMat(u, ker = phi, basis)
  Z_coef <- Z$coefs
  d <- nrow(Z_coef)
  N <- ncol(Z_coef)
  X_coef <- Z_coef
  for (i in (seasonal + 1L):N) X_coef[, i] <- X_coef[, (i - seasonal)] %*% k_mat + Z_coef[, i]
  return(fd(X_coef, basis))
}
