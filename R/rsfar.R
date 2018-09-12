#' Simulation Functional Seasonal Auto regressive SARH(1) process.
#'
#' Simulation Functional Seasonal Auto regressive SARH(1) process,
#' where H is Hilbert space of L2[0,1].
#'
#' Seasonal autoregressive integrated moving average (ARIMA) models, which are considered
#' by Box and Jenkins (1976) and fully described by Brockwel and Davis (1991), are widely
#' applied. Another formulation of these models was proposed by Salas et al. (1982), who
#' developed an estimation procedure for ARIMA models with periodic coefficients.
#' The SARIMA model is developed from the ARMA model. This model is based on the application of ARMA models to transformed time series, where the seasonal and non-stationary
#' behavior has been eliminated. The SARIMA model, reflecting the feature of seasonal variation in time series, can be devided into simple and multiple models.
#' @return A sample of Functional time series from SARH(1) model of the calss'fd'.
#'
#' @param phi a kernel function corresponding to seasonal autoregressive operator. the.
#' @param seasonal an positive integer variable specifying the seasonality paprameter.
#' @param Z the functional noise object of the class 'fd'.
#'
#' @examples
#' require(sfar)
#' require(fda)
#' # Compute the standardaized constant of a kerel function with resoect to given HS norm.
#' gamma0 <- function(norm,kr){
#'   f <- function(x){
#'     g <- function(y) kr(x,y)^2
#'     return(integrate(g,0,1)$value) #return into f.
#'   }
#'   f <- Vectorize(f)
#'   A <- integrate(f,0,1)$value
#'   return(norm/A) #return into gamma.
#' }
#' #_________________________________________________
#' # Definition of parabolic integral kernel:
#' norm <- 0.99
#' kr <- function(x,y) 2-(2*x-1)^2-(2*y-1)^2
#' c0 <- gamma0(norm,kr)
#' phi <- function(x,y) c0*kr(x,y)
#'
#'
#' # Simulation a path from SFAR(1) model ________________________
#' N <- 300 # the length of the series
#' n <- 200 # the sample rate that each function will be sampled
#' s <- 5 # the period number
#' u <- seq(0, 1,length.out = n) # argvalues of the functions
#' d <- 15
#' # the basises number
#' basis <- create.fourier.basis(c(0,1),d) # the basis system
#'
#' sigma <- 0.05 # the std of noise norm
#' # Brownian motion noise generation:________________________________________
#' set.seed(10)
#' Z0 <- matrix(rnorm(N*n,0,sigma),nr=n,nc=N)
#' Z0[,1] <- 0
#' Z_mat <- apply(Z0,2,cumsum) # N standard Brownian motion
#' Z <- smooth.basis(u,Z_mat,basis)$fd
#' X <- rsfar(phi,s,Z)
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
