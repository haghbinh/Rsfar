#' Estimation of an SFAR(1) Model
#'
#' Estimate a seasonal functional autoregressive (SFAR) model of order 1 for a given functional time series.
#' @return A matrix of size p*p.
#'
#' @param X a functional time series.
#' @param seasonal a positive integer variable specifying the seasonality parameter.
#' @param cpv a numeric with values in [0,1] which determines the cumulative proportion variance explained by the first kn eigencomponents.
#' @param kn an integer variable specifying the number of eigencomponents.
#' @param method a character string giving the method of estimation. The following values are possible:
#'   "MME" for Method of Moments, "ULSE" for Unconditional Least Square Estimation Method, and "KOE" for Kargin-Ontaski Estimation.
#' @param a a numeric with value in [0,1].
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
#'	 (2 - (2 * x - 1)^2 - (2 * y - 1)^2) / 2
#' }
#' s <- 5 # the period number
#' X <- rsfar(kr, s, Z)
#' plot(X)
#'
#' # SFAR(1) model parameter estimation:
#' Model1 <- sfar(X, seasonal = s, kn = 1)
#' @importFrom fda fd smooth.basis inprod pca.fd bifd
#' @export
sfar <- function(X, seasonal, cpv = 0.85, kn = NULL, method = c("MME","ULSE","KOE"), a = ncol(Coefs)^(-1/6)) {
  method <- match.arg(method)
  Coefs <- X$coefs
  d <- nrow(X$coefs)
  N <- ncol(Coefs)
  basis <- X$basis
  G <- inprod(basis, basis)
  Cl <- function(l) {
    s0 <- matrix(0L, nrow= d, ncol= d)
    for (i in (l + 1L):N) {
      s0 <- s0 + Coefs[, (i - l)] %*% t(Coefs[, i])
    }
    return(G %*% s0 / (N - l))
  }
  C0 <- Cl(0L)
  CS <- Cl(seasonal)
  Q <- pca.fd(fdobj = X, nharm = d, centerfns = FALSE)
  lambda0 <- Q$values
  if (is.null(kn)) kn <- which(cumsum(lambda0) / sum(lambda0) > cpv)[1L]
  lambda <- lambda0[1:kn]
  cpv1 <- sum(lambda) / sum(lambda0)
  nu <- Q$harmonics[1:kn] # FPC's components
  Vhat <- inprod(basis, nu)
  if (method == "MME") {
    V <- nu$coefs
    V0 <- inprod(basis, nu)
    La_inv <- diag(1/lambda,nrow=kn,ncol = kn)
    Phi <- V0 %*% La_inv %*% t(V) %*% CS %*% G # an d*d matrix corresponds to autoregressive operator.
  } else if (method == "ULSE") {
    X_scores <- as.matrix(Q$scores[, 1L:kn]) # an N*kn matrix of scores
    X0 <- matrix(c(X_scores[(seasonal + 1):N, ]), ncol = 1) # (N-s)kn by 1 matrix
    Z <- Bdiag(t(as.matrix(X_scores[1L, ])), kn)
    for (i in seq(2L, N - seasonal)) {
      Z1 <- Bdiag(t(as.matrix(X_scores[i, ])), kn)
      Z <- rbind(Z, Z1)
    }
    Phi_00 <- solve(t(Z) %*% Z) %*% t(Z) %*% X0
    Phi_hat0 <- matrix(c(Phi_00), byrow = TRUE, nrow= kn)
    Phi <- Vhat %*% Phi_hat0 %*% t(Vhat)
  } else if (method == "KOE") {
    Q00 <- eigen(C0)
    Ca <- C0 + diag(a, nrow = d, ncol = d)
    C12 <- invsqrt(Ca)
    B <- C12 %*% CS %*% t(CS) %*% C12
    Q <- eigen(B)
    lambda00 <- Q$values[1:kn]
    Va <- as.matrix(Q$vectors[, 1L:kn])
    V <- as.matrix(Q00$vectors[, 1L:kn])
    Phi <- C12 %*% G %*% Va %*% t(Va) %*% C12 %*% G %*% CS  # an d*d matrix corresponds to autoregressive operator.
  }
  structure(list(
    Phi = Phi,
    kernel = bifd(solve(G) %*% Phi, basis, basis),
    cpv = cpv1,
    CS = CS,
    a = a,
    X = X,
    tau = Q$values[1:kn],
    seasonal = seasonal,
    method = method
  ), class="sfar")
}

