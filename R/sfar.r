#' Estimation of SFAR(1) Model  Parameters.
#'
#' sfar is a function for SFAR(1) Model of a functional time series.
#' @return A matrix of size p*p.
#'
#' @param X a functional time series.
#' @param seasonal a positive integer variable specifying the seasonality paprameter.
#' @param cpv a numeric with values in [0,1] which determines the cumulative proportion variance explained by the first kn eigencomponents.
#' @param kn an integer variable specifying the number of eigencomponents.
#' @param method a character string giving the method of estimation. The following values are possible  :"MME" for Method of Moments and "ULSE" for Unconditional Least Square Estimation Method.
#' @param a a numeric with values in [0,1].
#' @examples
#' kr <- function(x,y) (2-(2*x-1)^2-(2*y-1)^2)/2
#' N <- 300 # the length of the series
#' n <- 200 # the sample rate that each function will be sampled
#' s <- 5 # the period number
#' u <- seq(0, 1,length.out = n) # argvalues of the functions
#' d <- 45 # the basises number
#' basis <- create.fourier.basis(c(0,1),d) # the basis system
#' sigma <- 0.05 # the std of noise norm
#'
#' # Brownian motion noise generation:
#' set.seed(10)
#' Z0 <- matrix(rnorm(N*n,0,sigma),nrow=n,nc=N)
#' Z0[,1] <- 0
#' Z_mat <- apply(Z0,2,cumsum) # N standard Brownian motion
#' Z <- smooth.basis(u,Z_mat,basis)$fd
#' X <- rsfar(kr,s,Z)
#' # SFAR(1) model parameters estimation:
#' Model1 <- sfar(X,seasonal=s,kn=1)
#' @export
sfar <- function(X, seasonal, cpv = 0.85, kn = NULL, method = "MME", a = ncol(Coefs)^(-1/6)) {
  Coefs <- X$coefs
  d <- nrow(X$coefs)
  N <- ncol(Coefs)
  basis <- X$basis
  G <- fda::inprod(basis, basis)
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
  out <- vector(mode = "list", length = 4L)
  Vhat <- inprod(basis, nu)
  if (method == "MME") {
    V <- nu$coefs
    V0 <- inprod(basis, nu)
    La_inv <- diag(1/lambda,nrow=kn,ncol = kn)
    Phi <- V0 %*% La_inv %*% t(V) %*% CS %*% G # an d*d matrix corresponds to autoregressive operator.
  }
  if (method == "ULSE") {
    X_scores <- as.matrix(Q$scores[, 1L:kn]) # an N*kn matrix of scores
    X0 <- matrix(c(X_scores[(seasonal + 1):N, ]), ncol = 1) # (N-s)kn by 1 matrix
    Z <- sfar:::Bdiag(t(as.matrix(X_scores[1L, ])), kn)
    for (i in seq(2L, N - seasonal)) {
      Z1 <- sfar:::Bdiag(t(as.matrix(X_scores[i, ])), kn)
      Z <- rbind(Z, Z1)
    }
    Phi_00 <- solve(t(Z) %*% Z) %*% t(Z) %*% X0
    Phi_hat0 <- matrix(c(Phi_00), byrow = TRUE, nrow= kn)
    Phi <- Vhat %*% Phi_hat0 %*% t(Vhat)
  }
  if (method == "KOE") {
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
    out <- list(Phi = Phi)
    out$kernel <- bifd(solve(G) %*% Phi, basis, basis)
    out$cpv <- cpv1
    out$CS <- CS
    out$a <- a
    out$X <- X
    out$tau <- Q$values[1:kn]
    out$seasonal <- seasonal
    out$method <- method
    class(out) <- "sfar"
    return(out)
}

