# Create Block Diagonal matrix
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

# Setting wireframe plot panel in lattice:
panel.3d.levelplot <-
  function(x, y, z, rot.mat, distance, zlim.scaled, at,
             drape = TRUE, shade = FALSE, ...) {
    panel.3dwire(x, y, z,
      rot.mat = rot.mat, distance = distance,
      zlim.scaled = zlim.scaled, at = at,
      drape = FALSE, shade = shade, ...
    )
    zrng <- 0.001 * diff(zlim.scaled) ## vertical range of 2nd surface (>0)
    z.scaled <- (z - min(z)) / diff(range(z))
    at.scaled <- (at - min(z)) / diff(range(z))
    new.z <- zlim.scaled[2] + zrng * (z.scaled - 1)
    new.at <- zlim.scaled[2] + zrng * (at.scaled - 1)
    panel.3dwire(x, y, new.z,
      at = new.at,
      col = "transparent",
      rot.mat = rot.mat, distance = distance,
      shade = FALSE,
      drape = TRUE,
      zlim.scaled = zlim.scaled,
      alpha = 0.8,
      ...
    )
  }

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
