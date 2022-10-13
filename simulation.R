
phi_inv <- function(z) {
  p <- c(exp(z), 1)
  p / sum(p)
}

#' @param X An N x p matrix giving design setting for N samples
#' @param beta K x p matrix giving effects for each of K species
lnm_simulator <- function(N, X, beta, sigma = 1, depth = 1000) {
  K <- nrow(X)
  Z <- rnorm(K - 1, X %*% t(beta), sigma)
  result <- matrix(nrow = N, ncol = K)
  
  for (i in seq_len(N)) {
    z <- rnorm(K - 1, X[i, ] %*% t(beta), sigma)
    result[i, ] <- rmultinom(1, depth, phi_inv(z))
  }
  
  result
}

composition_estimates <- function(Y) {
  Y / rowSums(Y)
}

#' @param X An N x p matrix giving design setting for N samples
#' @param beta K x p matrix giving effects for each of K species
nonparametric_simulator <- function(N, X, beta, p_hats) {
  exponential_tilt <- exp(X %*% t(beta))
  resample_ix <- sample(nrow(p_hats), replace = TRUE)
  p_sim <- p_hats[resample_ix, ] * exponential_tilt
  result <- matrix(nrow = N, ncol = K)
  
  for (i in seq_len(N)) {
    result[i, ] <- rmultinom(1, depth, p_sim[i, ])
  }
  
  result
}