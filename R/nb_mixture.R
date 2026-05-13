# Negative-binomial mixture model -----------------------------------------
#' @importFrom stats dnbinom dnorm
#' @importFrom utils getFromNamespace tail
NBLogpmf <- function(x, mu, size, l) {
  stats::dnbinom(x, size = size, mu = l * mu, log = TRUE)
}

Estep <- function(X, l_vec, pis, mus, size) {
  n <- length(X)
  K <- length(pis)
  gamma <- matrix(0, n, K)

  for (k in seq_len(K)) {
    gamma[, k] <- pis[k] * exp(NBLogpmf(X, mus[k], size, l_vec))
  }

  rs <- rowSums(gamma)
  gamma <- gamma / rs
  gamma
}

MstepComponent <- function(X, l_vec, weights) {
  sum(weights * (X / l_vec)) / sum(weights)
}

PrecomputeDispersionData <- function(X, l_vec, mus, gamma) {
  n <- length(X)
  K <- length(mus)

  list(
    X_j = rep(X, times = K),
    mu_j = rep(l_vec, times = K) * rep(mus, each = n),
    w_j = as.vector(gamma)
  )
}

NegLogLikelihoodFast <- function(size, X_j, mu_j, w_j) {
  if (size <= 0) return(Inf)
  ll <- w_j * stats::dnbinom(X_j, size = size, mu = mu_j, log = TRUE)
  -sum(ll)
}

NegLogLikGradient <- function(size, X_j, mu_j, w_j) {
  if (size <= 0) return(Inf)

  term <- w_j * (
    base::digamma(X_j + size) - base::digamma(size) +
      log(size) - log(size + mu_j) +
      1 - (X_j + size) / (size + mu_j)
  )

  -sum(term)
}

MstepOptim <- function(X, l_vec, mus, gamma, size_init = 10) {
  pre <- PrecomputeDispersionData(X, l_vec, mus, gamma)
  res <- stats::optim(
    par = size_init,
    fn = NegLogLikelihoodFast,
    gr = NegLogLikGradient,
    X_j = pre$X_j,
    mu_j = pre$mu_j,
    w_j = pre$w_j,
    method = "L-BFGS-B",
    lower = 1e-5
  )
  res$par
}

EMmix <- function(X, l_vec, K, max_iter = 10, tol = 1e-4) {
  keep <- !is.na(X) & !is.na(l_vec) & l_vec > 0
  X <- X[keep]
  l_vec <- l_vec[keep]

  n <- length(X)
  if (n == 0) {
    return(list(
      pis = rep(NA_real_, K),
      mus = rep(NA_real_, K),
      size = NA_real_,
      gamma = matrix(NA_real_, 0, K),
      loglikelihood_trace = NA_real_,
      loglikelihood = NA_real_
    ))
  }

  pis <- rep(1 / K, K)
  x_scaled <- X / l_vec
  if (length(unique(x_scaled)) <= 1) {
    mus <- rep(mean(x_scaled), K)
  } else {
    mus <- seq(min(x_scaled), max(x_scaled), length.out = K + 2)[-c(1, K + 2)]
  }
  mus <- pmax(mus, .Machine$double.eps)
  size <- 1

  loglikelihood_trace <- numeric(max_iter)

  for (iter in seq_len(max_iter)) {
    gamma <- Estep(X, l_vec, pis, mus, size)

    logL_i <- numeric(n)
    for (i in seq_len(n)) {
      component_probs <- pis * stats::dnbinom(
        X[i],
        size = size,
        mu = l_vec[i] * mus
      )
      logL_i[i] <- log(sum(component_probs))
    }
    loglikelihood_trace[iter] <- sum(logL_i)

    pis_query <- colSums(gamma) / n
    mus_query <- numeric(K)
    for (k in seq_len(K)) {
      mus_query[k] <- MstepComponent(X, l_vec, gamma[, k])
    }
    mus_query <- pmax(mus_query, .Machine$double.eps)

    size_query <- MstepOptim(X, l_vec, mus_query, gamma, size_init = size)

    delta <- sum(abs(pis_query - pis)) + sum(abs(mus_query - mus)) + abs(size_query - size)
    pis <- pis_query
    mus <- mus_query
    size <- size_query

    if (is.finite(delta) && delta < tol) {
      loglikelihood_trace <- loglikelihood_trace[seq_len(iter)]
      break
    }
  }

  list(
    pis = pis,
    mus = mus,
    size = size,
    gamma = gamma,
    loglikelihood_trace = loglikelihood_trace,
    loglikelihood = tail(loglikelihood_trace, 1)
  )
}

#' Reweight an NB mixture on query data
#'
#' Keeps reference component means and dispersion fixed, and updates only the
#' query mixture proportions by a small number of EM-style iterations.
Reweight <- function(X_query, l_query, fit_ref) {
  # Extract ref model
  mus     <- fit_ref$mus
  size    <- fit_ref$size
  K       <- length(fit_ref$pis)
  n_query <- length(X_query)
  # Reweight using E-step on query data
  gamma_query <- matrix(0, n_query, K)
  for (k in 1:K) {
    gamma_query[, k] <- dnbinom(X_query, size = size, mu = l_query * mus[k])
  }
  gamma_query <- gamma_query / rowSums(gamma_query)
  pis_query <- colSums(gamma_query) / n_query
  # Compute log-likelihood using pis_query
  loglik <- 0
  for (i in 1:n_query) {
    loglik <- loglik + log(sum(pis_query * dnbinom(X_query[i], size = size, mu = l_query[i] * mus)))
  }
  
  list(
    loglik_reweighted = loglik,
    pis_query = pis_query,
    mus_reweighted = mus,
    size_reweighted = size
  )
}
