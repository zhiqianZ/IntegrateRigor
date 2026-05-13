# Gaussian mixture model helpers ------------------------------------------
#' @importFrom mclust Mclust mclustBIC
EMmix_normal <- function(x, K_max, cl = NULL) {
  x <- x[!is.na(x)]

  if (!is.null(cl)) {
    gmm <- mclust::Mclust(
      x,
      G = length(unique(cl)),
      initialization = list(class = cl)
    )
  } else {
    gmm <- mclust::Mclust(x, G = K_max)
  }

  param <- gmm$parameters
  K_fit <- length(param$pro)

  if (gmm$modelName == "E") {
    variance <- rep(sqrt(param$variance$sigmasq), K_fit)
  } else {
    variance <- sqrt(param$variance$sigmasq)
  }

  list(
    pis = param$pro,
    mus = param$mean,
    sizes = pmax(variance, .Machine$double.eps),
    gamma = gmm$z,
    loglikelihood = gmm$loglik,
    mclust = gmm
  )
}

#' Reweight a 1D Gaussian mixture on query data
#'
#' Keeps reference component means and standard deviations fixed, and updates
#' only the query mixture proportions.
Reweight_normal <- function(X_query, fit_ref) {
  # Extract ref model
  mus     <- fit_ref$mus
  sizes    <- fit_ref$sizes
  K       <- length(fit_ref$pis)
  n_query <- length(X_query)
  log_gamma_query <- matrix(NA_real_, n_query, K)
  for (k in 1:K) {
    log_gamma_query[, k] <- dnorm(X_query, mean = mus[k], sd = sizes[k], log = TRUE)
  }
  row_max <- apply(log_gamma_query, 1, max)
  log_gamma_query_centered <- log_gamma_query - row_max
  gamma_query <- exp(log_gamma_query_centered) 
  gamma_query <- gamma_query / rowSums(gamma_query)
  
  pis_query <- colSums(gamma_query) / n_query
  
  log_sum_exp <- function(z) {
    m <- max(z)
    m + log(sum(exp(z - m)))
  }
  
  loglik_vec <- vapply(
    X_query,
    function(x) {
      log_sum_exp(
        log(pis_query) + dnorm(x, mean = mus, sd = sizes, log = TRUE)
      )
    },
    numeric(1L)
  )
  
  total_loglik <- sum(loglik_vec)
  
  
  list(
    loglik_reweighted = total_loglik,
    pis_query = pis_query,
    mus_reweighted = mus,
    size_reweighted = sizes
  )
}

EMmix_normal_joint <- function(X, K_max, model = "VVI") {
  X <- as.matrix(X)
  X <- X[stats::complete.cases(X), , drop = FALSE]

  gmm <- mclust::Mclust(X, G = K_max, modelNames = model)
  param <- gmm$parameters

  list(
    pis = param$pro,
    mus = param$mean,
    covs = param$variance$sigma,
    gamma = gmm$z,
    loglik = gmm$loglik,
    mclust = gmm
  )
}

Reweight_normal_joint <- function(X_query, fit_ref) {
  X_query <- as.matrix(X_query)
  X_query <- X_query[stats::complete.cases(X_query), , drop = FALSE]

  pis <- fit_ref$pis
  mus <- fit_ref$mus
  covs <- fit_ref$covs

  K <- length(pis)
  n <- nrow(X_query)

  if (n == 0) {
    return(list(loglik_reweighted = NA_real_, pis_query = rep(NA_real_, K)))
  }

  log_gamma_query <- matrix(NA_real_, n, K)
  for (k in seq_len(K)) {
    log_gamma_query[, k] <- mvtnorm::dmvnorm(
      X_query,
      mean = mus[, k],
      sigma = covs[, , k],
      log = TRUE
    )
  }

  row_max <- apply(log_gamma_query, 1, max)
  gamma_centered <- exp(log_gamma_query - row_max)
  gamma_query <- gamma_centered / rowSums(gamma_centered)
  pis_query <- colSums(gamma_query) / n

  loglik_vec <- vapply(seq_len(n), function(i) {
    ll_comp <- numeric(K)
    for (k in seq_len(K)) {
      ll_comp[k] <- log(pis_query[k]) +
        mvtnorm::dmvnorm(
          X_query[i, , drop = FALSE],
          mean = mus[, k],
          sigma = covs[, , k],
          log = TRUE
        )
    }
    .log_sum_exp(ll_comp)
  }, numeric(1))

  list(
    loglik_reweighted = sum(loglik_vec),
    pis_query = pis_query
  )
}
