# Integration scores -------------------------------------------------------

#' Score an integrated embedding dimension by dimension
#' @export
IntegrationScore <- function(seurat.obj, reduction, ref.batch,
                             ndims = NULL,
                             batch = "batch",
                             K = 5, n.cores = 10, seed = 42) {
  set.seed(seed)

  emb <- Seurat::Embeddings(seurat.obj, reduction)
  if (is.null(ndims)) {
    ndims <- ncol(emb)
  }

  X <- emb[, seq_len(ndims), drop = FALSE]
  batch_label <- .get_batch_vector(seurat.obj, batch)
  batches <- unique(batch_label)
  query_batches <- batches[batches != ref.batch]

  batchE_by_batch <- matrix(NA_real_, nrow = ncol(X), ncol = length(query_batches))
  celltypeE_by_batch <- matrix(NA_real_, nrow = ncol(X), ncol = length(query_batches))

  message(paste0("Computing integration score on ", length(query_batches), " batches...."))

  for (i in seq_along(query_batches)) {
    qb <- query_batches[i]
    message(paste0("On batch: ", qb))

    res <- pbmcapply::pbmclapply(
      seq_len(ncol(X)),
      function(d) {
        set.seed(seed)
        tryCatch({
          x <- X[, d]
          qs <- stats::quantile(x, probs = c(0.25, 0.75), na.rm = TRUE)
          iqr <- qs[2] - qs[1]
          l <- qs[1] - 1.5 * iqr
          h <- qs[2] + 1.5 * iqr
          ids <- (l <= x & x <= h) | (x <= stats::quantile(x, 0.97) & -x <= stats::quantile(-x, 0.97))

          id_ref <- batch_label == ref.batch
          id_query <- batch_label == qb
          X_ref <- x[id_ref & ids]
          X_query <- x[id_query & ids]

          fit_ref <- EMmix_normal(X_ref, K_max = K)
          reweight_query <- Reweight_normal(X_query, fit_ref)
          fit_query <- EMmix_normal(X_query, K_max = K)

          list(
            batch = -(fit_query$loglikelihood - reweight_query$loglik_reweighted) / length(X_query),
            celltype = 1 / length(query_batches) *
              (fit_ref$loglikelihood - sum(stats::dnorm(X_ref, mean = mean(X_ref), sd = stats::sd(X_ref), log = TRUE))) / length(X_ref) +
              (fit_query$loglikelihood - sum(stats::dnorm(X_query, mean = mean(X_query), sd = stats::sd(X_query), log = TRUE))) / length(X_query)
          )
        }, error = function(e) {
          NA_real_
        })
      },
      mc.cores = n.cores
    )

    batchE_by_batch[, i] <- vapply(res, function(x) if (is.list(x)) x[[1]] else NA_real_, numeric(1))
    celltypeE_by_batch[, i] <- vapply(res, function(x) if (is.list(x) && length(x) >= 2) x[[2]] else NA_real_, numeric(1))
  }

  effects <- data.frame(
    BatchAlignment = rowMeans(batchE_by_batch, na.rm = TRUE),
    CellIdentity = rowSums(celltypeE_by_batch, na.rm = TRUE) / (length(query_batches) + 1)
  )
  rownames(effects) <- seq_len(ndims)

  seurat.obj@misc$integration_effects[[reduction]] <- effects
  seurat.obj
}

#' Score an integrated embedding using a joint multivariate GMM
#' @export
IntegrationScore.Joint <- function(seurat.obj, reduction, ref.batch,
                                   ndims = NULL,
                                   batch = "batch",
                                   K = 20, seed = 42, n.cores = 1) {
  set.seed(seed)

  emb <- Seurat::Embeddings(seurat.obj, reduction)
  if (is.null(ndims)) {
    ndims <- ncol(emb)
  }

  X <- emb[, seq_len(ndims), drop = FALSE]
  batch_label <- .get_batch_vector(seurat.obj, batch)
  batches <- unique(batch_label)
  query_batches <- batches[batches != ref.batch]

  message(paste0("Computing integration score on ", 1 + length(query_batches), " batches...."))
  message("Pre-calculating Reference model...")

  id_ref <- which(batch_label == ref.batch)
  X_ref <- X[id_ref, , drop = FALSE]
  fit_ref <- EMmix_normal_joint(X_ref, K, model = "VVI")

  n_ref <- nrow(X_ref)
  p <- ncol(X_ref)
  mu_ref <- colMeans(X_ref)
  Xc_ref <- sweep(X_ref, 2, mu_ref, "-")
  vars_ref_mle <- colSums(Xc_ref^2) / n_ref

  mvn1_ref_vals <- mvtnorm::dmvnorm(
    X_ref,
    mean = mu_ref,
    sigma = diag(vars_ref_mle, p),
    log = TRUE
  )

  ref_loglik_diff <- (fit_ref$loglik - sum(mvn1_ref_vals)) / nrow(X_ref)

  if (n.cores > 1) {
    results <- pbmcapply::pbmclapply(query_batches, function(qb) {
      .score_joint_query_batch(qb, X, batch_label, fit_ref)
    }, mc.cores = n.cores)
  } else {
    results <- lapply(query_batches, function(qb) {
      .score_joint_query_batch(qb, X, batch_label, fit_ref)
    })
  }

  res_mat <- do.call(rbind, results)
  batchE_by_batch <- res_mat[, "bE"]
  cE_parts <- res_mat[, "cE_part"]

  n_query <- length(query_batches)
  avg_celltypeE <- ref_loglik_diff + sum(cE_parts)

  effects <- data.frame(
    BatchAlignment = mean(batchE_by_batch, na.rm = TRUE),
    CellIdentity = avg_celltypeE / (n_query + 1)
  )

  seurat.obj@misc$integration_effects[[reduction]] <- effects
  seurat.obj
}

.score_joint_query_batch <- function(qb, X, batch_label, fit_ref) {
  id_query <- which(batch_label == qb)
  X_query <- X[id_query, , drop = FALSE]

  n <- nrow(X_query)
  p <- ncol(X_query)

  fit_query <- EMmix_normal_joint(X_query, length(fit_ref$pis), model = "VVI")
  reweight_query <- Reweight_normal_joint(X_query, fit_ref)

  mu_query <- colMeans(X_query)
  Xc_query <- sweep(X_query, 2, mu_query, "-")
  vars_query_mle <- colSums(Xc_query^2) / n

  mvn1_query_vals <- mvtnorm::dmvnorm(
    X_query,
    mean = mu_query,
    sigma = diag(vars_query_mle, p),
    log = TRUE
  )

  bE <- -(fit_query$loglik - reweight_query$loglik_reweighted) / nrow(X_query)
  cE_part <- (fit_query$loglik - sum(mvn1_query_vals)) / nrow(X_query)

  c(bE = bE, cE_part = cE_part)
}
