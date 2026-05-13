# Parameter search ---------------------------------------------------------

#' Search integration parameters using marginal integration scores
#' @export
IntegrateRigor.ParameterS <- function(seurat.obj, parameter.df, method,
                                      batch = "batch", ndims = 30,
                                      ndims.score = min(30, ndims),
                                      ref.batch = NULL, use.prev.ref = TRUE,
                                      force.run = TRUE,
                                      K = 5, n.cores = 5, seed = 42,
                                      return.all = TRUE,
                                      use.batch.stable.genes = TRUE,
                                      verbose = FALSE,
                                      subsample = NULL,
                                      run.optimal = TRUE,
                                      default.assay = "RNA") {
  .parameter_search_impl(
    seurat.obj = seurat.obj,
    parameter.df = parameter.df,
    method = method,
    method.name =  tolower(deparse(substitute(method))),
    score_fun = IntegrationScore,
    score_mode = "marginal",
    batch = batch,
    ndims = ndims,
    ndims.score = ndims.score,
    ref.batch = ref.batch,
    use.prev.ref = use.prev.ref,
    force.run = force.run,
    K = K,
    n.cores = n.cores,
    seed = seed,
    return.all = return.all,
    use.batch.stable.genes = use.batch.stable.genes,
    verbose = verbose,
    subsample = subsample,
    run.optimal = run.optimal,
    default.assay = default.assay
  )
}

#' Search integration parameters using a joint multivariate GMM score
#' @export
IntegrateRigor.ParameterS.Joint <- function(seurat.obj, parameter.df, method,
                                            batch = "batch", ndims = 30,
                                            ndims.score = min(30, ndims),
                                            ref.batch = NULL, use.prev.ref = TRUE,
                                            force.run = TRUE,
                                            K = 20, seed = 42, n.cores = 1,
                                            return.all = TRUE,
                                            use.batch.stable.genes = TRUE,
                                            verbose = FALSE,
                                            subsample = NULL,
                                            run.optimal = TRUE,
                                            default.assay = "RNA") {
  .parameter_search_impl(
    seurat.obj = seurat.obj,
    parameter.df = parameter.df,
    method = method,
    method.name =  tolower(deparse(substitute(method))),
    score_fun = IntegrationScore.Joint,
    score_mode = "joint",
    batch = batch,
    ndims = ndims,
    ndims.score = ndims.score,
    ref.batch = ref.batch,
    use.prev.ref = use.prev.ref,
    force.run = force.run,
    K = K,
    n.cores = n.cores,
    seed = seed,
    return.all = return.all,
    use.batch.stable.genes = use.batch.stable.genes,
    verbose = verbose,
    subsample = subsample,
    run.optimal = run.optimal,
    default.assay = default.assay
  )
}

.parameter_search_impl <- function(seurat.obj, parameter.df, method, score_fun,
                                   method.name,
                                   score_mode = c("marginal", "joint"),
                                   batch = "batch", ndims = 30,
                                   ndims.score = ndims,
                                   ref.batch = NULL, use.prev.ref = TRUE,
                                   force.run = TRUE,
                                   K = 5, n.cores = 5, seed = 42,
                                   return.all = TRUE,
                                   use.batch.stable.genes = TRUE,
                                   verbose = FALSE,
                                   subsample = NULL,
                                   run.optimal = TRUE,
                                   default.assay = "RNA") {
  score_mode <- match.arg(score_mode)

  BAS <- numeric(nrow(parameter.df))
  CIS <- numeric(nrow(parameter.df))
  nms <- character(nrow(parameter.df))
  reduction_names <- character(nrow(parameter.df))

  if (!is.null(subsample)) {
    stopifnot(subsample > 0, subsample <= 1)
    set.seed(seed)
    obj.raw <- seurat.obj
    id <- sample(seq_len(ncol(seurat.obj)), round(subsample * ncol(seurat.obj)))
    seurat.obj <- seurat.obj[, id]
  } else {
    batch_vec0 <- .get_batch_vector(seurat.obj, batch)
    if (ncol(seurat.obj) > 20000 || length(unique(batch_vec0)) > 6) {
      message("Executing individual integration method may be time-consuming due to the large number of batches or cells. Consider using the subsample parameter to accelerate the process.")
    }
  }

  if (is.null(ref.batch)) {
    if (use.prev.ref && !is.null(seurat.obj@misc$reference)) {
      ref.batch <- seurat.obj@misc$reference
      message(paste0("Using previously recorded reference batch: ", ref.batch))
    } else {
      message("Finding reference batch")
      ref.batch <- FindReference(seurat.obj, batch)
      seurat.obj@misc$reference <- ref.batch
      message(paste0("Setting reference batch as: ", ref.batch))
    }
  }

  for (i in seq_len(nrow(parameter.df))) {
    values <- parameter.df[i, , drop = FALSE]
    param_list <- .row_to_param_list(values)

    msg <- paste(
      vapply(names(param_list), function(nm) {
        val <- param_list[[nm]]
        paste0(nm, " = ", if (is.null(val)) "NULL" else val)
      }, character(1)),
      collapse = ", "
    )
    message("Running with ", msg)

    param_suffix <- .make_param_suffix(param_list)
    reduction_prefix <- if (use.batch.stable.genes) "integrated.bsg." else "integrated."
    new.reduction <- paste0(reduction_prefix, method.name, ".", param_suffix)

    if (!new.reduction %in% names(seurat.obj@reductions) || force.run) {
      set.seed(seed)
      seurat.obj <- rlang::inject(
        Integration(
          seurat.obj = seurat.obj,
          batch = batch,
          method = method,
          new.reduction = new.reduction,
          ndims = ndims,
          verbose = verbose,
          default.assay = default.assay,
          !!!param_list
        )
      )
      message("Integration Done!")
    }

    message("Starting calculating score...")
    seurat.obj <- score_fun(
      seurat.obj,
      reduction = new.reduction,
      batch = batch,
      ref.batch = ref.batch,
      K = K,
      n.cores = n.cores,
      seed = seed,
      ndims = ndims.score
    )

    effects <- seurat.obj@misc$integration_effects[[new.reduction]]
    if (score_mode == "joint") {
      BAS[i] <- effects[1, "BatchAlignment"]
      CIS[i] <- effects[1, "CellIdentity"]
    } else {
      BAS[i] <- mean(effects[, "BatchAlignment"], na.rm = TRUE)
      CIS[i] <- mean(effects[, "CellIdentity"], na.rm = TRUE)
    }

    if (!is.null(subsample)) {
      obj.raw@misc$integration_effects[[new.reduction]] <- effects
    }

    nms[i] <- param_suffix
    reduction_names[i] <- new.reduction
  }

  names(BAS) <- nms
  names(CIS) <- nms
  IS <- BAS + CIS
  names(IS) <- nms

  best_idx <- which.max(IS)
  best_name <- names(IS)[best_idx]

  search.df <- data.frame(
    Parameter = nms,
    Reduction = reduction_names,
    BatchAlignment = as.numeric(BAS),
    CellIdentity = as.numeric(CIS),
    IntegrationScore = as.numeric(IS),
    row.names = NULL
  )

  message("Optimal parameter: ", best_name)

  optimal.reduction <- if (use.batch.stable.genes) {
    paste0("integrated.bsg.optimal.", method.name)
  } else {
    paste0("integrated.optimal.", method.name)
  }

  if (!is.null(subsample)) {
    if (run.optimal) {
      message("Running integration with selected optimal parameter...")
      best_values <- parameter.df[best_idx, , drop = FALSE]
      best_param_list <- .row_to_param_list(best_values)

      obj.raw <- rlang::inject(
        Integration(
          seurat.obj = obj.raw,
          batch = batch,
          method = method,
          new.reduction = optimal.reduction,
          ndims = ndims,
          verbose = verbose,
          default.assay = default.assay,
          !!!best_param_list
        )
      )
    }

    obj.raw@misc$optimal.reduc[[optimal.reduction]] <- best_name
    obj.raw@misc$parameter.search[[optimal.reduction]] <- search.df
    return(obj.raw)
  }

  if (!return.all) {
    for (red in reduction_names) {
      seurat.obj@reductions[[red]] <- NULL
      seurat.obj@misc$integration_effects[[red]] <- NULL
    }
  }

  seurat.obj@reductions[[optimal.reduction]] <- seurat.obj@reductions[[reduction_names[best_idx]]]
  seurat.obj@misc$optimal.reduc[[optimal.reduction]] <- best_name
  seurat.obj@misc$parameter.search[[optimal.reduction]] <- search.df
  seurat.obj
}

.row_to_param_list <- function(values) {
  lapply(values, function(x) {
    x <- x[[1]]
    if ((length(x) == 1 && is.na(x)) || identical(x, "NA")) NULL else x
  })
}

.make_param_suffix <- function(param_list) {
  paste(
    vapply(names(param_list), function(nm) {
      val <- param_list[[nm]]
      paste0(tolower(nm), ".", if (is.null(val)) "NULL" else val)
    }, character(1)),
    collapse = "."
  )
}
