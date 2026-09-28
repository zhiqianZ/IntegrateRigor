# Parameter search ---------------------------------------------------------

#' Search integration parameters using marginal integration scores
#' @param method.args Named list of extra arguments passed unchanged to
#'   `method` for every parameter set, e.g.
#'   `list(group.by.vars = c("sample", "tech"))` for `Harmony`.
#' @export
IntegrateRigor.ParameterS <- function(seurat.obj, parameter.df, method,
                                      batch = "batch", ndims = 30,
                                      ndims.score = min(30, ndims),
                                      ref.batch = NULL, use.prev.ref = TRUE,
                                      force.run = TRUE,
                                      K = 10, n.cores = 5, seed = 42,
                                      return.all = TRUE,
                                      use.batch.stable.genes = TRUE,
                                      verbose = FALSE,
                                      subsample = NULL,
                                      run.optimal = TRUE,
                                      default.assay = "RNA",
                                      method.args = list()) {
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
    default.assay = default.assay,
    method.args = method.args
  )
}

#' Search integration parameters using a joint multivariate GMM score
#' @param method.args Named list of extra arguments passed unchanged to
#'   `method` for every parameter set.
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
                                            default.assay = "RNA",
                                            method.args = list()) {
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
    default.assay = default.assay,
    method.args = method.args
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
                                   default.assay = "RNA",
                                   method.args = list()) {
  score_mode <- match.arg(score_mode)
  if (!is.list(method.args) || (length(method.args) > 0 &&
                                (is.null(names(method.args)) || any(names(method.args) == "")))) {
    stop("method.args must be a named list.", call. = FALSE)
  }
  clash <- intersect(names(method.args), colnames(parameter.df))
  if (length(clash) > 0) {
    stop("Arguments given in both parameter.df and method.args: ",
         paste(clash, collapse = ", "), call. = FALSE)
  }
  
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
        paste0(nm, " = ", .format_param_value(val, sep = ", ", wrap = TRUE))
      }, character(1)),
      collapse = ", "
    )
    message("Running with ", msg)
    
    param_suffix <- .make_param_suffix(param_list)
    reduction_prefix <- if (use.batch.stable.genes) "integrated.bsg." else "integrated."
    new.reduction <- paste0(reduction_prefix, method.name, ".", param_suffix)
    
    # Record identifiers up front so bookkeeping (cleanup, search.df) stays
    # consistent even if this parameter set fails below.
    nms[i] <- param_suffix
    reduction_names[i] <- new.reduction
    
    # Run integration + scoring for this parameter set. If anything fails, do
    # NOT abort the whole search: record -Inf for both metrics so this set can
    # never be selected as optimal, warn, and continue to the next set.
    metrics_i <- tryCatch({
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
            !!!method.args,
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
        bas_i <- effects[1, "BatchAlignment"]
        cis_i <- effects[1, "CellIdentity"]
      } else {
        bas_i <- stats::median(effects[, "BatchAlignment"], na.rm = TRUE)
        cis_i <- stats::median(effects[, "CellIdentity"], na.rm = TRUE)
      }
      
      # Treat a non-finite score (e.g. all dimensions returned NA) the same as
      # a hard failure, so it is never chosen as the optimum.
      if (!is.finite(bas_i) || !is.finite(cis_i)) {
        stop("integration score is not finite (all values NA/NaN/Inf)")
      }
      
      if (!is.null(subsample)) {
        obj.raw@misc$integration_effects[[new.reduction]] <- effects
      }
      
      c(bas_i, cis_i)
    }, error = function(e) {
      warning("Parameter set '", param_suffix, "' failed and will be skipped: ",
              conditionMessage(e), call. = FALSE, immediate. = TRUE)
      c(-Inf, -Inf)
    })
    
    BAS[i] <- metrics_i[1]
    CIS[i] <- metrics_i[2]
  }
  
  names(BAS) <- nms
  names(CIS) <- nms
  IS <- BAS + CIS
  names(IS) <- nms
  
  if (!any(is.finite(IS))) {
    stop("All parameter sets failed to produce a valid integration score; ",
         "no optimal parameter could be selected. See the warnings above for ",
         "the per-parameter errors.", call. = FALSE)
  }
  
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
          !!!method.args,
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
      paste0(tolower(nm), ".", .format_param_value(val, sep = "_"))
    }, character(1)),
    collapse = "."
  )
}

.format_param_value <- function(val, sep = "_", wrap = FALSE) {
  if (is.null(val)) return("NULL")
  if (length(val) <= 1) return(as.character(val))
  out <- paste(val, collapse = sep)
  if (wrap) paste0("c(", out, ")") else out
}