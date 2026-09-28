# Utility functions ---------------------------------------------------------

#' Make a parameter grid
#'
#' @param ... Parameter values, ideally named (e.g. `theta = c(1, 2, 4)`).
#'   Unnamed arguments are named after the expression passed in. A
#'   vector-valued parameter (e.g. Harmony's `theta` with several covariates)
#'   is given as a list of vectors, e.g. `theta = list(c(2, 1), c(4, 2))`,
#'   and becomes a list column in which each cell holds one whole vector.
#' @return A data.frame containing all parameter combinations.
#' @export
make.parameter.df <- function(...) {
  args <- list(...)
  exprs <- as.list(substitute(list(...)))[-1]
  arg_names <- names(args)
  if (is.null(arg_names)) arg_names <- rep("", length(args))
  unnamed <- arg_names == ""
  arg_names[unnamed] <- vapply(exprs[unnamed], function(e) {
    paste(deparse(e), collapse = "")
  }, character(1))
  if (anyDuplicated(arg_names)) {
    stop("Parameter names must be unique.", call. = FALSE)
  }
  
  # Expand over indices so list-valued parameters are never flattened.
  idx <- expand.grid(lapply(args, seq_along), KEEP.OUT.ATTRS = FALSE)
  out <- data.frame(row.names = seq_len(nrow(idx)))
  for (j in seq_along(args)) {
    vals <- args[[j]][idx[[j]]]
    out[[arg_names[j]]] <- if (is.list(vals)) I(unname(vals)) else vals
  }
  out
}

.check_package <- function(pkg, reason = NULL) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    msg <- paste0("Package '", pkg, "' is required")
    if (!is.null(reason)) {
      msg <- paste0(msg, " ", reason)
    }
    stop(msg, ".", call. = FALSE)
  }
  invisible(TRUE)
}

.get_batch_vector <- function(seurat.obj, batch) {
  if (!batch %in% colnames(seurat.obj@meta.data)) {
    stop("Batch variable '", batch, "' was not found in seurat.obj@meta.data.", call. = FALSE)
  }
  seurat.obj@meta.data[[batch]]
}

.log_sum_exp <- function(z) {
  m <- max(z)
  m + log(sum(exp(z - m)))
}