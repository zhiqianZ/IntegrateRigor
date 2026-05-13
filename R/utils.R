# Utility functions ---------------------------------------------------------

#' Make a parameter grid
#'
#' @param ... Named vectors of parameter values.
#' @return A data.frame containing all parameter combinations.
#' @export
make.parameter.df <- function(...) {
  args <- list(...)
  arg_names <- vapply(substitute(list(...))[-1], deparse, character(1))
  names(args) <- arg_names

  expand.grid(
    args,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
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
