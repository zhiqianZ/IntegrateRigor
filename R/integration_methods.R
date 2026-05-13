# Integration method wrappers ---------------------------------------------

#' Run an integration backend on a Seurat object
#' @export
Integration <- function(seurat.obj, batch = "batch", method,
                        new.reduction = NULL, ndims = 30,
                        verbose = FALSE, default.assay = "RNA", ...) {
  Seurat::DefaultAssay(seurat.obj) <- default.assay

  if (is.null(new.reduction)) {
    method_name <- deparse(substitute(method))
    new.reduction <- paste0("integrated.", tolower(method_name))
  }

  method(
    seurat.obj = seurat.obj,
    new.reduction = new.reduction,
    ndims = ndims,
    default.assay = default.assay,
    verbose = verbose,
    ...
  )
}

#' CCA integration wrapper
#' @export
CCA <- function(seurat.obj, new.reduction, ndims, default.assay, verbose = FALSE, ...) {
  Seurat::IntegrateLayers(
    object = seurat.obj,
    method = Seurat::CCAIntegration,
    orig.reduction = "pca",
    new.reduction = new.reduction,
    features = Seurat::VariableFeatures(object = seurat.obj, assay = default.assay),
    verbose = verbose,
    ...
  )
}

#' RPCA integration wrapper
#' @export
RPCA <- function(seurat.obj, new.reduction, ndims, default.assay, verbose = FALSE, ...) {
  Seurat::IntegrateLayers(
    object = seurat.obj,
    method = Seurat::RPCAIntegration,
    orig.reduction = "pca",
    new.reduction = new.reduction,
    features = Seurat::VariableFeatures(object = seurat.obj, assay = default.assay),
    verbose = verbose,
    ...
  )
}

#' Harmony integration wrapper
#' @export
Harmony <- function(seurat.obj, new.reduction, ndims, default.assay, verbose = FALSE, ...) {
  Seurat::IntegrateLayers(
    object = seurat.obj,
    method = Seurat::HarmonyIntegration,
    orig.reduction = "pca",
    new.reduction = new.reduction,
    features = Seurat::VariableFeatures(object = seurat.obj, assay = default.assay),
    verbose = verbose,
    ...
  )
}

#' scVI integration wrapper
#' @export
scVI <- function(seurat.obj, new.reduction, ndims, default.assay, verbose = FALSE, ...) {
  .check_package("reticulate", "for scVI integration")
  .check_package("SeuratWrappers", "for scVI batch extraction helpers")

  features <- Seurat::VariableFeatures(object = seurat.obj, assay = default.assay)

  Seurat::IntegrateLayers(
    object = seurat.obj,
    method = scVIIntegration_custom,
    new.reduction = new.reduction,
    verbose = verbose,
    conda_env = "scvi",
    features = features,
    layers = "counts",
    orig.reduction = NULL,
    scale.layer = NULL,
    assay = default.assay,
    ndimss = ndims,
    ...
  )
}

#' FastMNN integration wrapper
#' @export
FastMNN <- function(seurat.obj, new.reduction, ndims, default.assay, verbose = FALSE, ...) {
  .check_package("SeuratWrappers", "for FastMNN integration")
  fastmnn_fun <- getFromNamespace("FastMNNIntegration", "SeuratWrappers")

  Seurat::IntegrateLayers(
    object = seurat.obj,
    method = fastmnn_fun,
    orig.reduction = "pca",
    new.reduction = new.reduction,
    d = ndims,
    features = Seurat::VariableFeatures(object = seurat.obj, assay = default.assay),
    verbose = verbose,
    ...
  )
}

#' LIGER integration wrapper
#' @export
LIGER <- function(seurat.obj, new.reduction, ndims, default.assay, verbose = FALSE, ...) {
  .check_package("rliger", "for LIGER integration")

  seurat.obj <- rliger::runIntegration(seurat.obj, k = ndims, verbose = verbose, ...)
  seurat.obj <- rliger::alignFactors(seurat.obj, method = "centroidAlign")
  seurat.obj[[new.reduction]] <- seurat.obj[["inmfNorm"]]
  seurat.obj
}

scVIIntegration_custom <- function(
    object,
    features = NULL,
    layers = "counts",
    conda_env = NULL,
    new.reduction = "integrated.dr",
    ndimss = 30,
    nlayers = 2,
    nhidden = 128,
    dropout.rate = 0.1,
    gene_likelihood = "nb",
    max_epochs = NULL,
    ...) {
  .check_package("reticulate", "for scVI integration")
  .check_package("SeuratWrappers", "for scVI batch extraction helpers")

  reticulate::use_condaenv(conda_env, required = TRUE)
  sc <- reticulate::import("scanpy", convert = FALSE)
  scvi <- reticulate::import("scvi", convert = FALSE)
  scipy <- reticulate::import("scipy", convert = FALSE)

  if (is.null(max_epochs)) {
    max_epochs <- reticulate::r_to_py(max_epochs)
  } else {
    max_epochs <- as.integer(max_epochs)
  }

  if (inherits(object, what = "SCTAssay")) {
    find_sct_batches <- getFromNamespace(".FindSCTBatches", "SeuratWrappers")
    batches <- find_sct_batches(object)
  } else {
    find_batches <- getFromNamespace(".FindBatches", "SeuratWrappers")
    batches <- find_batches(object, layers = layers)
    object <- SeuratObject::JoinLayers(object = object, layers = "counts")
  }

  adata <- sc$AnnData(
    X = scipy$sparse$csr_matrix(
      Matrix::t(SeuratObject::LayerData(object, layer = "counts")[features, ])
    ),
    obs = batches,
    var = object[[]][features, ]
  )

  scvi$model$SCVI$setup_anndata(adata, batch_key = "batch")

  model <- scvi$model$SCVI(
    adata = adata,
    n_latent = as.integer(ndimss),
    n_layers = as.integer(nlayers),
    n_hidden = as.integer(nhidden),
    dropout_rate = dropout.rate,
    gene_likelihood = gene_likelihood
  )
  model$train(max_epochs = max_epochs)

  latent <- model$get_latent_representation()
  latent <- as.matrix(latent)
  rownames(latent) <- reticulate::py_to_r(adata$obs$index$values)
  colnames(latent) <- paste0(new.reduction, "_", seq_len(ncol(latent)))

  suppressWarnings(
    latent.dr <- SeuratObject::CreateDimReducObject(
      embeddings = latent,
      key = new.reduction
    )
  )

  output.list <- list(latent.dr)
  names(output.list) <- new.reduction
  output.list
}
