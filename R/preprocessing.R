# Preprocessing helpers ----------------------------------------------------

#' Preprocess a Seurat object for PCA-based integration methods
#' @export
Preprocess <- function(seurat.obj, batch = "batch", ndims = 30,
                       genes = NULL, ngenes = 2000,
                       default.assay = "RNA", verbose = FALSE) {
  Seurat::DefaultAssay(seurat.obj) <- default.assay
  seurat.obj <- .split_by_batch(seurat.obj, batch, default.assay)
  
  message("Preprocessing")
  seurat.obj <- Seurat::NormalizeData(seurat.obj, verbose = verbose)
  
  if (is.null(genes)) {
    if (nrow(seurat.obj) <= ngenes) {
      Seurat::VariableFeatures(seurat.obj) <- rownames(seurat.obj)
    } else {
      seurat.obj <- Seurat::FindVariableFeatures(
        seurat.obj,
        verbose = verbose,
        nfeatures = ngenes
      )
    }
  } else {
    if (length(genes) <= ngenes) {
      Seurat::VariableFeatures(seurat.obj) <- genes
    } else {
      tmp <- seurat.obj[genes, ]
      tmp <- Seurat::FindVariableFeatures(tmp, nfeatures = ngenes, verbose = verbose)
      Seurat::VariableFeatures(seurat.obj) <- Seurat::VariableFeatures(tmp)
      rm(tmp)
    }
  }
  
  seurat.obj <- Seurat::ScaleData(seurat.obj, verbose = verbose)
  message("Running PCA")
  seurat.obj <- Seurat::RunPCA(seurat.obj, npcs = ndims, verbose = verbose)
  seurat.obj
}

#' Preprocess a Seurat object for scVI integration
#' @export
Preprocess.scVI <- function(seurat.obj, batch = "batch", genes = NULL,
                            ngenes = 2000, default.assay = "RNA",
                            verbose = FALSE) {
  Seurat::DefaultAssay(seurat.obj) <- default.assay
  seurat.obj <- .split_by_batch(seurat.obj, batch, default.assay)
  
  message("Preprocessing")
  
  if (is.null(genes)) {
    if (nrow(seurat.obj) <= ngenes) {
      Seurat::VariableFeatures(seurat.obj) <- rownames(seurat.obj)
    } else {
      seurat.obj <- Seurat::NormalizeData(seurat.obj, verbose = verbose)
      seurat.obj <- Seurat::FindVariableFeatures(
        seurat.obj,
        verbose = verbose,
        nfeatures = ngenes
      )
    }
  } else {
    seurat.obj <- Seurat::NormalizeData(seurat.obj, verbose = verbose)
    if (length(genes) <= ngenes) {
      Seurat::VariableFeatures(seurat.obj) <- genes
    } else {
      tmp <- seurat.obj[genes, ]
      tmp <- Seurat::FindVariableFeatures(tmp, nfeatures = ngenes, verbose = verbose)
      Seurat::VariableFeatures(seurat.obj) <- Seurat::VariableFeatures(tmp)
      rm(tmp)
    }
  }
  
  seurat.obj
}

#' Preprocess a Seurat object for LIGER integration
#' @export
Preprocess.LIGER <- function(seurat.obj, batch = "batch", genes = NULL,
                             ngenes = 2000, default.assay = "RNA",
                             verbose = FALSE) {
  .check_package("rliger", "for LIGER preprocessing")
  
  Seurat::DefaultAssay(seurat.obj) <- default.assay
  seurat.obj <- .split_by_batch(seurat.obj, batch, default.assay)
  
  message("Preprocessing")
  
  if (is.null(genes)) {
    seurat.obj <- Seurat::FindVariableFeatures(seurat.obj, verbose = FALSE, nfeatures = ngenes)
    features <- Seurat::VariableFeatures(seurat.obj)
  } else {
    if (length(genes) <= ngenes) {
      features <- genes
    } else {
      tmp <- seurat.obj[genes, ]
      tmp <- Seurat::FindVariableFeatures(tmp, nfeatures = ngenes, verbose = FALSE)
      features <- Seurat::VariableFeatures(tmp)
      rm(tmp)
    }
  }
  
  seurat.obj <- rliger::normalize(seurat.obj, verbose = verbose)
  Seurat::VariableFeatures(seurat.obj) <- features
  seurat.obj <- rliger::scaleNotCenter(seurat.obj, verbose = verbose)
  seurat.obj
}

# Make sure the assay's layers are split by exactly the levels of `batch`.
# An existing split is kept only if every counts layer holds cells of a single
# batch level and each level has its own layer; otherwise (unsplit, or split by
# a different variable) the layers are joined and re-split by `batch`.
.split_by_batch <- function(seurat.obj, batch, assay) {
  batch_vec <- .get_batch_vector(seurat.obj, batch)
  if (anyNA(batch_vec)) {
    stop("Batch variable '", batch, "' contains NA values.", call. = FALSE)
  }
  batch_vec <- as.character(batch_vec)
  names(batch_vec) <- colnames(seurat.obj)
  
  assay_obj <- seurat.obj[[assay]]
  count_layers <- SeuratObject::Layers(assay_obj, search = "counts")
  
  if (length(count_layers) > 1) {
    per_layer <- lapply(count_layers, function(lyr) {
      unique(batch_vec[SeuratObject::Cells(assay_obj, layer = lyr)])
    })
    matches <- all(lengths(per_layer) == 1) &&
      !anyDuplicated(unlist(per_layer)) &&
      setequal(unlist(per_layer), unique(batch_vec))
    if (matches) return(seurat.obj)
    message("Layers of assay '", assay, "' are split by a different variable ",
            "than '", batch, "'; re-splitting by '", batch, "'.")
    assay_obj <- SeuratObject::JoinLayers(assay_obj)
  }
  
  seurat.obj[[assay]] <- split(assay_obj, f = batch_vec[colnames(assay_obj)])
  seurat.obj
}