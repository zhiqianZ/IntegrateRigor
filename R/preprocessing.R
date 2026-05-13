# Preprocessing helpers ----------------------------------------------------

#' Preprocess a Seurat object for PCA-based integration methods
#' @export
Preprocess <- function(seurat.obj, batch = "batch", ndims = 30,
                       genes = NULL, ngenes = 2000,
                       default.assay = "RNA", verbose = FALSE) {
  Seurat::DefaultAssay(seurat.obj) <- default.assay
  batch_vec <- .get_batch_vector(seurat.obj, batch)

  layer_names <- SeuratObject::Layers(seurat.obj[[default.assay]])
  if (length(grep("^data\\.|^counts\\.", layer_names)) == 0) {
    seurat.obj[[default.assay]] <- split(seurat.obj[[default.assay]], f = batch_vec)
  }

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
  batch_vec <- .get_batch_vector(seurat.obj, batch)
  seurat.obj[[default.assay]] <- split(seurat.obj[[default.assay]], f = batch_vec)

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
  batch_vec <- .get_batch_vector(seurat.obj, batch)
  seurat.obj[[default.assay]] <- split(seurat.obj[[default.assay]], f = batch_vec)

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
