# Batch-stable gene estimation --------------------------------------------

#' Estimate gene-level batch stability scores
#' @export
BatchStabilityEst <- function(seurat.obj, batch = "batch",
                              ref.batch = NULL, use.prev.ref = TRUE,
                              genes = NULL, ngenes = NULL,
                              K = 5, max_iter = 10, n.cores = 10,
                              seed = 42,
                              subsample = NULL) {
  if (!is.null(subsample)) {
    stopifnot(subsample > 0, subsample <= 1)
    set.seed(seed)
    obj_raw <- seurat.obj
    id <- sample(seq_len(ncol(obj_raw)), round(subsample * ncol(obj_raw)))
    seurat.obj <- seurat.obj[, id]
  } else {
    batch_vec0 <- .get_batch_vector(seurat.obj, batch)
    if (ncol(seurat.obj) > 20000 || length(unique(batch_vec0)) > 6) {
      message("Execution may be time-consuming due to the large number of batches or cells. Consider using the subsample parameter to accelerate the process.")
    }
  }

  batch_label <- .get_batch_vector(seurat.obj, batch)

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

  if (is.null(genes) && is.null(ngenes)) {
    genes <- rownames(seurat.obj)
    message("Calculating only on genes with higher expression (expressed in > 0.5% of cells)")
    counts <- as.matrix(Seurat::GetAssayData(seurat.obj, layer = "counts"))
    print(dim(counts))
    filter_genes <- which(rowSums(counts != 0) > ncol(counts) * 0.5 / 100)
    genes <- genes[filter_genes]
    message(paste0("On ", length(genes), " genes..."))
  } else if (!is.null(ngenes)) {
    seurat.list <- Seurat::SplitObject(seurat.obj, split.by = batch)
    seurat.list <- lapply(seurat.list, function(so) {
      so <- Seurat::NormalizeData(so, verbose = FALSE)
      Seurat::FindVariableFeatures(so, verbose = FALSE)
    })
    genes <- Seurat::SelectIntegrationFeatures(seurat.list, nfeatures = ngenes)
  }

  ls <- seurat.obj$nCount_RNA
  X <- Seurat::GetAssayData(seurat.obj, layer = "counts")[genes, , drop = FALSE]

  batches <- unique(batch_label)
  query_batches <- batches[batches != ref.batch]
  effects_by_batch <- matrix(NA_real_, nrow = length(genes), ncol = length(query_batches))

  message(paste0("Computing batch stability score on ", length(query_batches), " batches...."))
  for (i in seq_along(query_batches)) {
    qb <- query_batches[i]
    message(paste0("On batch: ", qb))

    res <- pbmcapply::pbmclapply(
      seq_along(genes),
      function(g) {
        set.seed(seed)
        tryCatch({
          id_ref <- which(batch_label == ref.batch)
          id_query <- which(batch_label == qb)

          X_ref <- round(X[g, id_ref])
          X_query <- round(X[g, id_query])
          l_ref <- ls[id_ref]
          l_query <- ls[id_query]

          fit_ref <- EMmix(X_ref, l_ref, K = K, max_iter = max_iter)
          reweight_query <- Reweight(X_query, l_query, fit_ref)
          fit_query <- EMmix(X_query, l_query, K = K, max_iter = max_iter)

          -(fit_query$loglikelihood - reweight_query$loglik_reweighted) / length(l_query)
        }, error = function(e) {
          NA_real_
        })
      },
      mc.cores = n.cores
    )

    effects_by_batch[, i] <- unlist(res, use.names = FALSE)
  }

  effects <- data.frame(batch.stability.score = rowMeans(effects_by_batch, na.rm = TRUE))
  rownames(effects) <- genes

  if (!is.null(subsample)) {
    obj_raw[["RNA"]] <- Seurat::AddMetaData(obj_raw[["RNA"]], metadata = effects)
    return(obj_raw)
  }

  seurat.obj[["RNA"]] <- Seurat::AddMetaData(seurat.obj[["RNA"]], metadata = effects)
  seurat.obj
}

#' Select batch-stable genes using an elbow rule
#' @export
BatchStableGenes <- function(obj, plot = FALSE, include.na = T) {
  batch_score <- obj[["RNA"]]@meta.data[["batch.stability.score"]]
  names(batch_score) <- rownames(obj@assays$RNA)

  genes_brs <- batch_score[!is.na(batch_score)]
  csum <- sort(genes_brs, decreasing = TRUE)

  elbow_res <- elbow::elbow(data.frame(seq_along(genes_brs), csum), plot = plot)
  n <- elbow_res$seq_along.genes_brs._selected

  genes.sensitive <- names(which(csum <= csum[n]))
  if(include.na){
    genes.stable <- rownames(obj)[!rownames(obj) %in% genes.sensitive]
  }else{
    genes.stable <- names(which(csum > csum[n]))
  }

  
  obj@misc$batch.stable.genes <- genes.stable
  obj@misc$batch.unstable.genes <- genes.sensitive
  obj@misc$batch.stability.score <- batch_score
  obj
}
