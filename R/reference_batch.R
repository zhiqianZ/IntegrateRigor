# Reference batch selection ------------------------------------------------

FindReference <- function(seurat.obj, batch, min_count = 10) {
  batch_vec <- .get_batch_vector(seurat.obj, batch)

  seurat.obj <- SeuratObject::JoinLayers(seurat.obj)
  seurat.obj <- Seurat::NormalizeData(seurat.obj, verbose = FALSE)
  seurat.obj <- Seurat::FindVariableFeatures(seurat.obj, verbose = FALSE)
  seurat.obj <- Seurat::ScaleData(seurat.obj, verbose = FALSE)
  seurat.obj <- Seurat::RunPCA(seurat.obj, verbose = FALSE)
  seurat.obj <- Seurat::FindNeighbors(seurat.obj, verbose = FALSE)
  seurat.obj <- Seurat::FindClusters(seurat.obj, verbose = FALSE)

  cluster_table <- table(seurat.obj$seurat_clusters, batch_vec)
  cluster_df <- as.data.frame.matrix(cluster_table)
  coverage_count <- colSums(cluster_df > min_count)
  names(which.max(coverage_count))
}
