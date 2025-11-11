Monocle3Simple <- function(
    srt,
    assay = NULL,
    slot  = "counts",
    reduction = Seurat::DefaultReduction(srt),  
    reuse_umap = TRUE,      
    num_dim = 50,           
    use_partition = TRUE,   
    root_cells = NULL,      
    seed = 11,
    plot = FALSE            
){
  set.seed(seed)
  stopifnot(requireNamespace("Seurat", quietly = TRUE),
            requireNamespace("monocle3", quietly = TRUE),
            requireNamespace("SingleCellExperiment", quietly = TRUE),
            requireNamespace("Matrix", quietly = TRUE),
            packageVersion("monocle3") >= "1.2.0")
  
  assay <- assay %||% Seurat::DefaultAssay(srt)
  expr <- Seurat::GetAssayData(srt, assay = assay, slot = slot)
  if (!inherits(expr, "dgCMatrix")) expr <- as(expr, "dgCMatrix")
  cell_md <- srt@meta.data[ colnames(expr), , drop = FALSE ]
  gene_md <- data.frame(gene_short_name = rownames(expr), row.names = rownames(expr))
  cds <- monocle3::new_cell_data_set(expr, cell_metadata = cell_md, gene_metadata = gene_md)
  sf_col <- paste0("nCount_", assay)
  if (!"Size_Factor" %in% colnames(SummarizedExperiment::colData(cds)) &&
      sf_col %in% colnames(srt@meta.data)) {
    cds$Size_Factor <- srt@meta.data[[sf_col]]
  }
  ok_to_reuse <- isTRUE(reuse_umap) && reduction %in% Seurat::Reductions(srt)
  if (ok_to_reuse) {
    umap_seu <- Seurat::Embeddings(srt, reduction)[colnames(cds), , drop = FALSE]
    SingleCellExperiment::reducedDims(cds)$UMAP <- umap_seu
  } else {
    cds <- monocle3::preprocess_cds(cds, num_dim = num_dim, method = "PCA")
    cds <- monocle3::reduce_dimension(cds, reduction_method = "UMAP", preprocess_method = "PCA")
  }
  
  cds <- monocle3::cluster_cells(cds, reduction_method = "UMAP")
  cds <- monocle3::learn_graph(cds, use_partition = use_partition, close_loop = TRUE)
  
  if (is.null(root_cells)) {
    cl <- SummarizedExperiment::colData(cds)$cluster
    major <- names(sort(table(cl), decreasing = TRUE))[1]
    root_cells <- names(which(cl == major))[1]
  }
  cds <- monocle3::order_cells(cds, root_cells = root_cells)
  srt[["Monocle3_Pseudotime"]] <- monocle3::pseudotime(cds)
  srt[["Monocle3_clusters"]]   <- SummarizedExperiment::colData(cds)$cluster
  srt@tools$Monocle3 <- list(cds = cds)
  
  if (isTRUE(plot)) {
    print(monocle3::plot_cells(cds, color_cells_by = "pseudotime",
                               show_trajectory_graph = TRUE,
                               label_groups_by_cluster = FALSE,
                               label_branch_points = FALSE,
                               label_leaves = FALSE))
  }
  return(srt)
}

