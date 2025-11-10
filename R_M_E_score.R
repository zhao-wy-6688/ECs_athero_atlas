##RM/E_score
library(Seurat)
library(AUCell)
scRNA<-all
DimPlot(all)
all_gene_sets <- list(Endothelial = endothelial_genes, Mesenchymal = mesenchymal_genes)
cells_rankings <- AUCell_buildRankings(
  scRNA@assays$SCT@data,
  plotStats = FALSE
)
cells_AUC <- AUCell_calcAUC(
  all_gene_sets,
  cells_rankings,
  aucMaxRank = nrow(cells_rankings) * 0.1
)
SE_aucs <- as.numeric(getAUC(cells_AUC)["Endothelial", ])
scRNA$Endothelial_Score <- SE_aucs # 
SM_aucs <- as.numeric(getAUC(cells_AUC)["Mesenchymal", ])
scRNA$Mesenchymal_Score <- SM_aucs # 对应描述中的 SM
scRNA$RME <- (scRNA$Endothelial_Score + 1) / (scRNA$Mesenchymal_Score + 1)

