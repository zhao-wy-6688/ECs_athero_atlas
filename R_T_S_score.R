##RT/S_score
library(Seurat)
library(AUCell)
set.seed(1234)
geneSets <- list(Tip = tip_genes, Stalk = stalk_genes)
DefaultAssay(cap) <- DefaultAssay(scRNA)   
slot_use <- "data"                       
expr <- GetAssayData(cap, slot = slot_use)
cells_rankings <- AUCell_buildRankings(
  expr, nCores = 1, plotStats = FALSE, verbose = FALSE
)
aucMaxRank <- ceiling(0.1 * nrow(cells_rankings))  
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings, aucMaxRank = aucMaxRank)
aucMat <- getAUC(cells_AUC)  
# 取 ST / SS
ST <- as.numeric(aucMat["Tip",   colnames(cap)])
SS <- as.numeric(aucMat["Stalk", colnames(cap)])
names(ST) <- names(SS) <- colnames(cap)
RTS <- (ST + 1) / (SS + 1)