##hdWGCNA
library(SCP)

sc_endo <- SetupForWGCNA(
  sc_endo,
  features = VariableFeatures(sc_endo),
  wgcna_name = "sc_endo_SCT"
)
sc_endo <- MetacellsByGroups(
  seurat_obj =  sc_endo,
  group.by = c("cluster_annotation","orig.ident"),
  k = 10,
  max_shared=5,
  min_cells = 300,
  #reduction = 'pca',
  reduction = 'umap',
  ident.group = 'cluster_annotation',
  slot = 'data',
  assay = 'SCT' 
)
Idents(sc_endo) = "cluster_annotation"
# set expression matrix for hdWGCNA
sc_endo <- SetDatExpr(sc_endo,assay = "SCT",slot = "data")
# Test different soft powers:
sc_endo <- TestSoftPowers(
  sc_endo,
  networkType = 'signed' # you can also use "unsigned" or "signed hybrid"
)
sc_endo <- ConstructNetwork(
  sc_endo,
  tom_name = " sc_endo_SCT_data",
  overwrite_tom = TRUE
)
# compute module eigengenes and connectivity
sc_endo <- ModuleEigengenes(sc_endo)
sc_endo <- ModuleConnectivity( sc_endo)
# rename the modules
sc_endo <- ResetModuleNames(
  sc_endo,
  new_name = "M"
)
