# ----------------------------------------------------
# AddModuleScore.R
# Methode: score_genes / AddModuleScore (Seurat)
# Ergebnis: Score ist die Differenz zwischen Signatur-Genen und Kontroll-Genen
# ----------------------------------------------------
library(Seurat)
library(SeuratData)

print("Lade und verarbeite pbmc3k (ca. 2700 Zellen)...")
data("pbmc3k.final")
pbmc <- pbmc3k.final 

DefaultAssay(pbmc) <- "RNA"
pbmc <- DietSeurat(pbmc, assays = "RNA", misc = FALSE, images = FALSE) 

# Preprocessing
pbmc <- NormalizeData(pbmc)
pbmc <- FindVariableFeatures(pbmc)
pbmc <- ScaleData(pbmc) 
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
pbmc <- RunUMAP(pbmc, dims = 1:10)
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.8)

# Signatur Definition (anpassen)
# T.cell.signatures <- list(
# T_Aktiv_Sig = c("CD3E", "CD8A", "IFNG", "IL2RA") 
# )

sig_other <- make_signature(genes, "Platelet_Signature", "other")

print(pbmc)

# SCORING 

print("--- Berechne AddModuleScore ---")
pbmc <- AddModuleScore(
  pbmc, 
  features = sig_other, 
  name = "Seurat_Score"
)
score_name <- "Seurat_Score1" 

# extracting high ranked genes for extend gene set method
high_cells_Seurat <- names(pbmc[[score_name, drop = TRUE]])[
  pbmc[[score_name, drop = TRUE]] >
    quantile(pbmc[[score_name, drop = TRUE]], 0.9)
]

pbmc$SeuratScore_group <- ifelse(
  colnames(pbmc) %in% high_cells_Seurat,
  "high",
  "low"
)

markers_Seurat <- FindMarkers(
  pbmc,
  ident.1 = "high",
  ident.2 = "low",
  group.by = "SeuratScore_group",
  logfc.threshold = 0,  
  min.pct = 0.05,
  test.use = "wilcox"
)

top_genes_Seurat <- rownames(
  markers_Seurat[
    order(markers_Seurat$avg_log2FC, decreasing = TRUE),
  ]
)[1:50]

extended_gene_set <- extend_gene_set(genes, top_genes_Seurat)
# print(extended_gene_set)

# Visualization: per class and UMAP
print("--- Visualisierung ---")

png(filename = "plots/01_AddModuleScore_VlnPlot.png", width = 800, height = 600) 
VlnPlot(pbmc, features = score_name, group.by = "seurat_clusters", pt.size = 0.5, ncol = 1)
dev.off()


pbmc <- RunUMAP(pbmc, dims = 1:10)
#FeaturePlot(pbmc, features = score_name, reduction = "umap")

png(filename = "plots/01_AddModuleScore_UMAP.png", width = 800, height = 700) 
FeaturePlot(pbmc, features = score_name, reduction = "umap")
dev.off() 
