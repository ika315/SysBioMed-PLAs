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
T.cell.signatures <- list(
  T_Aktiv_Sig = c("CD3E", "CD8A", "IFNG", "IL2RA") 
)

print(pbmc)

# SCORING 

print("--- Berechne AddModuleScore ---")
pbmc <- AddModuleScore(
  pbmc, 
  features = T.cell.signatures, 
  name = "Seurat_Score"
)
score_name <- "Seurat_Score1" 

# Visualization: per class and UMAP
print("--- Visualisierung ---")

png(filename = "01_AddModuleScore_VlnPlot.png", width = 800, height = 600) 
VlnPlot(pbmc, features = score_name, group.by = "seurat_clusters", pt.size = 0.5, ncol = 1)
dev.off()


pbmc <- RunUMAP(pbmc, dims = 1:10)
#FeaturePlot(pbmc, features = score_name, reduction = "umap")

png(filename = "01_AddModuleScore_UMAP.png", width = 800, height = 700) 
FeaturePlot(pbmc, features = score_name, reduction = "umap")
dev.off() 
