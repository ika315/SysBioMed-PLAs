# ----------------------------------------------------
# 03_AUCell.R
# Methode: AUCell (Ursprüngliche Rang-Methode)
# Ergebnis: Area Under the Curve (AUC)-Score basierend auf Gen-Rankings
# ----------------------------------------------------
library(Seurat)
library(SeuratData)
library(AUCell)

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
print("--- Berechne AUCell Score ---")
expression_matrix <- GetAssayData(pbmc, slot = "data") 

cells_rankings <- AUCell_buildRankings(expression_matrix, plotStats=FALSE)
cells_AUC <- AUCell_calcAUC(T.cell.signatures, cells_rankings)

#
pbmc$AUCell_T_Aktiv <- as.numeric(getAUC(cells_AUC)[1, ])
score_name <- "AUCell_T_Aktiv" 
# Visulaization
print("--- Visualisierung ---")

# a) UMAP des Scores
pbmc <- RunUMAP(pbmc, dims = 1:10)
png(filename = "03_AUCell_UMAP.png", width = 800, height = 700) 
FeaturePlot(pbmc, features = score_name, reduction = "umap")
dev.off()

# b) Per-Cluster Verteilung
png(filename = "03_AUCell_VlnPlot.png", width = 800, height = 600) 
VlnPlot(pbmc, features = score_name, group.by = "seurat_clusters", pt.size = 0.5, ncol = 1)
dev.off()