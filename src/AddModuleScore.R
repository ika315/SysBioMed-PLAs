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

# platelet version
T.cell.signatures <- list(
  T_Aktiv_Sig = c(
    "ACTB", "ACT", "ACTA", "FIBB", "VASP", "ITA2B", "ITB3", "GPV", "GP1BA", "GP1BB",
    "GPIX", "GELS", "PECA1", "LYAM3", "KPCB", "GPVI", "ITA2", "ICAM2", "CD63", "TSN9",
    "PRIO", "FCG2A", "LAMB1", "CD92", "CD41", "CD61", "CD62P", "CD9", "CD23", "CD31",
    "CD36", "CD42a", "CD42b", "CD42c", "CD42d", "CD49b", "CD49f", "CD51", "CD84", "CD109",
    "CD110", "CD147", "CD151", "CD226", "CD107a", "CD107b", "ITGA2B", "ITGB3", "GP9", "GP5",
    "SELP", "VWF", "PF4", "PPBP", "RAB27B", "VAMP8", "SNAP23", "CLIC3", "THBS1", "FERMT3",
    "TREML1", "TUBB1", "ANXA3", "ANXA5", "GNAS", "SPARC", "F13A1"
  ))

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

png(filename = "plots/01_AddModuleScore_VlnPlot.png", width = 800, height = 600) 
VlnPlot(pbmc, features = score_name, group.by = "seurat_clusters", pt.size = 0.5, ncol = 1)
dev.off()


pbmc <- RunUMAP(pbmc, dims = 1:10)
#FeaturePlot(pbmc, features = score_name, reduction = "umap")

png(filename = "plots/01_AddModuleScore_UMAP.png", width = 800, height = 700) 
FeaturePlot(pbmc, features = score_name, reduction = "umap")
dev.off() 
