# ----------------------------------------------------
# UCell.R
# Methode: UCell (Rang-basierter Score)
# Ergebnis: Area Under the Curve (AUC)-Score basierend auf Rang-Statistik
# ----------------------------------------------------
library(Seurat)
library(SeuratData)
library(UCell)
library(ggplot2)

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
print("--- Berechne UCell Score ---")
pbmc <- AddModuleScore_UCell(
  pbmc, 
  features = T.cell.signatures, 
  assay = "RNA"
)
score_name <- "T_Aktiv_Sig_UCell"

# Visualization
print("--- Visualisierung ---")

png(filename = "plots/02_UCell_DensityPlot.png", width = 800, height = 600) 
plot_data <- data.frame(Score = pbmc[[score_name, drop = TRUE]])

ggplot(plot_data, aes(x = Score)) + 
  geom_density(fill = "#00BFC4", alpha = 0.7) +
  labs(title = "Density Plot of UCell Score") +
  theme_classic()
dev.off()

png(filename = "plots/02_UCell_VlnPlot.png", width = 800, height = 600) 
VlnPlot(pbmc, features = score_name, group.by = "seurat_clusters", pt.size = 0.5, ncol = 1)
dev.off()
