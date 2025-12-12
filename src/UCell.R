# ----------------------------------------------------
# UCell.R
# Methode: UCell (Rang-basierter Score)
# Ergebnis: Area Under the Curve (AUC)-Score basierend auf Rang-Statistik
# ----------------------------------------------------
library(Seurat)
library(SeuratData)
library(UCell)
library(ggplot2)

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
base_dir <- getwd()
gene_csv <- file.path(base_dir, "data", "updated_gene_list.csv")
genes <- read_gene_list(gene_csv)
sig_other <- make_signature(genes, "Platelet_Signature", "other")

print(pbmc)

# SCORING
print("--- Berechne UCell Score ---")
pbmc <- AddModuleScore_UCell(
  pbmc, 
  features = sig_other, 
  assay = "RNA"
)
score_name <- "Platelet_Signature_UCell"

# extracting high ranked genes for extend gene set method
res_ucell <- extend_gene_set(
  pbmc = pbmc,
  base_genes = genes,
  score_name = "Platelet_Signature_UCell"
)

extended_gene_set_ucell <- res_ucell$extended_genes
extended_gene_set_ucell

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
