# ----------------------------------------------------
# AddModuleScore.R
# Methode: score_genes / AddModuleScore (Seurat)
# Ergebnis: Score ist die Differenz zwischen Signatur-Genen und Kontroll-Genen
# ----------------------------------------------------
library(Seurat)
<<<<<<< HEAD
=======
#library(SeuratData)
>>>>>>> abec8e1 (Updated AddModuleScore.R and plots)

ds_name <- "seu_sx_final"
path = "~/SysBioMed-PLAs/data/seu_sx_final.rds"

if (!file.exists(path)) {
    stop(paste("Fehler: Datensatz nicht gefunden unter", path))
}

print(paste("Lade prozessiertes Seurat-Objekt:", path))
pbmc <- readRDS(path)
DefaultAssay(pbmc) <- "RNA"

# Preprocessing
#pbmc <- NormalizeData(pbmc, verbose = FALSE)
#pbmc <- FindVariableFeatures(pbmc, verbose = FALSE)
#pbmc <- ScaleData(pbmc, verbose = FALSE) 
#pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc), verbose = FALSE)
#pbmc <- RunUMAP(pbmc, dims = 1:10)
#pbmc <- FindNeighbors(pbmc, dims = 1:10, verbose = FALSE)
#pbmc <- FindClusters(pbmc, resolution = 0.8)

# Signatur Definition (anpassen)
# T.cell.signatures <- list(
# T_Aktiv_Sig = c("CD3E", "CD8A", "IFNG", "IL2RA") 
# )

base_dir <- getwd()
source(file.path(base_dir, "src","read_and_extend_gene_list.R"))
gene_csv <- file.path(base_dir, "data", "updated_gene_list.csv")
genes <- read_gene_list(gene_csv)
sig_other <- make_signature(genes, "Platelet_Signature", "other")

print(pbmc)

# SCORING 

print("--- Berechne AddModuleScore ---")
pbmc <- AddModuleScore(
  pbmc, 
  features = sig_other, 
  name = "AddModule_Score"
)
score_name <- "AddModule_Score1" 

summary(pbmc[[score_name]])

# extracting high ranked genes for extend gene set method
#res_seurat <- extend_gene_set(
#  pbmc = pbmc,
#  base_genes = genes,
#  score_name = "AddModule_Score1"
#)

#extended_gene_set_seurat <- res_seurat$extended_genes
#extended_gene_set_seurat

print("Führe Z-Score Normalisierung durch...")

pbmc$AddModuleScore_ZScore <- scale(pbmc[[score_name]][, 1])[, 1]

summary(pbmc$AddModuleScore_ZScore)


# Visualization: per class and UMAP
print("--- Visualisierung ---")

# RAW Score Violin Plot
png(filename = paste0(base_dir, "/plots/04_AddModuleScore_VlnPlot_", ds_name, "_RAW.png"),width = 800, height = 600)
VlnPlot(pbmc, features = score_name,group.by = "seurat_clusters", pt.size = 0.1, ncol = 1)
dev.off()

# RAW Score UMAP
png(
  filename = paste0(base_dir, "/plots/04_AddModuleScore_UMAP_", ds_name, "_RAW.png"),
  width = 800,
  height = 700
)
FeaturePlot(pbmc, features = score_name, reduction = "umap")
dev.off()

# Z-Score UMAP
png(
  filename = paste0(base_dir, "/plots/04_AddModuleScore_UMAP_", ds_name, "_ZScore.png"),
  width = 800,
  height = 700
)
FeaturePlot(pbmc, features = "AddModuleScore_ZScore", reduction = "umap")
dev.off()

# Violin Z-Score
png(
  filename = paste0(base_dir, "/plots/04_AddModuleScore_VlnPlot_", ds_name, "_ZScore.png"),
  width = 800,
  height = 600
)
VlnPlot(
  pbmc,
  features = "AddModuleScore_ZScore",
  group.by = "seurat_clusters",
  pt.size = 0.1
)
dev.off()
