# ----------------------------------------------------
# AddModuleScore.R
# Methode: score_genes / AddModuleScore (Seurat)
# Ergebnis: Score ist die Differenz zwischen Signatur-Genen und Kontroll-Genen
# ----------------------------------------------------
library(Seurat)

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

base_dir <- getwd()
gene_csv <- file.path(base_dir, "data", "updated_gene_list.csv")
genes <- read_gene_list(gene_csv)
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
res_seurat <- extend_gene_set(
  pbmc = pbmc,
  base_genes = genes,
  score_name = "Seurat_Score1"
)

extended_gene_set_seurat <- res_seurat$extended_genes
extended_gene_set_seurat

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
