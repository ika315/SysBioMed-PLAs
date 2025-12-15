# ----------------------------------------------------
# 03_AUCell.R
# Methode: AUCell (Ursprüngliche Rang-Methode)
# Ergebnis: Area Under the Curve (AUC)-Score basierend auf Gen-Rankings
# ----------------------------------------------------
library(Seurat)
library(AUCell)

ds_name <- "seu_sx_final"
path = "~/SysBioMed-PLAs/data/seu_sx_final.rds"

if (!file.exists(path)) {
    stop(paste("Fehler: Datensatz nicht gefunden unter", path))
}

print(paste("Lade prozessiertes Seurat-Objekt:", path))
pbmc <- readRDS(path)

DefaultAssay(pbmc) <- "RNA"
pbmc <- DietSeurat(
  pbmc,
  assays = "RNA",
  reductions = c("pca", "umap"),
  graphs = c("RNA_nn", "RNA_snn"),
  misc = FALSE,
  images = FALSE
)

# Preprocessing
pbmc <- NormalizeData(pbmc, verbose = FALSE)
pbmc <- FindVariableFeatures(pbmc, verbose = FALSE)
pbmc <- ScaleData(pbmc, verbose = FALSE) 
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc), verbose = FALSE)
#pbmc <- RunUMAP(pbmc, dims = 1:10, verbose = FALSE)
pbmc <- FindNeighbors(pbmc, dims = 1:10, verbose = FALSE)
pbmc <- RunUMAP(pbmc, dims = 1:10, verbose = FALSE)
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
print("--- Berechne AUCell Score ---")

#features_to_use <- VariableFeatures(pbmc)
expression_matrix <- GetAssayData(pbmc, layer = "data")

cells_rankings <- AUCell_buildRankings(expression_matrix, plotStats=FALSE)
cells_AUC <- AUCell_calcAUC(sig_other, cells_rankings)

pbmc$AUCell_Platelet_Signature <- as.numeric(getAUC(cells_AUC)[1, ])
score_name <- "AUCell_Platelet_Signature" 

str(pbmc$AUCell_Platelet_Signature)
summary(pbmc$AUCell_Platelet_Signature)

#print("--- Extracting high ranked genes ---")
# extracting high ranked genes for extend gene set method
#res_auc <- extend_gene_set(
#  pbmc = pbmc,
#  base_genes = genes,
#  score_name = "AUCell_Platelet_Signature"
#)

#extended_gene_set_auc <- res_auc$extended_genes
#extended_gene_set_auc

# intuitive version
# extracting high ranked genes for extend gene set method
# high_cells_AUC <- names(pbmc$AUCell_Platelet_Signature)[pbmc$AUCell_Platelet_Signature > quantile(pbmc$AUCell_Platelet_Signature, 0.9)]
# matrix of gene expression ranks where each col=cell, each row=gene, cell=how highly is gene expressed in cell
#rank_mat <- assay(cells_rankings)

# collapse per-cell rankings to one global ranking (each gene has one number)
# across platelet-enriched cells, where does this gene usually appear in the expression ranking?
#mean_rank <- rowMeans(rank_mat[, high_cells_AUC, drop=F]) # the smaller the rank the higher its position so higher expression

#top_genes_AUC <- names(sort(mean_rank))[1:50]
#extended_gene_set <- extend_gene_set(genes, top_genes_AUC)
#print(extended_gene_set)

# Visulaization
print("--- Visualisierung ---")

# a) UMAP des Scores
#pbmc <- RunUMAP(pbmc, dims = 1:10) 
plot_filename <- paste0(base_dir, "plots","03_AUCell_UMAP_", ds_name, "_RAW_Score.png")
png(filename = plot_filename, width = 800, height = 700)
FeaturePlot(pbmc, features = score_name, reduction = "umap") 
dev.off()

# b) Per-Cluster Verteilung
plot_filename <- paste0(base_dir, "plots","03_AUCell_VlnPlot_", ds_name, "_RAW_Score.png")
png(filename = plot_filename, width = 800, height = 600)
VlnPlot(pbmc, features = score_name, group.by = "seurat_clusters", pt.size = 0.5, ncol = 1)
dev.off()

# Z-Score
print("Führe Z-Score Normalisierung durch...")

pbmc$AUCell_ZScore <- scale(pbmc[[score_name]][, 1])[, 1]
score_plot_name <- "AUCell_ZScore"

summary(pbmc$AUCell_ZScore)

# UMAP Plot 
#pbmc <- RunUMAP(pbmc, dims = 1:10)
plot_filename_umap <- paste0(base_dir, "plots","03_AUCell_UMAP_", ds_name, "_ZScore.png")
png(filename = plot_filename_umap, width = 1000, height = 800) 
print(FeaturePlot(pbmc, features = score_plot_name, reduction = "umap", pt.size = 0.5))
dev.off()
print(paste("UMAP Z-Score Plot gespeichert:", plot_filename_umap))


# Violin Plot
plot_filename_vln <- paste0(base_dir, "plots", "03_AUCell_VlnPlot_", ds_name, "_ZScore.png")
png(filename = plot_filename_vln, width = 800, height = 600) 
print(VlnPlot(pbmc, features = score_plot_name, group.by = "seurat_clusters", pt.size = 0.1, ncol = 1))
dev.off()
print(paste("Violin Z-Score Plot gespeichert:", plot_filename_vln))
