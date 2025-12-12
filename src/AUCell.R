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
# T.cell.signatures <- list(
  # T_Aktiv_Sig = c("CD3E", "CD8A", "IFNG", "IL2RA") 
# )

sig_other <- make_signature(genes, "Platelet_Signature", "other")

print(pbmc)

# SCORING 
print("--- Berechne AUCell Score ---")
expression_matrix <- GetAssayData(pbmc, slot = "data") 

cells_rankings <- AUCell_buildRankings(expression_matrix, plotStats=FALSE)
cells_AUC <- AUCell_calcAUC(sig_other, cells_rankings)

#
pbmc$AUCell_Platelet_Signature <- as.numeric(getAUC(cells_AUC)[1, ])
score_name <- "AUCell_Platelet_Signature" 

# extracting high ranked genes for extend gene set method
high_cells_AUC <- names(pbmc$AUCell_Platelet_Signature)[pbmc$AUCell_Platelet_Signature > quantile(pbmc$AUCell_Platelet_Signature, 0.9)]

# matrix of gene expression ranks where each col=cell, each row=gene, cell=how highly is gene expressed in cell
rank_mat <- assay(cells_rankings)

# collapse per-cell rankings to one global ranking (each gene has one number)
# across platelet-enriched cells, where does this gene usually appear in the expression ranking?
mean_rank <- rowMeans(rank_mat[, high_cells_AUC, drop=F]) # the smaller the rank the higher its position so higher expression

top_genes_AUC <- names(sort(mean_rank))[1:200]
extended_gene_set <- extend_gene_set(genes, top_genes_AUC)
#print(extended_gene_set)

# Visulaization
print("--- Visualisierung ---")

# a) UMAP des Scores
pbmc <- RunUMAP(pbmc, dims = 1:10)
png(filename = "plots/03_AUCell_UMAP.png", width = 800, height = 700) 
FeaturePlot(pbmc, features = score_name, reduction = "umap")
dev.off()

# b) Per-Cluster Verteilung
png(filename = "plots/03_AUCell_VlnPlot.png", width = 800, height = 600) 
VlnPlot(pbmc, features = score_name, group.by = "seurat_clusters", pt.size = 0.5, ncol = 1)
dev.off()