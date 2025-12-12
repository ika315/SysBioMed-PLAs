# ----------------------------------------------------
# Vision.R
# Methode: Vision (neighborhood-based signature scoring)
# Ergebnis: AUC + local autocorrelation across the KNN graph
# ----------------------------------------------------

library(Seurat)
library(SeuratData)
library(VISION)

print("Lade und verarbeite pbmc3k (ca. 2700 Zellen)...")
data("pbmc3k.final")
pbmc <- pbmc3k.final 

DefaultAssay(pbmc) <- "RNA"
pbmc <- DietSeurat(pbmc, assays = "RNA", misc = FALSE, images = FALSE)

# preprocessing (copied from aurelias script)
pbmc <- NormalizeData(pbmc)
pbmc <- FindVariableFeatures(pbmc)
pbmc <- ScaleData(pbmc)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
pbmc <- RunUMAP(pbmc, dims = 1:10)
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.8)

# gene signature 
# need to put signature genes into different format for vision
base_dir <- getwd()
gene_csv <- file.path(base_dir, "data", "updated_gene_list.csv")
genes <- read_gene_list(gene_csv)
sig_vision <- make_signature(genes = genes, sig_name = "Platelet_Signature", method = "vision")

print(pbmc)

# ----------------------------------------------------
# SCORING WITH VISION
# ----------------------------------------------------

# convert seurat → vision object
# not this: 
# expr_mat <- as.matrix(GetAssayData(pbmc, slot = "data"))
# because vision does not want log-scaled data?

expr_mat <- as.matrix(GetAssayData(pbmc, slot = "counts"))

vision_obj <- Vision(expr_mat, signatures = list(sig_vision))

# if needed, set number of threads
# options(mc.cores = 2)

# run vision analysis (autocorrelation, AUC scoring, KNN smoothing)
vision_obj <- analyze(vision_obj)

# view results
# viewResults(vision_obj)

# extract vision AUC score
vision_auc <- getSignatureScores(vision_obj)
vision_score <- vision_auc[, "Platelet_Signature"]
vision_score <- as.numeric(vision_auc[, "Platelet_Signature"])

# store in seurat object
pbmc$Vision_Platelet_Signature <- vision_score
score_name <- "Vision_Platelet_Signature"

# extracting high ranked genes for extend gene set method
res_vision <- extend_gene_set(
  pbmc = pbmc,
  base_genes = genes,
  score_name = "Vision_Platelet_Signature"
)

extended_gene_set_vision <- res_vision$extended_genes
extended_gene_set_vision

# ----------------------------------------------------
# VISUALIZATION
# ----------------------------------------------------

# UMAP colored by Vision score
png(filename = "plots/04_Vision_UMAP.png", width = 800, height = 700)
FeaturePlot(pbmc, features = score_name, reduction = "umap")
dev.off()

# Violin plot by cluster
png(filename = "plots/04_Vision_VlnPlot.png", width = 800, height = 600)
VlnPlot(pbmc, features = score_name,
        group.by = "seurat_clusters",
        pt.size = 0.5, ncol = 1)
dev.off()