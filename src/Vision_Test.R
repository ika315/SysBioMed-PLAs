# ----------------------------------------------------
# Vision.R
# Methode: Vision (neighborhood-based signature scoring)
# Ergebnis: AUC + local autocorrelation across the KNN graph
# ----------------------------------------------------

library(VISION)
library(Seurat)
library(ggplot2)
library(SeuratObject)

#options(mc.cores = as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", 1)))
options(mc.cores = 1)

base_dir <- getwd()

ds_name <- "seu_sx_final"
path = "~/SysBioMed-PLAs/data/seu_sx_final.rds"

if (!file.exists(path)) {
  stop(paste("Fehler: Datensatz nicht gefunden unter", path))
}

print(paste("Lade prozessiertes Seurat-Objekt:", path))
pbmc <- readRDS(path)

#DefaultAssay(pbmc) <- "RNA"
#pbmc <- DietSeurat(
#  pbmc,
#  assays = "RNA",
#  reductions = c("pca", "umap"),
#  graphs = c("RNA_nn", "RNA_snn"),
#  misc = FALSE,
#  images = FALSE
#)

# Preprocessing
#pbmc <- NormalizeData(pbmc, verbose = FALSE)
#pbmc <- FindVariableFeatures(pbmc, verbose = FALSE)
#pbmc <- ScaleData(pbmc, verbose = FALSE) 
#pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc), verbose = FALSE)
#pbmc <- RunUMAP(pbmc, dims = 1:10, verbose = FALSE)
#pbmc <- FindNeighbors(pbmc, dims = 1:10, verbose = FALSE)
#pbmc <- RunUMAP(pbmc, dims = 1:10, verbose = FALSE)
#pbmc <- FindClusters(pbmc, resolution = 0.8)

# Signatur Definition (anpassen)
# T.cell.signatures <- list(
# T_Aktiv_Sig = c("CD3E", "CD8A", "IFNG", "IL2RA") 
# )
# gene signature 
# need to put signature genes into different format for vision
path_gene_list = "~/SysBioMed-PLAs/data/naive_b_cells.csv"
source("~/SysBioMed-PLAs/src/read_and_extend_gene_list.R")
genes <- read_gene_list(path_gene_list)
sig_vision <- make_signature(genes = genes, sig_name = "B_Cell_Signature", method = "vision")

print(pbmc)

# ----------------------------------------------------
# SCORING WITH VISION
# ----------------------------------------------------

# convert seurat → vision object
#counts_mat <- GetAssayData(pbmc, assay = "RNA", layer = "counts")
#pbmc[["RNA_vision"]] <- CreateAssayObject(counts = counts_mat)
#DefaultAssay(pbmc) <- "RNA_vision"

pbmc$celltype.l2 <- factor(pbmc$celltype.l2)

pbmc[["RNA_v3"]] <- as(pbmc[["RNA"]], "Assay")
DefaultAssay(pbmc) <- "RNA_v3"

vision_obj <- Vision(
  pbmc,
  assay = "RNA_v3",
  signatures = list(sig_vision),
  pool = TRUE
)

# 4) run vision analysis (autocorrelation, AUC scoring, KNN smoothing)
vision_obj <- analyze(vision_obj)

# view results
# viewResults(vision_obj)

# extract vision AUC score
vision_auc <- getSignatureScores(vision_obj)
vision_score <- vision_auc[, "B_Cell_Signature"]
vision_score <- as.numeric(vision_auc[, "B_Cell_Signature"])

# store in seurat object
pbmc$Vision_B_Cell_Signature <- vision_score
score_name <- "Vision_B_Cell_Signature"

# extracting high ranked genes for extend gene set method
#res_vision <- extend_gene_set(
#  pbmc = pbmc,
#  base_genes = genes,
#  score_name = "Vision_B_Cell_Signature"
#)

#extended_gene_set_vision <- res_vision$extended_genes
#extended_gene_set_vision

# ground truth vs pred

## ----------------------------------------------------
## Z-SCORE + CONFUSION MATRIX (ALIGNED LOGIC)
## ----------------------------------------------------

# Ground truth
pbmc$GT_Response <- ifelse(pbmc$celltype.l2 == "B naive", 1, 0)

# Z-score (same as benchmarking script)
pbmc$Vision_ZScore <- as.vector(scale(pbmc$Vision_B_Cell_Signature))

THRESHOLD_Z <- 1.5

# Prediction
pbmc$Prediction <- ifelse(pbmc$Vision_ZScore > THRESHOLD_Z, "Positive", "Negative")
pbmc$GT_Class   <- ifelse(pbmc$GT_Response == 1, "Positive", "Negative")

# Confusion matrix
ConfusionMatrix <- table(
  Predicted = pbmc$Prediction,
  Actual    = pbmc$GT_Class
)

TP <- ConfusionMatrix["Positive", "Positive"]
FP <- ConfusionMatrix["Positive", "Negative"]
FN <- ConfusionMatrix["Negative", "Positive"]
TN <- ConfusionMatrix["Negative", "Negative"]

# Metrics
Precision   <- TP / (TP + FP)
Recall      <- TP / (TP + FN)
Specificity <- TN / (TN + FP)
Accuracy    <- (TP + TN) / sum(ConfusionMatrix)
F1_Score    <- 2 * (Precision * Recall) / (Precision + Recall)

# Error type (same naming as new script)
pbmc$confusion <- dplyr::case_when(
  pbmc$Prediction == "Positive" & pbmc$GT_Class == "Positive" ~ "TP",
  pbmc$Prediction == "Positive" & pbmc$GT_Class == "Negative" ~ "FP",
  pbmc$Prediction == "Negative" & pbmc$GT_Class == "Positive" ~ "FN",
  pbmc$Prediction == "Negative" & pbmc$GT_Class == "Negative" ~ "TN"
)

## ----------------------------------------------------
## PRINT METRICS
## ----------------------------------------------------

message("=== Vision Benchmarking Metrics (Naive B) ===")

print(ConfusionMatrix)

metrics_df <- data.frame(
  Metric = c("Precision", "Recall", "Specificity", "Accuracy", "F1"),
  Value  = c(Precision, Recall, Specificity, Accuracy, F1_Score)
)

print(metrics_df)

message(
  sprintf(
    "Precision: %.3f | Recall: %.3f | F1: %.3f | Accuracy: %.3f",
    Precision, Recall, F1_Score, Accuracy
  )
)

## ----------------------------------------------------
## SAVE FINAL PBMC OBJECT
## ----------------------------------------------------

dir.create("results", recursive = TRUE, showWarnings = FALSE)

pbmc_out <- file.path(
  "results",
  paste0("pbmc_", ds_name, "_Vision_NaiveB.rds")
)

saveRDS(pbmc, pbmc_out)

message("Saved final pbmc object to:")
message(pbmc_out)


tryCatch({
  
  message("=== Starting plotting block ===")
  
  ## ------------------------------------------------
  ## Ensure output directories exist
  ## ------------------------------------------------
  plot_dir <- file.path(base_dir, "plots", "seurat_downgrade")
  dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
  
  ## ------------------------------------------------
  ## RAW SCORE PLOTS (DESCRIPTIVE)
  ## ------------------------------------------------
  
  ## UMAP (RAW)
  png(
    file.path(plot_dir, paste0("04_Vision_UMAP_", ds_name, "_RAW_Score.png")),
    width = 800, height = 700
  )
  print(
    FeaturePlot(
      pbmc,
      features = "Vision_B_Cell_Signature",
      reduction = "umap",
      raster = TRUE
    ) +
      scale_colour_viridis_c(option = "magma") +
      labs(
        title = "Vision Score auf UMAP",
        subtitle = "Raw score",
        color = "Vision Score"
      ) +
      theme_minimal()
  )
  dev.off()
  
  ## Violin (RAW)
  png(
    file.path(plot_dir, paste0("04_Vision_VlnPlot_byClusters_", ds_name, "_RAW_Score.png")),
    width = 800, height = 600
  )
  print(
    VlnPlot(
      pbmc,
      features = "Vision_B_Cell_Signature",
      group.by = "seurat_clusters",
      pt.size = 0
    ) +
      labs(title = "Vision RAW Score Verteilung über Cluster") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  )
  dev.off()
  
  ## ------------------------------------------------
  ## FACETED CONFUSION DISTRIBUTION (RAW SCORE)
  ## ------------------------------------------------
  
  raw_threshold <- mean(pbmc$Vision_B_Cell_Signature, na.rm = TRUE) +
    THRESHOLD_Z * sd(pbmc$Vision_B_Cell_Signature, na.rm = TRUE)
  
  faceted_raw <- ggplot(
    pbmc@meta.data,
    aes(x = Vision_B_Cell_Signature)
  ) +
    geom_density(fill = "steelblue", alpha = 0.6) +
    facet_wrap(~ confusion, scales = "free_y") +
    geom_vline(
      xintercept = raw_threshold,
      linetype = "dashed",
      color = "red"
    ) +
    theme_classic() +
    labs(
      title = "B naive Signatur – RAW Score Verteilungen",
      subtitle = "Faceted by TP / FP / FN / TN (Z-score decision rule)",
      x = "Vision RAW Score",
      y = "Density"
    )
  
  ggsave(
    filename = file.path(plot_dir, "04_Vision_Raw_Facetted_TP_FP_FN_TN.png"),
    plot = faceted_raw,
    width = 10,
    height = 6,
    dpi = 300
  )
  
  ## ------------------------------------------------
  ## Z-SCORE PLOTS (DECISION-RELEVANT)
  ## ------------------------------------------------
  
  ## UMAP (Z)
  png(
    file.path(plot_dir, paste0("04_Vision_UMAP_", ds_name, "_ZScore.png")),
    width = 1000, height = 800
  )
  print(
    FeaturePlot(
      pbmc,
      features = "Vision_ZScore",
      reduction = "umap",
      raster = TRUE
    ) +
      scale_colour_viridis_c(option = "magma") +
      labs(
        title = "Vision Z-Score auf UMAP",
        subtitle = paste("Threshold Z =", THRESHOLD_Z),
        color = "Vision Z-Score"
      ) +
      theme_minimal()
  )
  dev.off()
  
  ## Violin (Z)
  png(
    file.path(plot_dir, paste0("04_Vision_VlnPlot_byClusters_", ds_name, "_ZScore.png")),
    width = 800, height = 600
  )
  print(
    VlnPlot(
      pbmc,
      features = "Vision_ZScore",
      group.by = "seurat_clusters",
      pt.size = 0
    ) +
      labs(title = "Vision Z-Score Verteilung über Cluster") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  )
  dev.off()
  
  ## ------------------------------------------------
  ## FACETED CONFUSION DISTRIBUTION (Z-SCORE)
  ## ------------------------------------------------
  
  faceted_z <- ggplot(
    pbmc@meta.data,
    aes(x = Vision_ZScore)
  ) +
    geom_density(fill = "steelblue", alpha = 0.6) +
    facet_wrap(~ confusion, scales = "free_y") +
    geom_vline(
      xintercept = THRESHOLD_Z,
      linetype = "dashed",
      color = "red"
    ) +
    theme_classic() +
    labs(
      title = "B naive Signatur – Z-Score Verteilungen",
      subtitle = "Faceted by TP / FP / FN / TN",
      x = "Vision Z-Score",
      y = "Density"
    )
  
  ggsave(
    filename = file.path(plot_dir, "04_Vision_ZScore_Facetted_TP_FP_FN_TN.png"),
    plot = faceted_z,
    width = 10,
    height = 6,
    dpi = 300
  )
  
  message("=== Plotting block finished successfully ===")
  
}, error = function(e) {
  
  message("!!! Plotting failed — continuing pipeline !!!")
  message("Error message:")
  message(e$message)
  
})