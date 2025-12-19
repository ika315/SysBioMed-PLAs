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

base_dir <- getwd()

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

# 1) Counts aus Seurat v5
counts <- GetAssayData(pbmc, layer = "counts")

# 2) Lineare Library-Size-Normalisierung (laut tutorial)
n_umi <- colSums(counts)
vision_expr <- t(t(counts) / n_umi) * median(n_umi)

# 3) Vision-Objekt
vision_obj <- Vision(
  data = vision_expr,
  signatures = list(sig_vision),
  meta = pbmc@meta.data,
  pool = FALSE
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
# ground truth
pbmc$true_naive_b <- pbmc$celltype.l2 == "Naive B"
#pbmc$true_naive_b <- pbmc$seurat_annotations == "B" # for small dataset

# predicted from score
# logic:
# look only at true B cells, compute the median of their Vision signature score, use this as a decision threshold
threshold <- median(pbmc$Vision_B_Cell_Signature[pbmc$true_naive_b])
pbmc$pred_naive_b <- pbmc$Vision_B_Cell_Signature >= threshold # this is saying: “a cell must score at least as high as a typical true (naive) B cell”

# confusion classes
pbmc$confusion <- with(pbmc@meta.data,
                       ifelse(true_naive_b & pred_naive_b, "TP",
                              ifelse(!true_naive_b & pred_naive_b, "FP",
                                     ifelse(true_naive_b & !pred_naive_b, "FN",
                                            "TN")))
)

table(pbmc$confusion)

tryCatch({
  
  message("=== Starting plotting block ===")
  
  ## ------------------------------------------------
  ## Ensure output directories exist
  ## ------------------------------------------------
  dir.create(file.path(base_dir, "plots", "counts_matrix"),
             recursive = TRUE, showWarnings = FALSE)
  
  ## ------------------------------------------------
  ## RAW SCORE PLOTS
  ## ------------------------------------------------
  
  ## UMAP (RAW)
  fn_umap_raw <- file.path(
    base_dir, "plots", "counts_matrix",
    paste0("04_Vision_UMAP_", ds_name, "_RAW_Score.png")
  )
  
  png(fn_umap_raw, width = 800, height = 700)
  print(
    FeaturePlot(pbmc,
                features = "Vision_B_Cell_Signature",
                reduction = "umap",
                label = TRUE,
                label.size = 6,
                repel = TRUE,
                raster = TRUE) +
      scale_colour_viridis_c(option = "magma") +
      labs(title = "Vision Score auf UMAP",
           subtitle = "RAW score",
           color = "Vision Score") +
      theme(legend.position = "right")
  )
  dev.off()
  
  ## Violin (RAW)
  fn_vln_raw <- file.path(
    base_dir, "plots", "counts_matrix",
    paste0("04_Vision_VlnPlot_byClusters_", ds_name, "_RAW_Score.png")
  )
  
  png(fn_vln_raw, width = 800, height = 600)
  print(
    VlnPlot(pbmc,
            features = "Vision_B_Cell_Signature",
            group.by = "seurat_clusters",
            pt.size = 0,
            ncol = 1) +
      labs(title = "Vision RAW Score Verteilung über Cluster") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  )
  dev.off()
  
  ## ------------------------------------------------
  ## FACETED CONFUSION DISTRIBUTION (RAW SCORE)
  ## ------------------------------------------------
  
  fn_faceted <- file.path(
    base_dir, "plots", "counts_matrix",
    "04_Vision_Raw_Facetted_TP_FP_FN_TN.png"
  )
  
  faceted <- ggplot(
    pbmc@meta.data,
    aes(x = Vision_B_Cell_Signature)
  ) +
    geom_density(fill = "steelblue", alpha = 0.6) +
    facet_wrap(~ confusion, scales = "free_y") +
    geom_vline(xintercept = threshold,
               linetype = "dashed", color = "red") +
    theme_classic() +
    labs(
      title = "B cell naive signature raw score distributions",
      subtitle = "Faceted by TP / FP / FN / TN",
      x = "B cell naive signature score",
      y = "Density"
    )
  
  ggsave(
    filename = fn_faceted,
    plot = faceted,
    width = 10,
    height = 6,
    dpi = 300
  )
  
  ## ------------------------------------------------
  ## Z-SCORE COMPUTATION (explicit)
  ## ------------------------------------------------
  
  message("Computing Vision Z-score (explicit mean/sd)")
  
  scores <- pbmc@meta.data[["Vision_B_Cell_Signature"]]
  
  vision_mean <- mean(scores, na.rm = TRUE)
  vision_sd   <- sd(scores, na.rm = TRUE)
  
  pbmc$Vision_ZScore <- (scores - vision_mean) / vision_sd
  
  ## ------------------------------------------------
  ## Z-SCORE PLOTS
  ## ------------------------------------------------
  
  ## UMAP (Z)
  fn_umap_z <- file.path(
    base_dir, "plots", "counts_matrix",
    paste0("04_Vision_UMAP_", ds_name, "_ZScore.png")
  )
  
  png(fn_umap_z, width = 1000, height = 800)
  print(
    FeaturePlot(pbmc,
                features = "Vision_ZScore",
                reduction = "umap",
                label = TRUE,
                label.size = 6,
                repel = TRUE,
                raster = TRUE) +
      scale_colour_viridis_c(option = "magma") +
      labs(title = "Vision Z-Score auf UMAP",
           subtitle = "standardized score",
           color = "Vision Z-Score") +
      theme(legend.position = "right")
  )
  dev.off()
  
  ## Violin (Z)
  fn_vln_z <- file.path(
    base_dir, "plots", "counts_matrix",
    paste0("04_Vision_VlnPlot_byClusters_", ds_name, "_ZScore.png")
  )
  
  png(fn_vln_z, width = 800, height = 600)
  print(
    VlnPlot(pbmc,
            features = "Vision_ZScore",
            group.by = "seurat_clusters",
            pt.size = 0,
            ncol = 1) +
      labs(title = "Vision Z-Score Verteilung über Cluster") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  )
  dev.off()
  
  ## ------------------------------------------------
  ## FACETED CONFUSION DISTRIBUTION (Z-SCORE)
  ## ------------------------------------------------
  
  fn_faceted_z <- file.path(
    base_dir, "plots", "counts_matrix",
    "04_Vision_ZScore_Facetted_TP_FP_FN_TN.png"
  )
  
  z_threshold <- (threshold - vision_mean) / vision_sd
  
  faceted_z <- ggplot(
    pbmc@meta.data,
    aes(x = Vision_ZScore)
  ) +
    geom_density(fill = "steelblue", alpha = 0.6) +
    facet_wrap(~ confusion, scales = "free_y") +
    geom_vline(
      xintercept = z_threshold,
      linetype = "dashed",
      color = "red"
    ) +
    theme_classic() +
    labs(
      title = "B cell naive signature Z-score distributions",
      subtitle = "Faceted by TP / FP / FN / TN",
      x = "B cell naive signature Z-score",
      y = "Density"
    )
  
  ggsave(
    filename = fn_faceted_z,
    plot = faceted_z,
    width = 10,
    height = 6,
    dpi = 300
  )
  
  message("=== Plotting block finished successfully ===")
  
}, error = function(e) {
  
  message("!!! Plotting block failed — continuing pipeline !!!")
  message("Error message:")
  message(e$message)
  
})