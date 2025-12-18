# -------------------------------------------------------------------
# 04_Benchmarking.R
# Benchmarking von PLA-Methoden anhand B naive-Zell Signatur
# Für Vision
# -------------------------------------------------------------------

library(Seurat)
library(VISION)
library(pROC)
library(dplyr)
library(ggplot2)

options(mc.cores = as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", 1)))

PATH_DATA <- "~/SysBioMed-PLAs/data/seu_sx_final.rds"
DATASET_NAME <- "seu_sx_final"

print(paste("Lade prozessiertes Seurat-Objekt:", PATH_DATA))
pbmc <- readRDS(PATH_DATA)

# Ground truth (GT) Definition
gt_col_name <- "celltype.l2"

# Memory B-Zell Signatur
path_gene_list = "~/SysBioMed-PLAs/data/naive_b_cells.csv"
source("~/SysBioMed-PLAs/src/read_and_extend_gene_list.R")
genes <- read_gene_list(path_gene_list)
sig_vision <- make_signature(genes = genes, sig_name = "Naive_B_Cell_Signature", method = "vision")

pbmc$GT_Response <- ifelse(pbmc[[gt_col_name]] == "B naive", 1, 0)


# --- VISION SCORING ---
print("--- Berechne Vision Score ---")
pbmc$celltype.l2 <- factor(pbmc$celltype.l2)

pbmc[["RNA_v3"]] <- as(pbmc[["RNA"]], "Assay")
DefaultAssay(pbmc) <- "RNA_v3"

vision_obj <- Vision(
  pbmc,
  assay = "RNA_v3",
  signatures = list(sig_vision),
  pool = FALSE
)

vision_obj <- analyze(vision_obj)

vision_auc <- getSignatureScores(vision_obj)
vision_score <- vision_auc[, "Naive_B_Cell_Signature"]
vision_score <- as.numeric(vision_auc[, "Naive_B_Cell_Signature"])

pbmc$Vision_Naive_B_Cell_Signature <- vision_score
score_name <- "Vision_Naive_B_Cell_Signature"

pbmc$Vision_Raw <- pbmc$Vision_Naive_B_Cell_Signature

tryCatch({
  
  message("=== Starting plotting block ===")
  
  ## ------------------------------------------------
  ## Ensure output directories exist
  ## ------------------------------------------------
  dir.create(file.path("plots"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path("plots", "seurat_downgrade"), recursive = TRUE, showWarnings = FALSE)
  
  ## ------------------------------------------------
  ## UMAP + VIOLIN (RAW SCORE) for l1 / l2 / l3
  ## ------------------------------------------------
  
  annotation_levels <- c("celltype.l1", "celltype.l2", "celltype.l3")
  
  for (ann in annotation_levels) {
    
    if (!ann %in% colnames(pbmc@meta.data)) {
      message(paste("Skipping", ann, "- column not found"))
      next
    }
    
    message(paste("Plotting Vision RAW score for", ann))
    
    ## ---- UMAP ----
    fn_umap <- file.path(
      "plots",
      paste0("Vision_NaiveB_UMAP_", ann, ".png")
    )
    
    png(fn_umap, width = 900, height = 700)
    
    print(
      FeaturePlot(
        pbmc,
        features = "Vision_Raw",
        reduction = "umap",
        label = TRUE,
        label.size = 6,
        repel = TRUE,
        raster = TRUE
      ) +
        scale_colour_viridis_c(option = "magma") +
        labs(
          title = "Vision Score (Naive B) auf UMAP",
          subtitle = paste("Legende zeigt", ann),
          color = "Vision Score"
        ) +
        theme(legend.position = "right")
    )
    dev.off()
    
    ## ---- VIOLIN ----
    fn_vln <- file.path(
      "plots",
      paste0("Vision_NaiveB_Violin_", ann, ".png")
    )
    
    png(fn_vln, width = 1000, height = 600)
    
    print(
      VlnPlot(
        pbmc,
        features = "Vision_Raw",
        group.by = ann,
        pt.size = 0
      ) +
        labs(
          title = paste("Vision Score Verteilung über Zelllinie", ann)
        ) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
    )
    dev.off()
  }
  
  ## ------------------------------------------------
  ## PRECISION–RECALL CURVE
  ## ------------------------------------------------
  
  fn_pr <- file.path("plots", "Vision_NaiveB_PrecisionRecall_Curve.png")
  
  png(fn_pr, width = 800, height = 700)
  print(p_pr)
  dev.off()
  
  ## ------------------------------------------------
  ## Z-SCORE DISTRIBUTION (GT)
  ## ------------------------------------------------
  
  fn_zdist <- file.path("plots", "Vision_NaiveB_ZScore_Distribution.png")
  
  png(fn_zdist, width = 1500, height = 800)
  print(p_z)
  dev.off()
  
  ## ------------------------------------------------
  ## Z-SCORE 3-GROUP VIOLIN
  ## ------------------------------------------------
  
  fn_z3 <- file.path("plots", "Vision_NaiveB_ZScore_3Groups.png")
  
  png(fn_z3, width = 800, height = 600)
  print(p_z3)
  dev.off()
  
  ## ------------------------------------------------
  ## ERROR TYPE BOXPLOTS
  ## ------------------------------------------------
  
  fn_err1 <- file.path("plots", "04_Bench_ErrorTypes_NaiveB.png")
  png(fn_err1, width = 900, height = 650)
  print(p_error)
  dev.off()
  
  fn_err2 <- file.path("plots", "seurat_downgrade", "04_Vision_Bench_ErrorTypes_B.png")
  png(fn_err2, width = 900, height = 650)
  print(p_error)
  dev.off()
  
  ## ------------------------------------------------
  ## FACETED CONFUSION DISTRIBUTIONS
  ## ------------------------------------------------
  
  fn_fac_z <- file.path(
    "plots", "seurat_downgrade",
    "04_Vision_ZScore_Facetted_TP_FP_FN_TN.png"
  )
  
  ggsave(
    filename = fn_fac_z,
    plot = faceted_z,
    width = 10,
    height = 6,
    dpi = 300
  )
  
  fn_fac_raw <- file.path(
    "plots", "seurat_downgrade",
    "04_Vision_Raw_Facetted_TP_FP_FN_TN.png"
  )
  
  ggsave(
    filename = fn_fac_raw,
    plot = faceted_raw,
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