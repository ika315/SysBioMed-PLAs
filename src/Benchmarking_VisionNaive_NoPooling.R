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

#options(mc.cores = as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", 1)))
options(mc.cores = 1)

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

pbmc$celltype.l1 <- factor(pbmc$celltype.l1)
pbmc$celltype.l2 <- factor(pbmc$celltype.l2)
pbmc$celltype.l3 <- factor(pbmc$celltype.l3)


# --- VISION SCORING ---
print("--- Berechne Vision Score ---")

pbmc[["RNA_v3"]] <- as(pbmc[["RNA"]], "Assay")
DefaultAssay(pbmc) <- "RNA_v3"

vision_obj <- Vision(
  pbmc,
  assay = "RNA_v3",
  signatures = list(sig_vision)
  #,pool = TRUE
)

vision_obj <- analyze(vision_obj)

vision_scores <- getSignatureScores(vision_obj)[, "Naive_B_Cell_Signature"]

score_name <- "Vision_Naive_B_Cell_Signature"

pbmc$Vision_Raw <- as.numeric(vision_scores)
pbmc$Vision_ZScore <- as.vector(scale(pbmc$Vision_Raw))

THRESHOLD_Z <- 1.5

pbmc$Prediction <- ifelse(pbmc$Vision_ZScore > THRESHOLD_Z, "Positive", "Negative")
pbmc$GT_Class <- ifelse(pbmc$GT_Response == 1, "Positive", "Negative")

ConfusionMatrix <- table(
  Predicted = pbmc$Prediction,
  Actual = pbmc$GT_Class
)

TP <- ConfusionMatrix["Positive", "Positive"]
FP <- ConfusionMatrix["Positive", "Negative"]
FN <- ConfusionMatrix["Negative", "Positive"]
TN <- ConfusionMatrix["Negative", "Negative"]

Precision   <- TP / (TP + FP)
Recall      <- TP / (TP + FN)
Specificity <- TN / (TN + FP)
Accuracy    <- (TP + TN) / sum(ConfusionMatrix)
F1_Score    <- 2 * (Precision * Recall) / (Precision + Recall)

pbmc$Error_Type <- dplyr::case_when(
  pbmc$Prediction == "Positive" & pbmc$GT_Class == "Positive" ~ "TP",
  pbmc$Prediction == "Positive" & pbmc$GT_Class == "Negative" ~ "FP",
  pbmc$Prediction == "Negative" & pbmc$GT_Class == "Positive" ~ "FN",
  pbmc$Prediction == "Negative" & pbmc$GT_Class == "Negative" ~ "TN"
)

## ------------------------------------------------
## PRINT EVALUATION METRICS
## ------------------------------------------------
message("=== Vision Benchmarking Metrics (Naive B) ===")

print(ConfusionMatrix)

metrics_df <- data.frame(
  Metric = c("Precision", "Recall", "Specificity", "Accuracy", "F1"),
  Value = c(Precision, Recall, Specificity, Accuracy, F1_Score)
)

print(metrics_df)

message(
  sprintf(
    "Precision: %.3f | Recall: %.3f | F1: %.3f | Accuracy: %.3f",
    Precision, Recall, F1_Score, Accuracy
  )
)

## ------------------------------------------------
## SAVE FINAL PBMC OBJECT (WITH VISION SCORES)
## ------------------------------------------------
dir.create("results", recursive = TRUE, showWarnings = FALSE)

pbmc_out_path <- file.path(
  "results",
  paste0("pbmc_", DATASET_NAME, "_Vision_NaiveB.rds")
)

saveRDS(pbmc, pbmc_out_path)

message("Saved final pbmc object to:")
message(pbmc_out_path)

tryCatch({
  
  message("=== Starting plotting block ===")
  
  ## ------------------------------------------------
  ## Ensure output directories exist
  ## ------------------------------------------------
  dir.create("plots", recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path("plots", "benchmarking_vision"), recursive = TRUE, showWarnings = FALSE)
  
  ## ------------------------------------------------
  ## REFERENCE UMAP (GROUND TRUTH)
  ## ------------------------------------------------
  annotation_levels <- c("celltype.l1", "celltype.l2", "celltype.l3")
  
  for (ann in annotation_levels) {
    
    if (!ann %in% colnames(pbmc@meta.data)) {
      message(paste("Skipping reference UMAP for", ann, "- column not found"))
      next
    }
    
    png(
      filename = file.path("plots", "benchmarking_vision", paste0("04_Reference_UMAP_NaiveB_", ann, ".png")),
      width = 1200,
      height = 800
    )
    
    print(
      DimPlot(
        pbmc,
        group.by = ann,
        reduction = "umap",
        label = TRUE,
        repel = TRUE,
        label.size = 4,
        raster = TRUE
      ) +
        labs(
          title = paste("Referenz-UMAP:", ann),
          subtitle = "Basierend auf der originalen Annotation (Ground Truth)",
          color = ann
        ) +
        theme_minimal() +
        theme(legend.position = "right")
    )
    
    dev.off()
  }
  
  ## ------------------------------------------------
  ## UMAP + VIOLIN (RAW SCORE)
  ## ------------------------------------------------
  for (ann in annotation_levels) {
    
    if (!ann %in% colnames(pbmc@meta.data)) {
      message(paste("Skipping", ann, "- column not found"))
      next
    }
    
    message(paste("Plotting Vision RAW score for", ann))
    
    ## set identities
    Idents(pbmc) <- ann
    
    ## ---- UMAP ----
    png(
      filename = file.path(
        "plots", "benchmarking_vision",
        paste0("Vision_NaiveB_UMAP_", ann, ".png")
      ),
      width = 900, height = 700
    )
    
    p <- FeaturePlot(
      pbmc,
      features = "Vision_Raw",
      reduction = "umap",
      raster = TRUE
    ) +
      scale_colour_viridis_c(option = "magma") +
      labs(
        title = "VISION Score (Naive B) auf UMAP",
        subtitle = paste("Annotation:", ann),
        color = "VISION Score"
      ) +
      theme_minimal()
    
    p_labeled <- LabelClusters(
      plot  = p,
      id    = ann,
      repel = TRUE,
      size  = 4
    )
    
    print(p_labeled)
    
    dev.off()
    
    ## ---- VIOLIN ----
    png(
      filename = file.path(
        "plots", "benchmarking_vision",
        paste0("Vision_NaiveB_Violin_", ann, ".png")
      ),
      width = 1200, height = 800
    )
    
    print(
      VlnPlot(
        pbmc,
        features = "Vision_Raw",
        group.by = ann,
        pt.size = 0
      ) +
        labs(title = paste("VISION Score Verteilung (B Naive) –", ann)) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
    )
    
    dev.off()
  }
  
  ## ------------------------------------------------
  ## PRECISION–RECALL CURVE
  ## ------------------------------------------------
  eval_df <- data.frame(
    score = pbmc$Vision_ZScore,
    gt = pbmc$GT_Response
  ) %>%
    arrange(desc(score)) %>%
    mutate(
      tp_cum = cumsum(gt),
      fp_cum = cumsum(1 - gt),
      precision_vec = tp_cum / (tp_cum + fp_cum),
      recall_vec = tp_cum / sum(gt)
    )
  
  p_pr <- ggplot(eval_df, aes(x = recall_vec, y = precision_vec)) +
    geom_line(color = "#E41A1C", linewidth = 1.2) +
    labs(
      title = "Precision–Recall Kurve (VISION, B Naive)",
      x = "Recall",
      y = "Precision"
    ) +
    theme_minimal()
  
  png("plots/benchmarking_vision/Vision_NaiveB_PrecisionRecall_Curve.png", 800, 700)
  print(p_pr)
  dev.off()
  
  ## ------------------------------------------------
  ## Z-SCORE DISTRIBUTION (GT)
  ## ------------------------------------------------
  p_z <- ggplot(
    pbmc@meta.data,
    aes(x = .data[[gt_col_name]], y = Vision_ZScore, fill = .data[[gt_col_name]])
  ) +
    geom_violin(alpha = 0.7) +
    geom_boxplot(width = 0.1, outlier.shape = NA) +
    geom_hline(yintercept = THRESHOLD_Z, linetype = "dashed", color = "red") +
    labs(
      title = "VISION Z-Score Verteilung über Zelltypen (B Naive)",
      x = "Zelltyp",
      y = "VISION Z-Score"
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none")
  
  png("plots/benchmarking_vision/Vision_NaiveB_ZScore_Distribution.png", 1500, 800)
  print(p_z)
  dev.off()
  
  ## ------------------------------------------------
  ## Z-SCORE 3-GROUP VIOLIN
  ## ------------------------------------------------
  pbmc$ZScore_Groups <- factor(
    dplyr::case_when(
      pbmc$celltype.l2 == "B naive" ~ "B naive",
      pbmc$celltype.l1 == "B" ~ "Rest B-cells",
      TRUE ~ "Other"
    ),
    levels = c("B naive", "Rest B-cells", "Other")
  )
  
  p_z3 <- ggplot(
    pbmc@meta.data,
    aes(x = ZScore_Groups, y = Vision_ZScore, fill = ZScore_Groups)
  ) +
    geom_violin(alpha = 0.7) +
    geom_boxplot(width = 0.1, outlier.shape = NA) +
    geom_hline(yintercept = THRESHOLD_Z, linetype = "dashed", color = "red") +
    theme_minimal()
  
  png("plots/benchmarking_vision/Vision_NaiveB_ZScore_3Groups.png", 800, 600)
  print(p_z3)
  dev.off()
  
  ## ------------------------------------------------
  ## ERROR TYPE BOXPLOT
  ## ------------------------------------------------
  p_error <- ggplot(
    pbmc@meta.data,
    aes(x = Error_Type, y = Vision_ZScore, fill = Error_Type)
  ) +
    geom_boxplot(outlier.shape = NA) +
    theme_minimal() +
    labs(title = "VISION Z-Score nach Fehlerklassen")
  
  png("plots/benchmarking_vision/04_Vision_Bench_ErrorTypes_NaiveB.png", 900, 650)
  print(p_error)
  dev.off()
  
  ## ------------------------------------------------
  ## FACETED CONFUSION DISTRIBUTIONS
  ## ------------------------------------------------
  faceted_z <- ggplot(pbmc@meta.data, aes(x = Vision_ZScore)) +
    geom_density(fill = "steelblue", alpha = 0.6) +
    facet_wrap(~ Error_Type, scales = "free_y") +
    geom_vline(xintercept = THRESHOLD_Z, linetype = "dashed", color = "red") +
    theme_classic()
  
  faceted_raw <- ggplot(pbmc@meta.data, aes(x = Vision_Raw)) +
    geom_density(fill = "orange", alpha = 0.6) +
    facet_wrap(~ Error_Type, scales = "free_y") +
    theme_classic()
  
  ggsave(
    "plots/benchmarking_vision/04_Vision_ZScore_NaiveB_Facetted_TP_FP_FN_TN.png",
    faceted_z, width = 10, height = 6, dpi = 300
  )
  
  ggsave(
    "plots/benchmarking_vision/04_Vision_Raw_NaiveB_Facetted_TP_FP_FN_TN.png",
    faceted_raw, width = 10, height = 6, dpi = 300
  )
  
  message("=== Plotting block finished successfully ===")
  
}, error = function(e) {
  
  message("!!! Plotting failed — continuing pipeline !!!")
  message("Error message:")
  message(e$message)
  
})