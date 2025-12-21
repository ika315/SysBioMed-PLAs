# -------------------------------------------------------------------
# 04_Benchmarking.R
# Benchmarking von PLA-Methoden anhand Platelet-Zell Signatur
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
path_gene_list = "~/SysBioMed-PLAs/data/updated_gene_list.csv"
source("~/SysBioMed-PLAs/src/read_and_extend_gene_list.R")
genes <- read_gene_list(path_gene_list)
sig_vision <- make_signature(genes = genes, sig_name = "Platelet_Signature", method = "vision")

pbmc$GT_Response <- ifelse(pbmc[[gt_col_name]] == "Platelet", 1, 0)

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
  signatures = list(sig_vision),
  pool = TRUE
)

vision_obj <- analyze(vision_obj)

vision_scores <- getSignatureScores(vision_obj)[, "Platelet_Signature"]

score_name <- "Vision_Platelet_Signature"

pbmc$Vision_Raw <- as.numeric(vision_scores)
pbmc$Vision_ZScore <- as.vector(scale(pbmc$Vision_Raw))

message("=== Extending platelet gene set using VISION scores ===")

# (optional but recommended) restrict background
pbmc_ext <- subset(
  pbmc,
  celltype.l1 %in% c("Platelet", "Myeloid", "Erythroid")
)

res_ext <- extend_gene_set(
  pbmc       = pbmc_ext,
  base_genes = genes,
  score_name = "Vision_Raw",   # IMPORTANT
  top_n      = 50,
  high_quantile = 0.9,
  min.pct    = 0.05,
  test.use   = "wilcox"
)

extended_genes <- res_ext$extended_genes
top_genes      <- res_ext$top_genes
markers_ext    <- res_ext$markers

message(paste("Original genes:", length(genes)))
message(paste("Extended genes:", length(extended_genes)))

# Save extended gene list
write.csv(
  data.frame(geneName = extended_genes),
  "results/extended_platelet_gene_list.csv",
  row.names = FALSE,
  quote = FALSE
)

sig_vision_ext <- make_signature(
  genes = extended_genes,
  sig_name = "Platelet_Signature_Extended",
  method = "vision"
)

message("--- Berechne Vision Score (Extended Signature) ---")

vision_ext <- Vision(
  pbmc,
  assay = "RNA_v3",
  signatures = list(sig_vision_ext),
  pool = TRUE
)

vision_ext <- analyze(vision_ext)

pbmc$Vision_Raw_Ext <- getSignatureScores(vision_ext)[, "Platelet_Signature_Extended"]
pbmc$Vision_ZScore_Ext <- as.vector(scale(pbmc$Vision_Raw_Ext))

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

pbmc$Prediction_Ext <- ifelse(
  pbmc$Vision_ZScore_Ext > THRESHOLD_Z,
  "Positive",
  "Negative"
)

pbmc$Error_Type_Ext <- dplyr::case_when(
  pbmc$Prediction_Ext == "Positive" & pbmc$GT_Class == "Positive" ~ "TP",
  pbmc$Prediction_Ext == "Positive" & pbmc$GT_Class == "Negative" ~ "FP",
  pbmc$Prediction_Ext == "Negative" & pbmc$GT_Class == "Positive" ~ "FN",
  pbmc$Prediction_Ext == "Negative" & pbmc$GT_Class == "Negative" ~ "TN"
)

ConfusionMatrix_Ext <- table(
  Predicted = pbmc$Prediction_Ext,
  Actual = pbmc$GT_Class
)

print("=== Extended Signature Confusion Matrix ===")
print(ConfusionMatrix_Ext)

TP_ext <- ConfusionMatrix_Ext["Positive", "Positive"]
FP_ext <- ConfusionMatrix_Ext["Positive", "Negative"]
FN_ext <- ConfusionMatrix_Ext["Negative", "Positive"]
TN_ext <- ConfusionMatrix_Ext["Negative", "Negative"]

Precision_ext   <- TP_ext / (TP_ext + FP_ext)
Recall_ext      <- TP_ext / (TP_ext + FN_ext)
Specificity_ext <- TN_ext / (TN_ext + FP_ext)
Accuracy_ext    <- (TP_ext + TN_ext) / sum(ConfusionMatrix_Ext)
F1_Score_ext    <- 2 * (Precision_ext * Recall_ext) / (Precision_ext + Recall_ext)

## ------------------------------------------------
## PRINT EVALUATION METRICS
## ------------------------------------------------
message("=== Vision Benchmarking Metrics (Platelet) ===")

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

metrics_ext_df <- data.frame(
  Metric = c("Precision", "Recall", "Specificity", "Accuracy", "F1"),
  Value  = c(
    Precision_ext,
    Recall_ext,
    Specificity_ext,
    Accuracy_ext,
    F1_Score_ext
  )
)

print(metrics_ext_df)

message(
  sprintf(
    "EXTENDED | Precision: %.3f | Recall: %.3f | F1: %.3f | Accuracy: %.3f",
    Precision_ext,
    Recall_ext,
    F1_Score_ext,
    Accuracy_ext
  )
)

## ------------------------------------------------
## SAVE FINAL PBMC OBJECT (WITH VISION SCORES)
## ------------------------------------------------
dir.create("results", recursive = TRUE, showWarnings = FALSE)

pbmc_out_path <- file.path(
  "results",
  paste0("pbmc_", DATASET_NAME, "_Vision_Platelet.rds")
)

saveRDS(pbmc, pbmc_out_path)

message("Saved final pbmc object to:")
message(pbmc_out_path)

tryCatch({
  
  message("=== Starting plotting block ===")
  
  dir.create("plots", recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path("plots", "benchmarking_vision"), recursive = TRUE, showWarnings = FALSE)
  
  annotation_levels <- c("celltype.l1", "celltype.l2", "celltype.l3")
  
  score_variants <- list(
    original = list(raw = "Vision_Raw",     z = "Vision_ZScore"),
    extended = list(raw = "Vision_Raw_Ext", z = "Vision_ZScore_Ext")
  )
  
  ## ------------------------------------------------
  ## REFERENCE UMAP (GROUND TRUTH) — unchanged
  ## ------------------------------------------------
  for (ann in annotation_levels) {
    
    if (!ann %in% colnames(pbmc@meta.data)) next
    
    png(
      file.path("plots", "benchmarking_vision",
                paste0("04_Reference_UMAP_Platelet_", ann, ".png")),
      1200, 800
    )
    
    print(
      DimPlot(
        pbmc, group.by = ann, reduction = "umap",
        label = TRUE, repel = TRUE, label.size = 4, raster = TRUE
      ) +
        labs(
          title = paste("Referenz-UMAP:", ann),
          subtitle = "Ground truth annotation",
          color = ann
        ) +
        theme_minimal()
    )
    
    dev.off()
  }
  
  ## ------------------------------------------------
  ## SCORE-DEPENDENT PLOTS (ORIGINAL vs EXTENDED)
  ## ------------------------------------------------
  for (variant in names(score_variants)) {
    
    raw_col <- score_variants[[variant]]$raw
    z_col   <- score_variants[[variant]]$z
    
    message(paste("Plotting", variant, "VISION results"))
    
    for (ann in annotation_levels) {
      
      if (!ann %in% colnames(pbmc@meta.data)) next
      Idents(pbmc) <- ann
      
      ## ---- UMAP ----
      png(
        file.path(
          "plots", "benchmarking_vision",
          paste0("Vision_Platelet_", variant, "_UMAP_", ann, ".png")
        ),
        900, 700
      )
      
      p_umap <- FeaturePlot(
        pbmc,
        features = raw_col,
        reduction = "umap",
        raster = TRUE
      ) +
        scale_colour_viridis_c(option = "magma") +
        labs(
          title = paste("VISION Score (Platelet,", variant, ")"),
          subtitle = paste("Annotation:", ann),
          color = "VISION Score"
        ) +
        theme_minimal()
      
      print(
        p_umap +
          LabelClusters(pbmc, id = ann, repel = TRUE, size = 4)
      )
      
      dev.off()
      
      ## ---- VIOLIN ----
      png(
        file.path(
          "plots", "benchmarking_vision",
          paste0("Vision_Platelet_", variant, "_Violin_", ann, ".png")
        ),
        1200, 800
      )
      
      print(
        VlnPlot(
          pbmc,
          features = raw_col,
          group.by = ann,
          pt.size = 0
        ) +
          labs(
            title = paste("VISION Score (Platelet,", variant, ") –", ann)
          ) +
          theme(axis.text.x = element_text(angle = 45, hjust = 1))
      )
      
      dev.off()
    }
    
    ## ---- PRECISION–RECALL ----
    eval_df <- data.frame(
      score = pbmc[[z_col]],
      gt    = pbmc$GT_Response
    ) |>
      arrange(desc(score)) |>
      mutate(
        tp_cum = cumsum(gt),
        fp_cum = cumsum(1 - gt),
        precision = tp_cum / (tp_cum + fp_cum),
        recall    = tp_cum / sum(gt)
      )
    
    p_pr <- ggplot(eval_df, aes(recall, precision)) +
      geom_line(color = "#E41A1C", linewidth = 1.2) +
      theme_minimal() +
      labs(
        title = paste("Precision–Recall (Platelet,", variant, ")"),
        x = "Recall", y = "Precision"
      )
    
    ggsave(
      file.path(
        "plots", "benchmarking_vision",
        paste0("Vision_Platelet_", variant, "_PrecisionRecall.png")
      ),
      p_pr, width = 7, height = 6
    )
    
    ## ---- Z-SCORE DISTRIBUTION ----
    p_z <- ggplot(
      pbmc@meta.data,
      aes(x = .data[[gt_col_name]], y = .data[[z_col]],
          fill = .data[[gt_col_name]])
    ) +
      geom_violin(alpha = 0.7) +
      geom_boxplot(width = 0.1, outlier.shape = NA) +
      geom_hline(yintercept = THRESHOLD_Z,
                 linetype = "dashed", color = "red") +
      theme_minimal() +
      labs(
        title = paste("VISION Z-Score (Platelet,", variant, ")"),
        x = "Cell type", y = "Z-score"
      ) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            legend.position = "none")
    
    ggsave(
      file.path(
        "plots", "benchmarking_vision",
        paste0("Vision_Platelet_", variant, "_ZScore_Distribution.png")
      ),
      p_z, width = 12, height = 6
    )
    
    ## ---- ERROR TYPE ----
    error_col <- ifelse(variant == "original", "Error_Type", "Error_Type_Ext")
    
    p_err <- ggplot(
      pbmc@meta.data,
      aes(.data[[error_col]], .data[[z_col]], fill = .data[[error_col]])
    ) +
      geom_boxplot(outlier.shape = NA) +
      theme_minimal() +
      labs(
        title = paste("VISION Z-Score by Error Type (", variant, ")", sep = "")
      )
    
    ggsave(
      file.path(
        "plots", "benchmarking_vision",
        paste0("Vision_Platelet_", variant, "_ErrorTypes.png")
      ),
      p_err, width = 8, height = 6
    )
  }
  
  message("=== Plotting block finished successfully ===")
  
}, error = function(e) {
  
  message("!!! Plotting failed — continuing pipeline !!!")
  message(e$message)
  
})