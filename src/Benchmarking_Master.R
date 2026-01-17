# -------------------------------------------------------------------
# Benchmarking_Master.R
# Unified script for PLA benchmarking (AUCell, UCell, AddModuleScore)
# -------------------------------------------------------------------

library(Seurat)
library(AUCell)
library(UCell)
library(pROC)
library(dplyr)
library(ggplot2)


# Configuration
# Options: "AUCell", "UCell", "AddModuleScore"
METHOD_NAME  <- "AddModuleScore" 
TARGET_LABEL <- "PLA_Gating"
GT_COLUMN    <- "pla.status"
POSITIVE_VAL <- "PLA"
PATH_DATA    <- "~/SysBioMed-PLAs/data/seu_sx_final.rds"

# Load data and Meta data Sync
print(paste("Running Benchmark for Method:", METHOD_NAME))
pbmc <- readRDS(PATH_DATA)

new_metadata <- read.csv("~/SysBioMed-PLAs/data/external_dataset_pla-status_metatable.csv", row.names = 1)
rownames(new_metadata) <- new_metadata$barcodes_clean
common_cells <- intersect(Cells(pbmc), rownames(new_metadata))
pbmc <- subset(pbmc, cells = common_cells)
pbmc <- AddMetaData(pbmc, metadata = new_metadata[common_cells, ])

# Define Output Directory
OUT_DIR <- paste0("plots/Platelet_Main/", METHOD_NAME, "/")
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)
if (!dir.exists("results")) dir.create("results")

# Load Gene Lists
base_dir <- getwd()
source(file.path(base_dir, "src", "read_and_extend_gene_list.R"))
genes <- read_gene_list(file.path(base_dir, "data", "updated_gene_list.csv"))

pbmc$GT_Response <- ifelse(pbmc[[GT_COLUMN]] == POSITIVE_VAL, 1, 0)
pbmc$GT_Class    <- ifelse(pbmc$GT_Response == 1, "Positive", "Negative")

# Scoring Logic
print(paste("--- Calculating Scores using", METHOD_NAME, "---"))

if (METHOD_NAME == "AUCell") {
    expression_matrix <- GetAssayData(pbmc, layer = "data")
    rankings <- AUCell_buildRankings(expression_matrix, plotStats=FALSE)
    auc_orig <- AUCell_calcAUC(list(Platelet_Orig = genes), rankings)
    pbmc$Raw_Score_Original <- as.numeric(getAUC(auc_orig)[1, ])
} else if (METHOD_NAME == "UCell") {
    pbmc <- AddModuleScore_UCell(pbmc, features = list(Platelet_Orig = genes), name = NULL)
    pbmc$Raw_Score_Original <- pbmc$Platelet_Orig
} else if (METHOD_NAME == "AddModuleScore") {
    pbmc <- AddModuleScore(pbmc, features = list(genes), name = "AMS_Orig")
    pbmc$Raw_Score_Original <- pbmc$AMS_Orig1
}

# Geneset Extension
res_ext <- extend_gene_set(
    pbmc = pbmc, base_genes = genes, score_name = "Raw_Score_Original", 
    top_n = 50, high_quantile = 0.9, min.pct = 0.05, test.use = "wilcox"
)
extended_genes <- res_ext$extended_genes
write.csv(data.frame(geneName = extended_genes), 
          paste0("results/extended_platelet_gene_list_", METHOD_NAME, ".csv"), row.names = FALSE)

# Final Scoring
if (METHOD_NAME == "AUCell") {
    auc_final <- AUCell_calcAUC(list(Platelet_Score = extended_genes), rankings)
    pbmc$Raw_Score <- as.numeric(getAUC(auc_final)[1, ])
} else if (METHOD_NAME == "UCell") {
    pbmc <- AddModuleScore_UCell(pbmc, features = list(Platelet_Score = extended_genes), name = NULL, assay = "RNA")
    pbmc$Raw_Score <- pbmc$Platelet_Score
} else if (METHOD_NAME == "AddModuleScore") {
    pbmc <- AddModuleScore(pbmc, features = list(extended_genes), name = "AMS_Final")
    pbmc$Raw_Score <- pbmc$AMS_Final1
}

pbmc$Z_Score <- as.vector(scale(pbmc$Raw_Score))

# Metrics and Threshold
roc_obj <- roc(response = pbmc$GT_Response, predictor = pbmc$Z_Score, direction = "<", quiet = TRUE)
THRESHOLD_Z <- as.numeric(coords(roc_obj, x = "best", best.method = "youden")$threshold)
pbmc$Prediction <- ifelse(pbmc$Z_Score > THRESHOLD_Z, "Positive", "Negative")

pbmc$Error_Type <- case_when(
    pbmc$Prediction == "Positive" & pbmc$GT_Class == "Positive" ~ "TP",
    pbmc$Prediction == "Positive" & pbmc$GT_Class == "Negative" ~ "FP",
    pbmc$Prediction == "Negative" & pbmc$GT_Class == "Positive" ~ "FN",
    pbmc$Prediction == "Negative" & pbmc$GT_Class == "Negative" ~ "TN"
)
pbmc$Error_Type <- factor(pbmc$Error_Type, levels = c("TP", "FP", "FN", "TN"))

# Visualizations

# --- Reference UMAP: Gated Cell Types ---
png(filename = paste0(OUT_DIR, "Reference_UMAP_Gated_Celltypes.png"), width = 1200, height = 800)
p_ref <- DimPlot(pbmc, group.by = "celltype_clean", reduction = "umap", 
                  label = TRUE, repel = TRUE, raster = TRUE) +
    labs(title = "Reference UMAP: Gated Cell Types",
         subtitle = "Based on cleaned ground truth annotation",
         color = "Cell Type") +
    theme_minimal()
print(p_ref)
dev.off()

# --- UMAP Raw Score ---
png(filename = paste0(OUT_DIR, METHOD_NAME, "_UMAP_Raw_Scores.png"), width = 900, height = 700)
p_umap_score <- FeaturePlot(pbmc, features = "Raw_Score", reduction = "umap", raster = TRUE) +
    scale_colour_viridis_c(option = "magma") +
    labs(title = paste(METHOD_NAME, "PLA Signature Score"), 
         subtitle = "Intensity of platelet signature in leukocytes") +
    theme_minimal()
print(p_umap_score)
dev.off()

# --- UMAP Error Mapping ---
png(filename = paste0(OUT_DIR, METHOD_NAME, "_Error_UMAP.png"), width = 1000, height = 800)
print(DimPlot(pbmc, group.by = "Error_Type", reduction = "umap", raster = TRUE) +
    scale_color_manual(values = c("TP"="#228B22", "FP"="#FF4500", "FN"="#1E90FF", "TN"="#D3D3D3")) +
    labs(title = paste(METHOD_NAME, "Classification Error Mapping"), 
         subtitle = paste("Threshold Z =", round(THRESHOLD_Z, 3), "| GT:", TARGET_LABEL)) +
    theme_minimal())
dev.off()

# --- Violin: Raw Scores ---
png(filename = paste0(OUT_DIR, METHOD_NAME, "_Violin_Celltypes.png"), width = 1200, height = 600)
p_vln_clean <- VlnPlot(pbmc, 
                       features = "Raw_Score", 
                       group.by = "celltype_clean", 
                       pt.size = 0) + # No points for clarity
    labs(title = paste(METHOD_NAME, "Raw Score Distribution"),
         subtitle = paste("Grouped by celltype_clean | Ground Truth:", TARGET_LABEL),
         x = "Cell Type (Gated)",
         y = "Raw Score") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none")
print(p_vln_clean)
dev.off()

# Split Violin: PLA vs. Platelet-free
png(filename = paste0(OUT_DIR, METHOD_NAME, "_Split_Violin_PLA_Status.png"), width = 1400, height = 700)
p_split <- VlnPlot(pbmc, 
                   features = "Raw_Score", 
                   group.by = "celltype_clean", 
                   split.by = GT_COLUMN,   
                   pt.size = 0) +              
    scale_fill_manual(values = c("PLA" = "#FF4B4B", "platelet-free" = "#4B8BFF")) +
    labs(title = paste(METHOD_NAME, "Score Comparison: PLA vs. Platelet-free"),
         subtitle = "Internal score distribution within gated cell types",
         x = "Cell Type",
         y = "Raw Score",
         fill = "Gating Status") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
print(p_split)
dev.off()

# --- Density Plots: Score distribution by Error Class ---
# Z-Score Density
p_faceted_z <- ggplot(pbmc@meta.data, aes(x = Z_Score)) +
  geom_density(fill = "steelblue", alpha = 0.6) +
  facet_wrap(~ Error_Type, scales = "free_y") +
  geom_vline(xintercept = THRESHOLD_Z, linetype = "dashed", color = "red") +
  theme_classic() +
  labs(title = "Z-Score Density by Classification Error Type",
       subtitle = paste("Red Line = Automatic Threshold (Z =", round(THRESHOLD_Z, 2), ")"),
       x = "Z-Score", y = "Density")

ggsave(filename = paste0(OUT_DIR, METHOD_NAME, "_ZScore_Density_Errors.png"), 
       plot = p_faceted_z, width = 10, height = 6)

# --- ROC Curve ---
png(filename = paste0(OUT_DIR, METHOD_NAME, "_ROC_Curve.png"), width = 700, height = 700)
plot(roc_obj, col="#E41A1C", lwd=3, 
     main=paste(METHOD_NAME, "ROC Curve (AUC:", round(auc(roc_obj),3), ")"))
dev.off()

# --- PR Curve ---
eval_df <- data.frame(score = pbmc$Z_Score, gt = pbmc$GT_Response) %>% 
    arrange(desc(score)) %>%
    mutate(tp_cum = cumsum(gt), fp_cum = cumsum(1 - gt),
           precision_vec = tp_cum / (tp_cum + fp_cum), recall_vec = tp_cum / sum(gt))

# Calculate current metrics for annotation
tp <- sum(pbmc$Error_Type == "TP")
fp <- sum(pbmc$Error_Type == "FP")
fn <- sum(pbmc$Error_Type == "FN")

Prec_Val <- tp / (tp + fp)
Rec_Val <- tp / (tp + fn)
F1_Score <- 2 * Prec_Val * Rec_Val / (Prec_Val + Rec_Val)

png(filename = paste0(OUT_DIR, METHOD_NAME, "_PR_Curve.png"), width = 800, height = 700)
print(ggplot(eval_df, aes(x = recall_vec, y = precision_vec)) +
    geom_line(color = "#E41A1C", linewidth = 1.2) +
    annotate("point", x = Rec_Val, y = Prec_Val, color = "black", size = 4, shape = 18) +
    labs(title = paste("Precision-Recall Curve:", METHOD_NAME), 
         subtitle = paste("Best Youden Point - Precision:", round(Prec_Val, 2), "Recall:", round(Rec_Val, 2)),
         x = "Recall (Sensitivity)", y = "Precision") +
    theme_minimal())
dev.off()

# --- Clean Violin Plot ---
png(filename = paste0(OUT_DIR, METHOD_NAME, "_ZScore_Violin.png"), width = 1000, height = 700)
print(ggplot(pbmc@meta.data, aes(x = celltype_clean, y = Z_Score, fill = celltype_clean)) +
    geom_violin(alpha = 0.7) +
    geom_boxplot(width = 0.1, outlier.shape = NA) +
    geom_hline(yintercept = THRESHOLD_Z, linetype = "dashed", color = "red", linewidth = 1) +
    labs(title = paste(METHOD_NAME, "Z-Score Distribution across Cell Types"),
         subtitle = paste("Red Dashed Line = Optimal Threshold (Z =", round(THRESHOLD_Z, 2), ")"),
         x = "Cleaned Cell Types", y = "Z-Score") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none"))
dev.off()

print(paste("Done! All plots saved to:", OUT_DIR))

# --- SAVE PERFORMANCE METRICS FOR COMPARISON ---
performance_data <- data.frame(
    Method = METHOD_NAME,
    AUC = as.numeric(auc(roc_obj)),
    Precision = Prec_Val,
    Recall = Rec_Val,
    F1_Score = F1_Score,
    Threshold = THRESHOLD_Z,
    Timestamp = Sys.time()
)

write.csv(performance_data, 
          file = paste0("results/metrics_", METHOD_NAME, ".csv"), 
          row.names = FALSE)
