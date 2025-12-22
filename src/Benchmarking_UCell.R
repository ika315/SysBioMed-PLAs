# -------------------------------------------------------------------
# 04_Benchmarking_UCell.R
# -------------------------------------------------------------------

library(Seurat)
library(UCell)
library(pROC)
library(dplyr)
library(ggplot2)

METHOD <- "UCell"
PATH_DATA <- "~/SysBioMed-PLAs/data/seu_sx_final.rds"
pbmc <- readRDS(PATH_DATA)

# Ground truth & Signatur
gt_col_name <- "celltype.l2"
memory_b_genes <- c("RPS17", "MS4A1", "CD74", "HLA-DRA", "MIR1244-2", "RACK1", "IGKC", "RPL41", "LINC01857", "MT-ND3")
gene_set <- list(Memory_B_Score = memory_b_genes)

pbmc$GT_Response <- ifelse(pbmc[[gt_col_name]] == "B memory", 1, 0)
pbmc$GT_Class <- ifelse(pbmc$GT_Response == 1, "Positive", "Negative")

# --- SCORING ---
print(paste("--- Berechne", METHOD, "Score ---"))
pbmc <- AddModuleScore_UCell(
  pbmc, 
  features = gene_set, 
  assay = "RNA"
)

pbmc$Raw_Score <- pbmc$Memory_B_Score_UCell 
pbmc$Z_Score <- as.vector(scale(pbmc$Raw_Score))

# --- AUTOMATIC THRESHOLD (ROC) ---
roc_obj <- roc(response = pbmc$GT_Response, predictor = pbmc$Z_Score, direction = "<", quiet = TRUE)
best_coords <- coords(roc_obj, x = "best", best.method = "youden")
THRESHOLD_Z <- as.numeric(best_coords$threshold)

# Klassifizierung & Fehlerklassen
pbmc$Prediction <- ifelse(pbmc$Z_Score > THRESHOLD_Z, "Positive", "Negative")
ConfusionMatrix <- table(Predicted = pbmc$Prediction, Actual = pbmc$GT_Class)

TP <- ConfusionMatrix["Positive", "Positive"]; FP <- ConfusionMatrix["Positive", "Negative"]
FN <- ConfusionMatrix["Negative", "Positive"]; TN <- ConfusionMatrix["Negative", "Negative"]

Precision <- TP / (TP + FP); Recall <- TP / (TP + FN); Accuracy <- (TP + TN) / sum(ConfusionMatrix)

pbmc$Error_Type <- case_when(
    pbmc$Prediction == "Positive" & pbmc$GT_Class == "Positive" ~ "TP",
    pbmc$Prediction == "Positive" & pbmc$GT_Class == "Negative" ~ "FP",
    pbmc$Prediction == "Negative" & pbmc$GT_Class == "Positive" ~ "FN",
    pbmc$Prediction == "Negative" & pbmc$GT_Class == "Negative" ~ "TN"
)
pbmc$Error_Type <- factor(pbmc$Error_Type, levels = c("TP", "FP", "FN", "TN"))

# --- VISUALISIERUNGEN ---

# 1. UMAP Raw Score
png(filename = paste0("plots/", METHOD, "_Bmemory_UMAP_Raw.png"), width = 900, height = 700)
print(FeaturePlot(pbmc, features = "Raw_Score", reduction = "umap", raster = TRUE) +
    scale_colour_viridis_c(option = "magma") + labs(title = paste(METHOD, "Raw Score")))
dev.off()

# 2. Violin l2
png(filename = paste0("plots/", METHOD, "_Bmemory_Violin_l2.png"), width = 1000, height = 600)
print(VlnPlot(pbmc, features = "Raw_Score", group.by = "celltype.l2", pt.size = 0) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)))
dev.off()

# 3. Error UMAP
png(filename = paste0("plots/", METHOD, "_Bmemory_Error_UMAP.png"), width = 1000, height = 800)
print(DimPlot(pbmc, group.by = "Error_Type", reduction = "umap", raster = TRUE) +
    scale_color_manual(values = c("TP"="#228B22", "FP"="#FF4500", "FN"="#1E90FF", "TN"="#D3D3D3")) +
    labs(title = paste(METHOD, "Error Mapping"), subtitle = paste("Auto-Threshold Z =", round(THRESHOLD_Z, 2))))
dev.off()

# 4. Facetted Density (Z-Score)
png(filename = paste0("plots/", METHOD, "_ZScore_Facetted_Errors.png"), width = 10, height = 6, units="in", res=300)
print(ggplot(pbmc@meta.data, aes(x = Z_Score)) + geom_density(fill = "steelblue", alpha = 0.6) +
    facet_wrap(~ Error_Type, scales = "free_y") + geom_vline(xintercept = THRESHOLD_Z, linetype = "dashed", color = "red") +
    theme_classic() + labs(title = paste(METHOD, "Z-Score by Error Class")))
dev.off()

# 5. ROC Kurve
png(filename = paste0("plots/", METHOD, "_Bmemory_ROC.png"), width = 700, height = 700)
plot(roc_obj, col="#E41A1C", lwd=3, main=paste(METHOD, "ROC (AUC:", round(auc(roc_obj),3), ")"))
dev.off()

# 6. PR Kurve
eval_df <- data.frame(score = pbmc$Z_Score, gt = pbmc$GT_Response) %>% arrange(desc(score)) %>%
    mutate(tp_cum = cumsum(gt), fp_cum = cumsum(1 - gt),
           precision_vec = tp_cum / (tp_cum + fp_cum), recall_vec = tp_cum / sum(gt))
png(filename = paste0("plots/", METHOD, "_Bmemory_PR_Optimized.png"), width = 800, height = 700)
print(ggplot(eval_df, aes(x = recall_vec, y = precision_vec)) + geom_line(color = "#E41A1C", linewidth = 1.2) +
    annotate("point", x = Recall, y = Precision, color = "black", size = 4, shape = 18) +
    labs(title = paste(METHOD, "PR Curve"), x = "Recall", y = "Precision") + theme_minimal())
dev.off()

# 7. Z-Score Distribution (3 Groups)
pbmc$ZScore_Groups <- factor(case_when(
    pbmc$celltype.l2 == "B memory" ~ "B memory",
    pbmc$celltype.l1 == "B" ~ "Rest B-cells",
    TRUE ~ "Other"
), levels = c("B memory", "Rest B-cells", "Other"))
png(filename = paste0("plots/", METHOD, "_Bmemory_ZScore_3Groups.png"), width = 800, height = 600)
print(ggplot(pbmc@meta.data, aes(x = ZScore_Groups, y = Z_Score, fill = ZScore_Groups)) +
    geom_violin(alpha = 0.7) + geom_boxplot(width = 0.1, outlier.shape = NA) +
    geom_hline(yintercept = THRESHOLD_Z, linetype = "dashed", color = "red") + theme_minimal())
dev.off()

# --- Violin: Z-Score Verteilung über detaillierten Zelltypen ---
print(paste("--- Erstelle detaillierte Z-Score Verteilung für", METHOD, "---"))

png(filename = paste0("plots/", METHOD, "_Bmemory_ZScore_Distribution.png"), width = 1500, height = 800)

p_z_dist <- ggplot(pbmc@meta.data, aes(x = .data[[gt_col_name]], y = Z_Score, fill = .data[[gt_col_name]])) +
    geom_violin(alpha = 0.7) +
    geom_boxplot(width = 0.1, outlier.shape = NA) +
    geom_hline(yintercept = THRESHOLD_Z, linetype = "dashed", color = "red", linewidth = 1) +
    labs(title = paste(METHOD, ": Z-Score Verteilung über alle Zelltypen"),
         subtitle = paste("Roter Strich = Automatischer Threshold (Z =", round(THRESHOLD_Z, 2), ")"),
         x = "Zelltyp (Ground Truth)",
         y = "Z-Score") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none")

print(p_z_dist)
dev.off()



