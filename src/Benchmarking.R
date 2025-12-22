# -------------------------------------------------------------------
# 04_Benchmarking.R
# Benchmarking von PLA-Methoden anhand Memory B-Zell Signatur
# -------------------------------------------------------------------

library(Seurat)
library(AUCell)
library(UCell)
library(pROC)
library(dplyr)
library(ggplot2)

PATH_DATA <- "~/SysBioMed-PLAs/data/seu_sx_final.rds"
DATASET_NAME <- "seu_sx_final"

print(paste("Lade prozessiertes Seurat-Objekt:", PATH_DATA))
pbmc <- readRDS(PATH_DATA)

# Ground truth (GT) Definition
gt_col_name <- "celltype.l2"

# Memory B-Zell Signatur
memory_b_genes <- c("RPS17", "MS4A1", "CD74", "HLA-DRA", "MIR1244-2",
                    "RACK1", "IGKC", "RPL41", "LINC01857", "MT-ND3")
gene_set <- list(Memory_B_Score = memory_b_genes)

pbmc$GT_Response <- ifelse(pbmc[[gt_col_name]] == "B memory", 1, 0)

# --- AUCELL SCORING ---
print("--- Berechne AUCell Score ---")
expression_matrix <- GetAssayData(pbmc, layer = "data")
cells_rankings <- AUCell::AUCell_buildRankings(expression_matrix, plotStats=FALSE)
cells_AUC <- AUCell::AUCell_calcAUC(gene_set, cells_rankings)

pbmc$AUCell_Raw <- as.numeric(AUCell::getAUC(cells_AUC)[1, ])

# Z-Score berechnen
print("Berechne Z-Score für AUCell...")
pbmc$AUCell_ZScore <- as.vector(scale(pbmc$AUCell_Raw))


THRESHOLD_Z <- 1.5
pbmc$Prediction <- ifelse(pbmc$AUCell_ZScore > THRESHOLD_Z, "Positive", "Negative")
pbmc$GT_Class <- ifelse(pbmc$GT_Response == 1, "Positive", "Negative")

ConfusionMatrix <- table(Predicted = pbmc$Prediction, Actual = pbmc$GT_Class)
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

print("Confusion Matrix:")
print(ConfusionMatrix)
print(paste("Metriken: Prec:", round(Precision,4), "Rec:", round(Recall,4), "F1:", round(F1_Score,4), "Acc:", round(Accuracy,4)))

pbmc$Error_Type <- case_when(
    pbmc$Prediction == "Positive" & pbmc$GT_Class == "Positive" ~ "TP",
    pbmc$Prediction == "Positive" & pbmc$GT_Class == "Negative" ~ "FP",
    pbmc$Prediction == "Negative" & pbmc$GT_Class == "Positive" ~ "FN",
    pbmc$Prediction == "Negative" & pbmc$GT_Class == "Negative" ~ "TN"
)

# --- Visulaisierungen ---

# --- Refernz UMAP mit celltypes ---
print("--- Erstelle Referenz-UMAP der Zelltypen ---")

png(filename = "plots/04_Reference_UMAP_Celltypes.png", width = 1200, height = 800)

p_ref <- DimPlot(pbmc, 
                 group.by = "celltype.l1", 
                 reduction = "umap", 
                 label = TRUE, 
                 repel = TRUE, 
                 label.size = 5, 
                 raster = TRUE) +
    labs(title = "Referenz-UMAP: Zelltypen (l1)",
         subtitle = "Basiert auf der originalen Annotation (Ground Truth)",
         color = "Haupt-Zelltyp") +
    theme_minimal() +
    theme(legend.position = "right")

print(p_ref)
dev.off()


# --- UMAP: AUCell Raw Scores ---
png(filename = "plots/AUCell_Bmemory_UMAP_Raw.png", width = 900, height = 700)

p_umap_raw <- FeaturePlot(pbmc, features = "AUCell_Raw", reduction = "umap", 
                         label = TRUE, repel = TRUE, raster = TRUE) +
    scale_colour_viridis_c(option = "magma") +
    labs(title = "AUCell Score (Memory B) auf UMAP", color = "AUCell-Score") +
    theme_minimal()
dev.off()

# --- Violin: AUCell Raw Scores ---
png(filename = "plots/AUCell_Bmemory_Violin_l2.png", width = 1000, height = 600)
p_l1 <- VlnPlot(pbmc, features = "AUCell_Raw", group.by = "celltype.l2", pt.size = 0) +
    labs(title = "AUCell Score Verteilung über Haupt-Zelllinien (l2)") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
print(p_l1)
dev.off()

# --- UMAP: Error Mapping (TP, FP, FN, TN) ---
png(filename = "plots/AUCell_Bmemory_Error_UMAP.png", width = 1000, height = 800)
DimPlot(pbmc, group.by = "Error_Type", reduction = "umap", raster = TRUE) +
    scale_color_manual(values = c("TP"="#228B22", "FP"="#FF4500", "FN"="#1E90FF", "TN"="#D3D3D3")) +
    labs(title = "Mapping der Fehlerklassen auf UMAP", subtitle = paste("Threshold Z =", THRESHOLD_Z))
dev.off()

# --- Error Class Density Distribution ---
print("--- Erstelle Density Plots nach Fehlerklassen ---")
pbmc$Error_Type <- factor(pbmc$Error_Type, levels = c("TP", "FP", "FN", "TN"))

p_faceted_z <- ggplot(pbmc@meta.data, aes(x = AUCell_ZScore)) +
  geom_density(fill = "steelblue", alpha = 0.6) +
  facet_wrap(~ Error_Type, scales = "free_y") +
  geom_vline(xintercept = THRESHOLD_Z, linetype = "dashed", color = "red") +
  theme_classic() +
  labs(
    title = "AUCell Z-Score Verteilung nach Fehlerklassen",
    subtitle = paste("Rote Linie = Threshold Z =", THRESHOLD_Z),
    x = "AUCell Z-Score",
    y = "Dichte (Density)"
  )

ggsave(filename = "plots/AUCell_ZScore_Facetted_Errors.png", 
       plot = p_faceted_z, width = 10, height = 6)

p_faceted_raw <- ggplot(pbmc@meta.data, aes(x = AUCell_Raw)) +
  geom_density(fill = "orange", alpha = 0.6) +
  facet_wrap(~ Error_Type, scales = "free_y") +
  theme_classic() +
  labs(
    title = "AUCell Raw Score Verteilung nach Fehlerklassen",
    subtitle = "Einteilung basierend auf Z-Score Threshold",
    x = "AUCell Raw Score",
    y = "Dichte (Density)"
  )

ggsave(filename = "plots/AUCell_Raw_Facetted_Errors.png", 
       plot = p_faceted_raw, width = 10, height = 6)

# --- ROC Kurve ---
print("--- Erstelle ROC Plot ---")
roc_obj <- roc(response = pbmc$GT_Response, predictor = pbmc$AUCell_ZScore, direction = "<")
png(filename = "plots/AUCell_Bmemory_ROC.png", width = 700, height = 700)
plot(roc_obj, col="#E41A1C", lwd=3, main=paste("ROC Kurve (AUC:", round(auc(roc_obj),3), ")"))
dev.off()

# --- PR Kurve ---
print("--- Erstelle PR Plot ---")
eval_df <- data.frame(score = pbmc$AUCell_ZScore, gt = pbmc$GT_Response) %>% arrange(desc(score)) %>%
    mutate(tp_cum = cumsum(gt), fp_cum = cumsum(1 - gt),
           precision_vec = tp_cum / (tp_cum + fp_cum), recall_vec = tp_cum / sum(gt))

png(filename = "plots/AUCell_Bmemory_PR_Optimized.png", width = 800, height = 700)
ggplot(eval_df, aes(x = recall_vec, y = precision_vec)) +
    geom_line(color = "#E41A1C", linewidth = 1.2) +
    annotate("point", x = Recall, y = Precision, color = "black", size = 4, shape = 18) +
    annotate("text", x = Recall + 0.02, y = Precision + 0.02, 
             label = paste0("Aktueller Punkt\nPrec: ", round(Precision,2), "\nRec: ", round(Recall,2))) +
    labs(title = "Precision-Recall Kurve", x = "Recall", y = "Precision") +
    theme_minimal()
dev.off()

# --- Violin: Z-Score Verteilung (3 Gruppen) ---
pbmc$ZScore_Groups <- factor(case_when(
    pbmc$celltype.l2 == "B memory" ~ "B memory",
    pbmc$celltype.l1 == "B" ~ "Rest B-cells",
    TRUE ~ "Other"
), levels = c("B memory", "Rest B-cells", "Other"))

png(filename = "plots/AUCell_Bmemory_ZScore_3Groups.png", width = 800, height = 600)
ggplot(pbmc@meta.data, aes(x = ZScore_Groups, y = AUCell_ZScore, fill = ZScore_Groups)) +
    geom_violin(alpha = 0.7) + geom_boxplot(width = 0.1, outlier.shape = NA) +
    geom_hline(yintercept = THRESHOLD_Z, linetype = "dashed", color = "red") +
    scale_fill_manual(values = c("B memory"="#FF4B4B", "Rest B-cells"="#4B8BFF", "Other"="#BEBEBE")) +
    theme_minimal()
dev.off()


# --- Violin: Z-Score Verteilung über detaillierten Zelltypen ---
png(filename = "plots/AUCell_Bmemory_ZScore_Distribution.png", width = 1500, height = 800)
p_z <- ggplot(pbmc@meta.data, aes(x = .data[[gt_col_name]], y = AUCell_ZScore, fill = .data[[gt_col_name]])) +
    geom_violin(alpha = 0.7) +
    geom_boxplot(width = 0.1, outlier.shape = NA) +
    geom_hline(yintercept = THRESHOLD_Z, linetype = "dashed", color = "red", linewidth = 1) +
    labs(title = "AUCell Z-Score Verteilung über alle Zelltypen",
         subtitle = paste("Roter Strich = Aktueller Threshold (Z =", THRESHOLD_Z, ")"),
         x = "Zelltyp (Ground Truth)",
         y = "AUCell Z-Score") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none") # Legende hier aus, da X-Achse beschriftet ist
print(p_z)
dev.off()
