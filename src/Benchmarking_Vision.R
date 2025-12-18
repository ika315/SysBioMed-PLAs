# -------------------------------------------------------------------
# 04_Benchmarking.R
# Benchmarking von PLA-Methoden anhand Naive B-Zell Signatur
# Für Vision
# -------------------------------------------------------------------

library(Seurat)
library(VISION)
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
path_gene_list = "~/SysBioMed-PLAs/data/naive_b_cells.csv"
source("~/SysBioMed-PLAs/src/read_and_extend_gene_list.R")
genes <- read_gene_list(path_gene_list)
sig_vision <- make_signature(genes = genes, sig_name = "B_Cell_Signature", method = "vision")

pbmc$GT_Response <- ifelse(pbmc[[gt_col_name]] == "Naive B", 1, 0)


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
vision_score <- vision_auc[, "B_Cell_Signature"]
vision_score <- as.numeric(vision_auc[, "B_Cell_Signature"])

pbmc$Vision_B_Cell_Signature <- vision_score
score_name <- "Vision_B_Cell_Signature"

# UMAP
png(filename = "plots/Vision_NaiveB_UMAP_l2.png", width = 900, height = 700)
p_umap <- FeaturePlot(pbmc, features = "Vision_Raw", reduction = "umap", label = TRUE, label.size = 6, repel = TRUE, raster = TRUE) + 
  scale_colour_viridis_c(option = "magma") +
  labs(title = "Vision Score (Memory B) auf UMAP", 
       subtitle = "Farbskala: Gelb = Hoch, Dunkelblau = Niedrig | Legende zeigt l2-Klassen",
       color = "Vision Score") + 
  theme(legend.position = "right")	
print(p_umap)
dev.off()

# Violin Plot
png(filename = "plots/Vision_NaiveB_Violin_l2.png", width = 1000, height = 600)
p_l1 <- VlnPlot(pbmc, features = "Vision_Raw", group.by = "celltype.l2", pt.size = 0) +
  labs(title = "Vision Score Verteilung über Haupt-Zelllinien (l1)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
print(p_l1)
dev.off()

# Z-Score berechnen
print("Berechne Z-Score für Vision...")
pbmc$Vision_ZScore <- as.vector(scale(pbmc$Vision_Raw))

print("--- Führe AUC-ROC Benchmarking durch ---")
roc_obj <- pROC::roc(response = pbmc$GT_Response,
                     predictor = pbmc$Vision_ZScore,
                     levels = c(0, 1),
                     direction = "<")

print(paste("Vision ROC-AUC:", round(pROC::auc(roc_obj), 4)))

# Precision + Recall extrahieren und plotten
pr_coords <- coords(roc_obj, "all", ret = c("precision", "recall"), transpose = FALSE)

pr_coords <- pr_coords[complete.cases(pr_coords), ]
png(filename = "plots/Vision_NaiveB_PrecisionRecall_Curve.png", width = 800, height = 700)

p_pr <- ggplot(pr_coords, aes(x = recall, y = precision)) +
  geom_line(color = "#E41A1C", linewidth = 1.2) +
  geom_hline(yintercept = sum(pbmc$GT_Response == 1) / nrow(pbmc@meta.data), 
             linetype = "dashed", color = "grey") + 
  annotate("text", x = 0.5, y = 0.05, label = "Baseline", color = "grey") +
  labs(title = "Precision-Recall Kurve (B-Memory)",
       subtitle = paste("AUC-ROC:", round(pROC::auc(roc_obj), 4)),
       x = "Recall",
       y = "Precision") +
  ylim(0, 1) + xlim(0, 1) +
  theme_minimal()

print(p_pr)
dev.off()


png(filename = "plots/Vision_NaiveB_ZScore_Distribution.png", width = 1500, height = 800)

p_z <- ggplot(pbmc@meta.data, aes(x = .data[[gt_col_name]], y = Vision_ZScore, fill = .data[[gt_col_name]])) +
  geom_violin(alpha = 0.7) +
  geom_boxplot(width = 0.1, outlier.shape = NA) +
  geom_hline(yintercept = 1.5, linetype = "dashed", color = "red", linewidth = 1) +
  labs(title = "Vision Z-Score Verteilung über alle Zelltypen",
       subtitle = paste("Roter Strich = Aktueller Threshold (Z = 1.5)"),
       x = "Zelltyp (Ground Truth)",
       y = "Vision Z-Score") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")                          

print(p_z)
dev.off()

# --- FEHLERANALYSE (TP/FP/FN/TN) ---

# Hier Definition Threshold
THRESHOLD_Z <- 1.5 

pbmc$Prediction <- ifelse(pbmc$Vision_ZScore > THRESHOLD_Z, "Positive", "Negative")
pbmc$GT_Class <- ifelse(pbmc$GT_Response == 1, "Positive", "Negative")

# Confusion Matrix
ConfusionMatrix <- table(Predicted = pbmc$Prediction, Actual = pbmc$GT_Class)
print("Confusion Matrix (Positive = Naive B):")
print(ConfusionMatrix)

# Metriken extrahieren
TP <- ConfusionMatrix["Positive", "Positive"]
FP <- ConfusionMatrix["Positive", "Negative"] 
FN <- ConfusionMatrix["Negative", "Positive"]
TN <- ConfusionMatrix["Negative", "Negative"]

Precision <- TP / (TP + FP)
Accuracy <- (TP + TN) / sum(ConfusionMatrix)

print(paste("Metriken für B-Memory (Threshold Z =", THRESHOLD_Z, "):"))
print(paste("Precision:", round(Precision, 4)))
print(paste("Accuracy:", round(Accuracy, 4)))


# --- Visualisierung ---

# Violin Plot: Z-Score Verteilung B-Memory vs.B-general vs. Rest

pbmc$ZScore_Groups <- case_when(
  pbmc$celltype.l2 == "Naive B" ~ "Naive B",
  pbmc$celltype.l1 == "B" & pbmc$celltype.l2 != "Naive B" ~ "Rest B-cells",
  TRUE ~ "Other (Non-B)"
)

pbmc$ZScore_Groups <- factor(pbmc$ZScore_Groups, 
                             levels = c("Naive B", "Rest B-cells", "Other (Non-B)"))

png(filename = "plots/Vision_NaiveB_ZScore_3Groups.png", width = 800, height = 600)
p_z3 <- ggplot(pbmc@meta.data, aes(x = ZScore_Groups, y = Vision_ZScore, fill = ZScore_Groups)) +
  geom_violin(alpha = 0.7) +
  geom_boxplot(width = 0.1, outlier.shape = NA) +
  geom_hline(yintercept = THRESHOLD_Z, linetype = "dashed", color = "red", linewidth = 1) +
  scale_fill_manual(values = c("Naive B" = "#FF4B4B", "Rest B-cells" = "#4B8BFF", "Other (Non-B)" = "#BEBEBE")) +
  labs(title = "Spezifitäts-Check: Vision Z-Score",
       subtitle = paste("Rote Linie = Threshold Z =", THRESHOLD_Z),
       x = "Zellgruppen", y = "Vision Z-Score") +
  theme_minimal()
print(p_z3)
dev.off()


# Error Type Plot (TP, FP, FN, TN)
pbmc$Error_Type <- paste0(ifelse(pbmc$Prediction == "Positive", "P", "N"), 
                          ifelse(pbmc$GT_Class == "Positive", "T", "F"))

pbmc$Error_Type <- case_when(
  pbmc$Prediction == "Positive" & pbmc$GT_Class == "Positive" ~ "TP",
  pbmc$Prediction == "Positive" & pbmc$GT_Class == "Negative" ~ "FP",
  pbmc$Prediction == "Negative" & pbmc$GT_Class == "Positive" ~ "FN",
  pbmc$Prediction == "Negative" & pbmc$GT_Class == "Negative" ~ "TN"
)

png(filename = paste0("plots/04_Bench_ErrorTypes_NaiveB.png"), width = 900, height = 650)
p_error <- ggplot(pbmc@meta.data, aes(x = Error_Type, y = Vision_ZScore, fill = Error_Type)) +
  geom_boxplot() +
  geom_hline(yintercept = THRESHOLD_Z, linetype = "dashed", color = "red") +
  labs(title = "Vision Z-Score nach Fehlerklasse (Ziel: Naive B)") +
  theme_minimal()
print(p_error)
dev.off()