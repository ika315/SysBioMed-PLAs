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

print("--- Führe AUC-ROC Benchmarking durch ---")
roc_obj <- pROC::roc(response = pbmc$GT_Response,
                     predictor = pbmc$AUCell_ZScore,
                     levels = c(0, 1),
                     direction = "<")

print(paste("AUCell ROC-AUC:", round(pROC::auc(roc_obj), 4)))

# --- FEHLERANALYSE (TP/FP/FN/TN) ---

# Hier Definition Threshold
THRESHOLD_Z <- 2.0 

pbmc$Prediction <- ifelse(pbmc$AUCell_ZScore > THRESHOLD_Z, "Positive", "Negative")
pbmc$GT_Class <- ifelse(pbmc$GT_Response == 1, "Positive", "Negative")

# Confusion Matrix
ConfusionMatrix <- table(Predicted = pbmc$Prediction, Actual = pbmc$GT_Class)
print("Confusion Matrix (Positive = B memory):")
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

# Violin Plot: Z-Score Verteilung B-Memory vs. Rest
png(filename = paste0("plots/04_Bench_Violin_Bmemory_ZScore.png"), width = 800, height = 600)
p_violin <- ggplot(pbmc@meta.data, aes(x = GT_Class, y = AUCell_ZScore, fill = GT_Class)) +
    geom_violin(alpha = 0.5) +
    geom_boxplot(width = 0.1, outlier.shape = NA) +
    geom_hline(yintercept = THRESHOLD_Z, linetype = "dashed", color = "red") +
    labs(title = "AUCell Z-Score: B memory vs. Alle anderen",
         x = "Ist die Zelle eine B memory? (Ground Truth)",
         y = "AUCell Z-Score") +
    theme_minimal()
print(p_violin)
dev.off()


# Error Type Plot (TP, FP, FN, TN)
pbmc$Error_Type <- paste0(ifelse(pbmc$Prediction == "Positive", "P", "N"), 
                          ifelse(pbmc$GT_Class == "Positive", "T", "F"))
# Korrektur der Benennung für Klarheit:
pbmc$Error_Type <- case_when(
    pbmc$Prediction == "Positive" & pbmc$GT_Class == "Positive" ~ "TP",
    pbmc$Prediction == "Positive" & pbmc$GT_Class == "Negative" ~ "FP",
    pbmc$Prediction == "Negative" & pbmc$GT_Class == "Positive" ~ "FN",
    pbmc$Prediction == "Negative" & pbmc$GT_Class == "Negative" ~ "TN"
)

png(filename = paste0("plots/04_Bench_ErrorTypes_Bmemory.png"), width = 900, height = 650)
p_error <- ggplot(pbmc@meta.data, aes(x = Error_Type, y = AUCell_ZScore, fill = Error_Type)) +
    geom_boxplot() +
    geom_hline(yintercept = THRESHOLD_Z, linetype = "dashed", color = "red") +
    labs(title = "AUCell Z-Score nach Fehlerklasse (Ziel: B memory)") +
    theme_minimal()
print(p_error)
dev.off()

