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

#Settings
METHOD_NAME <- "AUCell"
TARGET_LABEL <- "PLA_Gating"
GT_COLUMN <- "pla.status"
POSITIVE_VAL <- "PLA"

PATH_DATA <- "~/SysBioMed-PLAs/data/seu_sx_final.rds"
DATASET_NAME <- "seu_sx_final"
print(paste("Lade prozessiertes Seurat-Objekt:", PATH_DATA))
pbmc <- readRDS(PATH_DATA)

new_metadata <- read.csv("~/SysBioMed-PLAs/data/external_dataset_pla-status_metatable.csv", row.names = 1)
print("Barcodes in Seurat:")
print(head(Cells(pbmc), 3))
print("Barcodes in CSV:")
print(head(new_metadata$barcodes_clean, 3))

# Barcodes als Zeilennamen setzen
rownames(new_metadata) <- new_metadata$barcodes_clean

# Nur die Zellen behalten, die in beiden Objekten sind
common_cells <- intersect(Cells(pbmc), rownames(new_metadata))

print(paste("Zellen im Seurat-Objekt:", length(Cells(pbmc))))
print(paste("Zellen in der CSV-Tabelle:", length(rownames(new_metadata))))
print(paste("Gemeinsame Zellen:", length(common_cells)))

if(length(common_cells) == 0) {
  stop("FEHLER: Keine Übereinstimmung der Barcodes")
}
pbmc <- subset(pbmc, cells = common_cells)
new_metadata <- new_metadata[Cells(pbmc), ]
pbmc <- AddMetaData(pbmc, metadata = new_metadata)

print("Metadaten erfolgreich hinzugefügt!")
print(colnames(pbmc@meta.data))

OUT_DIR <- paste0("plots/Platelet_Main/", METHOD_NAME, "/")
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

# Memory B-Zell Signatur
#memory_b_genes <- c("RPS17", "MS4A1", "CD74", "HLA-DRA", "MIR1244-2",
#                    "RACK1", "IGKC", "RPL41", "LINC01857", "MT-ND3")
#gene_set <- list(Memory_B_Score = memory_b_genes)
#pbmc$GT_Response <- ifelse(pbmc[[gt_col_name]] == "B memory", 1, 0)

# GT and Singature Genes for Platelets
pbmc$GT_Response <- ifelse(pbmc[[GT_COLUMN]] == POSITIVE_VAL, 1, 0)
pbmc$GT_Class    <- ifelse(pbmc$GT_Response == 1, "Positive", "Negative")
base_dir <- getwd()
source(file.path(base_dir, "src","read_and_extend_gene_list.R"))
genes <- read_gene_list(file.path(base_dir, "data", "updated_gene_list.csv"))

print("--- Berechne AUCell (Original) für Gen-Erweiterung ---")
expression_matrix <- GetAssayData(pbmc, layer = "data")
rankings <- AUCell_buildRankings(expression_matrix, plotStats=FALSE)
orig_set <- list(Platelet_Orig = genes)
auc_orig <- AUCell_calcAUC(orig_set, rankings)
pbmc$AUCell_Raw_Original <- as.numeric(getAUC(auc_orig)[1, ])

# Gen-Set Erweiterung
message("=== Erweitere Platelet-Gen-Set basierend auf AUCell Scores ===")
res_ext <- extend_gene_set(
  pbmc          = pbmc,
  base_genes    = genes,
  score_name    = "AUCell_Raw_Original", 
  top_n         = 50,
  high_quantile = 0.9,
  min.pct       = 0.05,
  test.use      = "wilcox"
)
extended_genes <- res_ext$extended_genes
gene_set_extended <- list(Platelet_Score = extended_genes)

#Speichern
if (!dir.exists("results")) dir.create("results")
write.csv( 
  data.frame(geneName = extended_genes), 
  "results/extended_platelet_gene_list_AUCell.csv", 
  row.names = FALSE, 
  quote = FALSE 
)

print(paste("Originale Gene:", length(genes)))
print(paste("Erweiterte Gene:", length(extended_genes)))

# Finaler AUCell Score
print("--- Berechne finales AUCell Scoring (Extended) ---")
cells_AUC <- AUCell_calcAUC(gene_set_extended, rankings)
pbmc$Raw_Score <- as.numeric(getAUC(cells_AUC)[1, ])

# Z-Score berechnen
print("Berechne Z-Score für AUCell...")
pbmc$Z_Score <- as.vector(scale(pbmc$Raw_Score))

#Automatisierter Threshold
print("--- Berechne optimalen Threshold via ROC (Youden-Index) ---")
roc_obj <- pROC::roc(response = pbmc$GT_Response, 
                     predictor = pbmc$Z_Score, 
                     direction = "<", 
                     quiet = TRUE)

best_coords <- pROC::coords(roc_obj, x = "best", best.method = "youden")
THRESHOLD_Z <- as.numeric(best_coords$threshold)

print(paste("Automatischer Threshold für Platelets festgelegt auf Z =", round(THRESHOLD_Z, 4)))

pbmc$Prediction <- ifelse(pbmc$Z_Score > THRESHOLD_Z, "Positive", "Negative")

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

png(filename = paste0(OUT_DIR, "Reference_UMAP_Gated_Celltypes.png"), width = 1200, height = 800)
p_ref <- DimPlot(pbmc, group.by = "celltype_clean", reduction = "umap", 
                 label = TRUE, repel = TRUE, raster = TRUE) +
    labs(title = "Referenz-UMAP: Gated Celltypes",
         subtitle = "Basierend auf der neuen sauberen Annotation",
         color = "Zelltyp") +
    theme_minimal()
print(p_ref)
dev.off()

# --- UMAP: AUCell Raw Scores ---
png(filename = paste0(OUT_DIR, METHOD_NAME, "_", TARGET_LABEL, "_UMAP_Scores.png"), width = 900, height = 700)
p_umap_score <- FeaturePlot(pbmc, features = "Raw_Score", reduction = "umap", raster = TRUE) +
    scale_colour_viridis_c(option = "magma") +
    labs(title = paste(METHOD_NAME, "PLA-Signatur Score"), 
         subtitle = "Intensität der Plättchen-Signatur in Leukozyten") +
    theme_minimal()
print(p_umap_score)
dev.off()

# --- Violin: AUCell Raw Scores ---
png(filename = paste0(OUT_DIR, METHOD_NAME, "_", TARGET_LABEL, "_Violin_Celltypes.png"), width = 1200, height = 600)
p_vln_clean <- VlnPlot(pbmc, 
                       features = "Raw_Score", 
                       group.by = "celltype_clean", 
                       pt.size = 0) +
    labs(title = paste(METHOD_NAME, "Raw Score Verteilung"),
         subtitle = paste("Gruppiert nach: celltype_clean | Ground Truth:", TARGET_LABEL),
         x = "Zelltyp (Gated)",
         y = "Raw Score (AUC)") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none")

print(p_vln_clean)
dev.off()

plot_path_split <- paste0(OUT_DIR, METHOD_NAME, "_", TARGET_LABEL, "_Split_Violin_PLA_Status.png")
png(filename = plot_path_split, width = 1400, height = 700)

p_split <- VlnPlot(pbmc, 
                   features = "Raw_Score", 
                   group.by = "celltype_clean", 
                   split.by = GT_COLUMN,   
                   pt.size = 0) +             
    scale_fill_manual(values = c("PLA" = "#FF4B4B", "platelet-free" = "#4B8BFF")) +
    labs(title = paste(METHOD_NAME, "Score: PLA vs. Platelet-free"),
         subtitle = "Vergleich des Scores innerhalb der gegateten Zelltypen",
         x = "Zelltyp (celltype_clean)",
         y = "Raw Score (AUC)",
         fill = "Gating Status") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(p_split)
dev.off()

# --- UMAP: Error Mapping (TP, FP, FN, TN) ---
png(filename = paste0(OUT_DIR, METHOD_NAME, "_", TARGET_LABEL, "_Error_UMAP.png"), width = 1000, height = 800)
print(DimPlot(pbmc, group.by = "Error_Type", reduction = "umap") +
    scale_color_manual(values = c("TP"="#228B22", "FP"="#FF4500", "FN"="#1E90FF", "TN"="#D3D3D3")) +
    labs(title = paste(METHOD_NAME, "Error Mapping:", TARGET_LABEL), subtitle = paste("Threshold Z =", THRESHOLD_Z))+
    theme_minimal())
dev.off()

# --- Error Class Density Distribution ---
print("--- Erstelle Density Plots nach Fehlerklassen ---")
pbmc$Error_Type <- factor(pbmc$Error_Type, levels = c("TP", "FP", "FN", "TN"))

# Z-Score Density
p_faceted_z <- ggplot(pbmc@meta.data, aes(x = Z_Score)) +
  geom_density(fill = "steelblue", alpha = 0.6) +
  facet_wrap(~ Error_Type, scales = "free_y") +
  geom_vline(xintercept = THRESHOLD_Z, linetype = "dashed", color = "red") +
  theme_classic() +
  labs(title = "Z-Score Dichte nach PLA-Fehlerklassen",
       subtitle = paste("Rote Linie = Automatischer Threshold (Z =", round(THRESHOLD_Z, 2), ")"),
       x = "Z-Score", y = "Dichte")

ggsave(filename = paste0(OUT_DIR, METHOD_NAME, "_", TARGET_LABEL, "_ZScore_Density_Errors.png"),
       plot = p_faceted_z, width = 10, height = 6)

p_faceted_raw <- ggplot(pbmc@meta.data, aes(x = Raw_Score)) +
  geom_density(fill = "orange", alpha = 0.6) +
  facet_wrap(~ Error_Type, scales = "free_y") +
  theme_classic() +
  labs(title = "Raw Score Dichte nach PLA-Fehlerklassen",
       x = "Raw Score", y = "Dichte")

ggsave(filename = paste0(OUT_DIR, METHOD_NAME, "_", TARGET_LABEL, "_RawScore_Density_Errors.png"),
       plot = p_faceted_raw, width = 10, height = 6)

# --- ROC Kurve ---
print("--- Erstelle ROC Plot ---")
png(filename = paste0(OUT_DIR, METHOD_NAME, "_", TARGET_LABEL, "_ROC.png"), width = 700, height = 700)
plot(roc_obj, col="#E41A1C", lwd=3, main=paste(METHOD_NAME, TARGET_LABEL, "ROC (AUC:", round(auc(roc_obj),3), ")"))
dev.off()

# --- PR Kurve ---
print("--- Erstelle PR Plot ---")
eval_df <- data.frame(score = pbmc$Z_Score, gt = pbmc$GT_Response) %>% 
    arrange(desc(score)) %>%
    mutate(tp_cum = cumsum(gt), 
           fp_cum = cumsum(1 - gt),
           precision_vec = tp_cum / (tp_cum + fp_cum), 
           recall_vec = tp_cum / sum(gt))
png(filename = paste0(OUT_DIR, METHOD_NAME, "_", TARGET_LABEL, "_PR_Optimized.png"), width = 800, height = 700)
p_pr <- ggplot(eval_df, aes(x = recall_vec, y = precision_vec)) +
    geom_line(color = "#E41A1C", linewidth = 1.2) +
    annotate("point", x = Recall, y = Precision, color = "black", size = 4, shape = 18) +
    annotate("text", x = Recall + 0.02, y = Precision + 0.02,
             label = paste0("Bester Youden-Punkt\nPrec: ", round(Precision,2), "\nRec: ", round(Recall,2))) +
    labs(title = paste("PR-Kurve: Vorhersage PLA-Status"), 
         subtitle = paste("Methode:", METHOD_NAME, "| Ground Truth:", TARGET_LABEL),
         x = "Recall (Sensitivität)", y = "Precision (Genauigkeit)") +
    theme_minimal()
print(p_pr)
dev.off()

#Z-Score 3 Groups
#pbmc$ZScore_Groups <- factor(case_when(
#    pbmc$celltype.l2 == "Platelet" ~ "Platelet",
#    pbmc$celltype.l1 == "Mono" ~ "Monocytes", # Ein Vergleichsanker
#    TRUE ~ "Other"
#), levels = c("Platelet", "Monocytes", "Other"))

#png(filename = "plots/AUCell_Platelets_ZScore_3Groups.png", width = 800, height = 600)
#ggplot(pbmc@meta.data, aes(x = ZScore_Groups, y = AUCell_ZScore, fill = ZScore_Groups)) +
#    geom_violin(alpha = 0.7) + geom_boxplot(width = 0.1, outlier.shape = NA) +
#    geom_hline(yintercept = THRESHOLD_Z, linetype = "dashed", color = "red") +
#    scale_fill_manual(values = c("Platelet"="#FF4B4B", "Monocytes"="#4B8BFF", "Other"="#BEBEBE")) +
#    theme_minimal()
#dev.off()

# Violin Plot (Über die gegateten Klassen)
png(filename = paste0(OUT_DIR, METHOD_NAME, "_", TARGET_LABEL, "_Violin_GT.png"), width = 800, height = 600)
print(VlnPlot(pbmc, features = "Z_Score", group.by = GT_COLUMN, pt.size = 0.1) +
    geom_hline(yintercept = THRESHOLD_Z, linetype = "dashed", color = "red") +
    labs(title = paste("Z-Score Verteilung nach", TARGET_LABEL)))
dev.off()

# --- Violin: Z-Score Verteilung über detaillierten Zelltypen ---
png(filename = paste0(OUT_DIR, METHOD_NAME, "_", TARGET_LABEL, "_ZScore_Distribution.png"), width = 1500, height = 800)

p_z <- ggplot(pbmc@meta.data, aes(x = celltype_clean, y = Z_Score, fill = celltype_clean)) +
    geom_violin(alpha = 0.7) +
    geom_boxplot(width = 0.1, outlier.shape = NA) +
    geom_hline(yintercept = THRESHOLD_Z, linetype = "dashed", color = "red", linewidth = 1) +
    labs(title = paste(METHOD_NAME, "Z-Score Verteilung: PLA-Gating"),
         subtitle = paste("Roter Strich = Optimaler Threshold aus ROC (Z =", round(THRESHOLD_Z, 2), ")"),
         x = "Gereinigte Zelltypen (celltype_clean)",
         y = "Z-Score") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none")

print(p_z)
dev.off()
