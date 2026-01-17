# -------------------------------------------------------------------
# 04_Benchmarking_AddModuleScore.R
# -------------------------------------------------------------------

library(Seurat)
library(pROC)
library(dplyr)
library(ggplot2)

METHOD_NAME <- "AddModuleScore"
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

# GT and Singature Genes for Platelets
pbmc$GT_Response <- ifelse(pbmc[[GT_COLUMN]] == POSITIVE_VAL, 1, 0)
pbmc$GT_Class    <- ifelse(pbmc$GT_Response == 1, "Positive", "Negative")
base_dir <- getwd()
source(file.path(base_dir, "src","read_and_extend_gene_list.R"))
genes <- read_gene_list(file.path(base_dir, "data", "updated_gene_list.csv"))

# AddModuleScore Original für Erweiterung
print(paste("--- Berechne", METHOD_NAME, "(Original) für Gen-Erweiterung ---"))
pbmc <- AddModuleScore(pbmc, features = list(genes), name = "AMS_Orig")
pbmc$Raw_Score_Original <- pbmc$AMS_Orig1

# Gen-set Erweiterung
message("=== Erweitere Platelet-Gen-Set basierend auf AddModuleScore ===")
res_ext <- extend_gene_set(
  pbmc          = pbmc,
  base_genes    = genes,
  score_name    = "Raw_Score_Original", 
  top_n         = 50,
  high_quantile = 0.9,
  min.pct       = 0.05,
  test.use      = "wilcox"
)

extended_genes <- res_ext$extended_genes
gene_set_extended <- list(Platelet_Score = extended_genes)

#Speichern
if (!dir.exists("results")) dir.create("results")
write.csv(data.frame(geneName = extended_genes), 
          "results/extended_platelet_gene_list_AddModuleScore.csv", row.names = FALSE, quote = FALSE)

# Final Scoring
print(paste("--- Berechne finales", METHOD_NAME, "Scoring (Extended) ---"))
pbmc <- AddModuleScore(pbmc, features = list(extended_genes), name = "AMS_Final")
pbmc$Raw_Score <- pbmc$AMS_Final1
pbmc$Z_Score <- as.vector(scale(pbmc$Raw_Score))

# --- AUTOMATIC THRESHOLD (ROC) ---
roc_obj <- roc(response = pbmc$GT_Response, predictor = pbmc$Z_Score, direction = "<", quiet = TRUE)
best_coords <- coords(roc_obj, x = "best", best.method = "youden")
THRESHOLD_Z <- as.numeric(best_coords$threshold)

# Klassifizierung & Metriken
pbmc$Prediction <- ifelse(pbmc$Z_Score > THRESHOLD_Z, "Positive", "Negative")
ConfusionMatrix <- table(Predicted = pbmc$Prediction, Actual = pbmc$GT_Class)
TP <- ConfusionMatrix["Positive", "Positive"]
FP <- ConfusionMatrix["Positive", "Negative"]
FN <- ConfusionMatrix["Negative", "Positive"]
TN <- ConfusionMatrix["Negative", "Negative"]

Precision <- TP / (TP + FP)
Recall <- TP / (TP + FN)
Accuracy <- (TP + TN) / sum(ConfusionMatrix)

pbmc$Error_Type <- case_when(
    pbmc$Prediction == "Positive" & pbmc$GT_Class == "Positive" ~ "TP",
    pbmc$Prediction == "Positive" & pbmc$GT_Class == "Negative" ~ "FP",
    pbmc$Prediction == "Negative" & pbmc$GT_Class == "Positive" ~ "FN",
    pbmc$Prediction == "Negative" & pbmc$GT_Class == "Negative" ~ "TN"
)
pbmc$Error_Type <- factor(pbmc$Error_Type, levels = c("TP", "FP", "FN", "TN"))

# --- VISUALISIERUNGEN ---

# 1. UMAP Raw
png(filename = paste0(OUT_DIR,METHOD_NAME,"_",TARGET_LABEL, "UMAP_Scores.png"), width = 900, height = 700)
print(FeaturePlot(pbmc, features = "Raw_Score", reduction = "umap", raster = TRUE) + 
    scale_colour_viridis_c(option = "magma") + labs(title = paste(METHOD_NAME, "PLA-Signatur Score"),
         subtitle = "Intensität der Plättchen-Signatur in Leukozyten") +
    theme_minimal())
dev.off()

# 2. Violin Raw Scores
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


# 3. Error UMAP
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

# 6. ROC
png(filename = paste0(OUT_DIR, METHOD_NAME, "_", TARGET_LABEL, "_ROC.png"), width = 700, height = 700)
plot(roc_obj, col="#E41A1C", lwd=3, main=paste(METHOD_NAME, TARGET_LABEL, "ROC (AUC:", round(auc(roc_obj),3), ")"))
dev.off()

# 7. PR Kurve
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

# 8. Z-Score 3 Groups
#pbmc$ZScore_Groups <- factor(case_when(
#    pbmc$celltype.l2 == "Platelet" ~ "Platelet",
#    pbmc$celltype.l1 == "Mono" ~ "Monocytes", 
#    TRUE ~ "Other"
#), levels = c("Platelet", "Monocytes", "Other"))
#png(filename = paste0(OUT_DIR, METHOD, "_Platelets_ZScore_3Groups.png"), width = 800, height = 600)
#print(ggplot(pbmc@meta.data, aes(x = ZScore_Groups, y = Z_Score, fill = ZScore_Groups)) +
#    geom_violin(alpha = 0.7) + geom_boxplot(width = 0.1, outlier.shape = NA) +
#    geom_hline(yintercept = THRESHOLD_Z, linetype = "dashed", color = "red") +
#    scale_fill_manual(values = c("Platelet"="#FF4B4B", "Monocytes"="#4B8BFF", "Other"="#BEBEBE")) +
#    theme_minimal())
#dev.off()

# Violin Plot (Über die gegateten Klassen)
png(filename = paste0(OUT_DIR, METHOD_NAME, "_", TARGET_LABEL, "_Violin_GT.png"), width = 800, height = 600)
print(VlnPlot(pbmc, features = "Z_Score", group.by = GT_COLUMN, pt.size = 0) +
    geom_hline(yintercept = THRESHOLD_Z, linetype = "dashed", color = "red") +
    labs(title = paste("Z-Score Verteilung nach", TARGET_LABEL)))
dev.off()

# --- Violin: Z-Score Verteilung über detaillierten Zelltypen ---
print(paste("--- Erstelle detaillierte Z-Score Verteilung für", METHOD_NAME, "---"))
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
