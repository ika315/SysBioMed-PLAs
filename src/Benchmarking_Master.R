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
library(pheatmap)
library(mclust)
# Configuration
# Options: "AUCell", "UCell", "AddModuleScore"
#METHOD_NAME  <- "AUCell"

args <- commandArgs(trailingOnly = TRUE)

if (length(args) >= 1) {
  METHOD_NAME <- args[1] 
} else {
  METHOD_NAME <- "AUCell" 
}

SIG_NAME     <- "MANNE_COVID19_DN"
SIGNATURE    <- "MANNE_COVID19_COMBINED_COHORT_VS_HEALTHY_DONOR_PLATELETS_DN"
TARGET_LABEL <- "PLA_Gating"
GT_COLUMN    <- "pla.status"
POSITIVE_VAL <- "PLA"
PATH_DATA    <- "~/SysBioMed-PLAs/data/seu_sx_integration.rds"
IMMUNE_SIG <- "GOBP_LEUKOCYTE_ACTIVATION_INVOLVED_IN_INFLAMMATORY_RESPONSE.v2025.1.Hs"

# Load data and Meta data Sync
print(paste("Running Benchmark for Method:", METHOD_NAME, " with Gene Signature:", SIG_NAME))
pbmc <- readRDS(PATH_DATA)

new_metadata <- read.csv("~/SysBioMed-PLAs/data/external_dataset_pla-status_metatable.csv", row.names = 1)
rownames(new_metadata) <- new_metadata$barcodes_clean
common_cells <- intersect(Cells(pbmc), rownames(new_metadata))
pbmc <- subset(pbmc, cells = common_cells)
pbmc <- AddMetaData(pbmc, metadata = new_metadata[common_cells, ])

# Define Output Directory
OUT_DIR <- paste0("plots/Platelet_Main_New/",SIG_NAME,"/", METHOD_NAME, "/")
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)
if (!dir.exists("results")) dir.create("results")

# Load Gene Lists
base_dir <- getwd()
source(file.path(base_dir, "src", "read_and_extend_gene_list.R"))
#genes <- read_gene_list(file.path(base_dir, "data", "updated_gene_list.csv"))
PATH_SIG <- file.path(base_dir, "data", paste0(SIGNATURE, ".v2025.1.Hs.csv"))
genes <- read_gene_list(PATH_SIG)

pbmc$GT_Response <- ifelse(pbmc[[GT_COLUMN]] == POSITIVE_VAL, 1, 0)
pbmc$GT_Class    <- ifelse(pbmc$GT_Response == 1, "Positive", "Negative")

# Load Immune Marker Genes
PATH_IMMUNE_SIG <- file.path(base_dir, "data", paste0(IMMUNE_SIG, ".csv"))
immune_genes <- read_gene_list(PATH_IMMUNE_SIG)
immune_genes <- intersect(immune_genes, rownames(GetAssayData(pbmc, layer = "data")))

 
# Grid Search
#source(file.path(base_dir, "src", "Grid_Search.R"))
#cell_rankings <- AUCell_buildRankings(GetAssayData(pbmc, layer = "data"), plotStats=FALSE)
# grid_search(pbmc, genes, cell_rankings)

# Scoring Logic
print(paste("--- Calculating Scores using", METHOD_NAME, "---"))

if (METHOD_NAME == "AUCell" || METHOD_NAME == "WeightedAUCell") {
    expression_matrix <- GetAssayData(pbmc, layer = "data")
    rankings <- AUCell_buildRankings(expression_matrix, plotStats=FALSE)
    auc_orig <- AUCell_calcAUC(list(Platelet_Orig = genes), rankings)
    pbmc$Raw_Score_Original <- as.numeric(getAUC(auc_orig)[1, ])

    auc_imm <- AUCell_calcAUC(list(Immune_Score= immune_genes), rankings)
    pbmc$Immune_Score <- as.numeric(getAUC(auc_imm)[1, ])
} else if (METHOD_NAME == "UCell") {
    pbmc <- AddModuleScore_UCell(pbmc, features = list(Platelet_Orig = genes), name = NULL)
    pbmc$Raw_Score_Original <- pbmc$Platelet_Orig
} else if (METHOD_NAME == "AddModuleScore") {
    pbmc <- AddModuleScore(pbmc, features = list(genes), name = "AMS_Orig")
    pbmc$Raw_Score_Original <- pbmc$AMS_Orig1
}

# Geneset Extension mit Skip-Logik
EXT_FILE <- paste0("results/extended_platelet_gene_list_", SIG_NAME, "_", METHOD_NAME, ".csv")

if (file.exists(EXT_FILE)) {
    message("Lade existierende erweiterte Genliste: ", EXT_FILE)
    ext_df <- read.csv(EXT_FILE)
    extended_genes <- ext_df$geneName
} else {
    message("Berechne neue Geneset Extension (das kann dauern)...")
    res_ext <- extend_gene_set(
        pbmc = pbmc, base_genes = genes, score_name = "Raw_Score_Original", 
        top_n = 50, high_quantile = 0.9, min.pct = 0.05, test.use = "wilcox"
    )
    extended_genes <- res_ext$extended_genes
    write.csv(data.frame(geneName = extended_genes), EXT_FILE, row.names = FALSE)
}


#res_ext <- extend_gene_set(
#    pbmc = pbmc, base_genes = genes, score_name = "Raw_Score_Original", 
#    top_n = 50, high_quantile = 0.9, min.pct = 0.05, test.use = "wilcox"
#)
#extended_genes <- res_ext$extended_genes
#write.csv(data.frame(geneName = extended_genes), 
#          paste0("results/extended_platelet_gene_list_", SIG_NAME, "_", METHOD_NAME, ".csv"), row.names = FALSE)

# Final Scoring
if(METHOD_NAME == "WeightedAUCell") {
    expression_matrix <- GetAssayData(pbmc, layer = "data")

    genes_in_data <- intersect(genes, rownames(expression_matrix))
    expr_subset <- expression_matrix[genes_in_data, , drop = FALSE]

    genes_auc <- sapply(rownames(expr_subset), function(g) {
        a <- as.numeric(pROC::roc(pbmc$GT_Response, expr_subset[g, ], quiet = TRUE)$auc)
        max(a, 1 - a)
    })
    names(genes_auc) <- rownames(expr_subset)

    boost_full <- setNames(rep(1, nrow(expression_matrix)), rownames(expression_matrix))
    known_auc <- genes_auc[intersect(names(genes_auc), names(boost_full))]
    boost_full[names(known_auc)] <- pmax(0.1, as.numeric(known_auc)^3)

    weight_matrix <- sweep(expression_matrix, 1, boost_full[rownames(expression_matrix)], `*`)

    rankings_weighted <- AUCell_buildRankings(weight_matrix, plotStats = FALSE)

    final_genes <- unique(c(genes_in_data, extended_genes))
    final_genes <- intersect(final_genes, rownames(rankings_weighted))

    auc_final <- AUCell_calcAUC(
        list(Platelet_Score = final_genes),
        rankings_weighted,
        aucMaxRank = ceiling(0.05 * nrow(rankings_weighted))
    )

    pbmc$Raw_Score <- as.numeric(getAUC(auc_final)[1, ])
}else if (METHOD_NAME == "AUCell") {
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
roc_obj <- roc(response = pbmc$GT_Response, predictor = pbmc$Z_Score, direction = "<", quiet = TRUE)

set.seed(42)

# observed AUC
auc_obs <- as.numeric(getAUC(auc_final)[1, ])
auc_imm <- as.numeric(getAUC(auc_imm)[1, ])

x <- auc_obs 

fit_plat <- Mclust(auc_obs, G = 2)
plat_high <- which.max(fit_plat$parameters$mean)
pbmc$Platelet_High <- fit_plat$classification == plat_high

idx <- which(pbmc$Platelet_High)

fit_imm <- Mclust(auc_imm[idx], G = 2)
imm_high <- which.max(fit_imm$parameters$mean)

pbmc$Immune_High <- FALSE
pbmc$Immune_High[idx] <- fit_imm$classification == imm_high

pbmc$Prediction <- factor(
  ifelse((pbmc$Platelet_High) & (pbmc$Immune_High), "Positive", "Negative"),
  levels = c("Negative","Positive")
)

table(pbmc$Prediction)
# Metrics and Threshold
#roc_obj <- roc(response = pbmc$GT_Response, predictor = pbmc$Z_Score, direction = "<", quiet = TRUE)
#THRESHOLD_Z <- as.numeric(coords(roc_obj, x = "best", best.method = "youden")$threshold)
#pbmc$Prediction <- ifelse(pbmc$Z_Score > THRESHOLD_Z, "Positive", "Negative")

pbmc$Error_Type <- case_when(
    pbmc$Prediction == "Positive" & pbmc$GT_Class == "Positive" ~ "TP",
    pbmc$Prediction == "Positive" & pbmc$GT_Class == "Negative" ~ "FP",
    pbmc$Prediction == "Negative" & pbmc$GT_Class == "Positive" ~ "FN",
    pbmc$Prediction == "Negative" & pbmc$GT_Class == "Negative" ~ "TN"
)
pbmc$Error_Type <- factor(pbmc$Error_Type, levels = c("TP", "FP", "FN", "TN"))

tp <- sum(pbmc$Error_Type == "TP", na.rm = TRUE)
fp <- sum(pbmc$Error_Type == "FP", na.rm = TRUE)
fn <- sum(pbmc$Error_Type == "FN", na.rm = TRUE)

Prec_Val <- if((tp + fp) > 0) tp / (tp + fp) else 0
Rec_Val  <- if((tp + fn) > 0) tp / (tp + fn) else 0
F1_Score <- if((Prec_Val + Rec_Val) > 0) 2 * Prec_Val * Rec_Val / (Prec_Val + Rec_Val) else 0

performance_data <- data.frame(
    Method = METHOD_NAME,
    Signature = SIG_NAME,
    AUC = as.numeric(auc(roc_obj)),
    Precision = Prec_Val,
    Recall = Rec_Val,
    F1_Score = F1_Score,
    Threshold = THRESHOLD_Z,
    Timestamp = Sys.time()
)

write.csv(performance_data, 
          file = paste0("results/metrics_", SIG_NAME, "_", METHOD_NAME, ".csv"),
          row.names = FALSE)

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
         subtitle = paste("Signature:", SIG_NAME),
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

png(filename = paste0(OUT_DIR, METHOD_NAME, "_PR_Curve.png"), width = 800, height = 700)
print(ggplot(eval_df, aes(x = recall_vec, y = precision_vec)) +
    geom_line(color = "#E41A1C", linewidth = 1.2) +
    annotate("point", x = Rec_Val, y = Prec_Val, color = "black", size = 4, shape = 18) +
    labs(title = paste("Precision-Recall Curve:", METHOD_NAME), 
         subtitle = paste("Precision:", round(Prec_Val, 2), "Recall:", round(Rec_Val, 2)),
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

# --- Heatmap: TP vs. FP Signature Gene Expression ---

genes_to_show <- head(intersect(extended_genes, rownames(pbmc)), 30)
cells_subset <- subset(pbmc, subset = Error_Type %in% c("TP", "FP"))

if(length(genes_to_show) > 0 && ncol(cells_subset) > 0) {
    avg_exp <- AverageExpression(cells_subset, 
                                 features = genes_to_show, 
                                 group.by = "Error_Type", 
                                 layer = "data")$RNA

    png(filename = paste0(OUT_DIR, SIG_NAME, "_", METHOD_NAME, "_Heatmap_TP_vs_FP.png"), width = 800, height = 1000)
    pheatmap(avg_exp, 
             scale = "row", 
             clustering_distance_cols = "euclidean",
             main = paste("Gene Profile: TP vs FP (", SIG_NAME, ")"),
             color = colorRampPalette(c("blue", "white", "red"))(100))
    dev.off()
}

# --- Plot: FP Count per Cell Type ---

fp_data <- pbmc@meta.data %>%
    filter(Error_Type == "FP") %>%
    group_by(celltype_clean) %>%
    tally() %>%
    arrange(desc(n))

if(nrow(fp_data) > 0) {
    png(filename = paste0(OUT_DIR, SIG_NAME, "_", METHOD_NAME, "_FP_Count_per_CellType.png"), width = 1000, height = 700)
    p_fp_bar <- ggplot(fp_data, aes(x = reorder(celltype_clean, -n), y = n, fill = celltype_clean)) +
        geom_bar(stat = "identity") +
        geom_text(aes(label = n), vjust = -0.5, size = 5) +
        labs(title = paste("False Positives per Cell Type"),
             subtitle = paste("Method:", METHOD_NAME, "| Signature:", SIG_NAME, "| Total FPs:", sum(fp_data$n)),
             x = "Cell Type", y = "Number of FP Cells") +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")
    print(p_fp_bar)
    dev.off()
} else {
    message("No False Positives found to plot.")
}

# --- Heatmap: All Error Types (TP, FP, FN, TN) ---
cells_subset <- subset(pbmc, subset = Error_Type %in% c("TP", "FP", "FN", "TN"))

if(length(genes_to_show) > 0 && ncol(cells_subset) > 0) {
    avg_exp <- AverageExpression(cells_subset, 
                                 features = genes_to_show, 
                                 group.by = "Error_Type", 
                                 layer = "data")$RNA
    
    avg_exp <- avg_exp[, c("TP", "FP", "FN", "TN"), drop = FALSE]

    png(filename = paste0(OUT_DIR, SIG_NAME, "_", METHOD_NAME, "_Heatmap_All_Errors.png"), width = 1000, height = 1000)
    pheatmap(avg_exp, 
             scale = "row", 
             cluster_cols = FALSE, # Reihenfolge TP, FP, FN, TN beibehalten
             main = paste("Gene Profile: Comparison of All Error Types (", SIG_NAME, ")"),
             color = colorRampPalette(c("blue", "white", "red"))(100))
    dev.off()
}

# --- DotPlot Zoom: TP vs FP across Celltypes ---
cells_to_compare <- subset(pbmc, subset = Error_Type %in% c("TP", "FP"))

cells_to_compare$Plot_Group <- paste(cells_to_compare$celltype_clean, cells_to_compare$Error_Type, sep = "_")

png(filename = paste0(OUT_DIR, SIG_NAME, "_", METHOD_NAME, "_DotPlot_TP_vs_FP.png"), width = 1200, height = 800)
p_dot <- DotPlot(cells_to_compare, 
                 features = genes_to_show, 
                 group.by = "Plot_Group") + 
    coord_flip() + # Gene auf der Y-Achse sind meist besser lesbar
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_colour_gradient2(low = "blue", mid = "white", high = "red") +
    labs(title = "Zoom into FP: Gene Expression Comparison",
         subtitle = "TP vs FP profiles per Cell Type",
         x = "Signature Genes", y = "Cell Type & Error Class")

print(p_dot)
dev.off()

# Barplot with percentage of error types for each celltype 
error_counts <- pbmc@meta.data %>%
    group_by(celltype_clean, Error_Type) %>%
    tally() %>%
    group_by(celltype_clean) %>%
    mutate(perc = n / sum(n) * 100)

png(filename = paste0(OUT_DIR, SIG_NAME, "_", METHOD_NAME, "_Error_Composition_by_Celltype.png"), width = 1200, height = 700)
p_comp <- ggplot(error_counts, aes(x = celltype_clean, y = perc, fill = Error_Type)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = c("TP"="#228B22", "FP"="#FF4500", "FN"="#1E90FF", "TN"="#D3D3D3")) +
    labs(title = "Error Composition per Cell Type",
         y = "Percentage (%)", x = "Cell Type") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
print(p_comp)
dev.off()

expr <- GetAssayData(pbmc, layer = "data")
det_rate <- Matrix::rowMeans(expr > 0)

df_det <- data.frame(
  det_rate = det_rate,
  group = "Universe"
)
df_det_sig <- data.frame(
  det_rate = det_rate[intersect(names(det_rate), extended_genes)],
  group = "Signature"
)
df <- rbind(df_det, df_det_sig)

p <- ggplot(df, aes(x = det_rate, fill = group)) +
  geom_density(alpha = 0.4) +
  theme_classic() +
  labs(title = "Gene detection rate: signature vs universe",
       x = "Detection rate per gene (% cells > 0)", y = "Density")

ggsave(filename = paste0(OUT_DIR, SIG_NAME, "_", METHOD_NAME, "_Gene_Detection_Rate_Comparison.png"), 
       plot = p, width = 8, height = 6)


mean_expr <- Matrix::rowMeans(expr)

df_mean <- data.frame(mean_expr = mean_expr, group = "Universe")
df_mean_sig <- data.frame(mean_expr = mean_expr[intersect(names(mean_expr), extended_genes)],
                          group = "Signature")
df2 <- rbind(df_mean, df_mean_sig)

p <- ggplot(df2, aes(x = mean_expr, fill = group)) +
  geom_density(alpha = 0.4) +
  theme_classic() +
  labs(title = "Mean expression: signature vs universe",
       x = "Mean expression per gene", y = "Density")

ggsave(filename = paste0(OUT_DIR, SIG_NAME, "_", METHOD_NAME, "_Gene_Mean_Expression_Comparison.png"), 
       plot = p, width = 8, height = 6)

df_scatter <- data.frame(
  gene = names(det_rate),
  det_rate = as.numeric(det_rate),
  mean_expr = as.numeric(mean_expr),
  group = "Universe"
)
df_scatter$group[df_scatter$gene %in% extended_genes] <- "Signature"

p <- ggplot(df_scatter, aes(x = det_rate, y = mean_expr, color = group)) +
  geom_point(alpha = 0.3, size = 0.6) +
  theme_classic() +
  labs(title = "Gene characteristics",
       x = "Detection rate", y = "Mean expression")

ggsave(filename = paste0(OUT_DIR, SIG_NAME, "_", METHOD_NAME, "_Gene_Detection_vs_Mean_Expression.png"), 
       plot = p, width = 8, height = 6)


df2 <- data.frame(AUC = auc_obs, cluster = pbmc$performance_data <- data.frame(
    Method = METHOD_NAME,
    Signature = SIG_NAME,
    AUC = as.numeric(auc(roc_obj)),
    Precision = Prec_Val,
    Recall = Rec_Val,
    F1_Score = F1_Score,
    Threshold = THRESHOLD_Z,
    Timestamp = Sys.time()
)

write.csv(performance_data, 
          file = paste0("results/metrics_", SIG_NAME, "_", METHOD_NAME, ".csv"),
          row.names = FALSE)

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
         subtitle = paste("Signature:", SIG_NAME),
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

png(filename = paste0(OUT_DIR, METHOD_NAME, "_PR_Curve.png"), width = 800, height = 700)
print(ggplot(eval_df, aes(x = recall_vec, y = precision_vec)) +
    geom_line(color = "#E41A1C", linewidth = 1.2) +
    annotate("point", x = Rec_Val, y = Prec_Val, color = "black", size = 4, shape = 18) +
    labs(title = paste("Precision-Recall Curve:", METHOD_NAME), 
         subtitle = paste("Precision:", round(Prec_Val, 2), "Recall:", round(Rec_Val, 2)),
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

# --- Heatmap: TP vs. FP Signature Gene Expression ---

genes_to_show <- head(intersect(extended_genes, rownames(pbmc)), 30)
cells_subset <- subset(pbmc, subset = Error_Type %in% c("TP", "FP"))

if(length(genes_to_show) > 0 && ncol(cells_subset) > 0) {
    avg_exp <- AverageExpression(cells_subset, 
                                 features = genes_to_show, 
                                 group.by = "Error_Type", 
                                 layer = "data")$RNA

    png(filename = paste0(OUT_DIR, SIG_NAME, "_", METHOD_NAME, "_Heatmap_TP_vs_FP.png"), width = 800, height = 1000)
    pheatmap(avg_exp, 
             scale = "row", 
             clustering_distance_cols = "euclidean",
             main = paste("Gene Profile: TP vs FP (", SIG_NAME, ")"),
             color = colorRampPalette(c("blue", "white", "red"))(100))
    dev.off()
}

# --- Plot: FP Count per Cell Type ---

fp_data <- pbmc@meta.data %>%
    filter(Error_Type == "FP") %>%
    group_by(celltype_clean) %>%
    tally() %>%
    arrange(desc(n))

if(nrow(fp_data) > 0) {
    png(filename = paste0(OUT_DIR, SIG_NAME, "_", METHOD_NAME, "_FP_Count_per_CellType.png"), width = 1000, height = 700)
    p_fp_bar <- ggplot(fp_data, aes(x = reorder(celltype_clean, -n), y = n, fill = celltype_clean)) +
        geom_bar(stat = "identity") +
        geom_text(aes(label = n), vjust = -0.5, size = 5) +
        labs(title = paste("False Positives per Cell Type"),
             subtitle = paste("Method:", METHOD_NAME, "| Signature:", SIG_NAME, "| Total FPs:", sum(fp_data$n)),
             x = "Cell Type", y = "Number of FP Cells") +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")
    print(p_fp_bar)
    dev.off()
} else {
    message("No False Positives found to plot.")
}

# --- Heatmap: All Error Types (TP, FP, FN, TN) ---
cells_subset <- subset(pbmc, subset = Error_Type %in% c("TP", "FP", "FN", "TN"))

if(length(genes_to_show) > 0 && ncol(cells_subset) > 0) {
    avg_exp <- AverageExpression(cells_subset, 
                                 features = genes_to_show, 
                                 group.by = "Error_Type", 
                                 layer = "data")$RNA
    
    avg_exp <- avg_exp[, c("TP", "FP", "FN", "TN"), drop = FALSE]

    png(filename = paste0(OUT_DIR, SIG_NAME, "_", METHOD_NAME, "_Heatmap_All_Errors.png"), width = 1000, height = 1000)
    pheatmap(avg_exp, 
             scale = "row", 
             cluster_cols = FALSE, # Reihenfolge TP, FP, FN, TN beibehalten
             main = paste("Gene Profile: Comparison of All Error Types (", SIG_NAME, ")"),
             color = colorRampPalette(c("blue", "white", "red"))(100))
    dev.off()
}

# --- DotPlot Zoom: TP vs FP across Celltypes ---
cells_to_compare <- subset(pbmc, subset = Error_Type %in% c("TP", "FP"))

cells_to_compare$Plot_Group <- paste(cells_to_compare$celltype_clean, cells_to_compare$Error_Type, sep = "_")

png(filename = paste0(OUT_DIR, SIG_NAME, "_", METHOD_NAME, "_DotPlot_TP_vs_FP.png"), width = 1200, height = 800)
p_dot <- DotPlot(cells_to_compare, 
                 features = genes_to_show, 
                 group.by = "Plot_Group") + 
    coord_flip() + # Gene auf der Y-Achse sind meist besser lesbar
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_colour_gradient2(low = "blue", mid = "white", high = "red") +
    labs(title = "Zoom into FP: Gene Expression Comparison",
         subtitle = "TP vs FP profiles per Cell Type",
         x = "Signature Genes", y = "Cell Type & Error Class")

print(p_dot)
dev.off()

# Barplot with percentage of error types for each celltype 
error_counts <- pbmc@meta.data %>%
    group_by(celltype_clean, Error_Type) %>%
    tally() %>%
    group_by(celltype_clean) %>%
    mutate(perc = n / sum(n) * 100)

png(filename = paste0(OUT_DIR, SIG_NAME, "_", METHOD_NAME, "_Error_Composition_by_Celltype.png"), width = 1200, height = 700)
p_comp <- ggplot(error_counts, aes(x = celltype_clean, y = perc, fill = Error_Type)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = c("TP"="#228B22", "FP"="#FF4500", "FN"="#1E90FF", "TN"="#D3D3D3")) +
    labs(title = "Error Composition per Cell Type",
         y = "Percentage (%)", x = "Cell Type") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
print(p_comp)
dev.off()seurat_clusters)

ggplot(df2, aes(x = AUC)) +
  geom_histogram(bins = 80) +
  facet_wrap(~ cluster, scales = "free_y") +
  theme_classic() +
  labs(title = "AUC distribution per cluster", x = "AUCell AUC", y = "Count")

df <- data.frame(AUC = auc_obs)

ggplot(df, aes(x = AUC)) +
  geom_histogram(aes(y = after_stat(density)), bins = 80, alpha = 0.4) +
  geom_density(adjust = 1) +
  theme_classic() +
  labs(title = "AUC distribution (histogram + density)", x = "AUCell AUC", y = "Density")


df2 <- data.frame(PlateletAUC = auc_obs, ImmuneAUC = auc_imm, Class = pbmc$GT_Class)

p <- ggplot(df2, aes(x = PlateletAUC, y = ImmuneAUC, color = Class)) +
  geom_point(alpha = 0.4, size = 0.6) +
  theme_classic() +
  labs(title = "Platelet vs Immune score separates Platelets from PLAs",
       x = "Platelet AUCell AUC", y = "Immune AUCell AUC")
ggsave(filename = paste0(OUT_DIR, SIG_NAME, "_", METHOD_NAME, "_Platelet_vs_Immune_AUC_Scatter.png"), 
       plot = p, width = 8, height = 6)



library(ggplot2)
library(tidyr)
library(dplyr)

df <- data.frame(
  Platelet_AUC = pbmc$Raw_Score,
  Immune_AUC   = pbmc$Immune_Score
) %>%
  pivot_longer(cols = c(Platelet_AUC, Immune_AUC),
               names_to = "Score", values_to = "Value")

p <- ggplot(df, aes(x = Value)) +
  geom_density(fill = "grey70", alpha = 0.5) +
  facet_wrap(~ Score, scales = "free") +
  theme_classic() +
  labs(title = "Score distributions (density)",
       x = "Score value", y = "Density")

ggsave(filename = paste0(OUT_DIR, SIG_NAME, "_", METHOD_NAME, "_Score_Distributions_Density.png"),
         plot = p, width = 8, height = 6)


p <- ggplot(df, aes(x = Value)) +
  geom_histogram(aes(y = after_stat(density)), bins = 80, alpha = 0.35) +
  geom_density(size = 1) +
  facet_wrap(~ Score, scales = "free") +
  theme_classic() +
  labs(title = "Score distributions (histogram + density)",
       x = "Score value", y = "Density")
ggsave(filename = paste0(OUT_DIR, SIG_NAME, "_", METHOD_NAME, "_Score_Distributions_Histogram_Density.png"),
         plot = p, width = 8, height = 6)


df2 <- data.frame(
  Platelet_AUC = pbmc$Raw_Score,
  Immune_AUC   = pbmc$Immune_Score,
  Prediction   = pbmc$Prediction
) %>%
  pivot_longer(cols = c(Platelet_AUC, Immune_AUC),
               names_to = "Score", values_to = "Value")

p <- ggplot(df2, aes(x = Prediction, y = Value, fill = Prediction)) +
  geom_violin(trim = FALSE) +
  facet_wrap(~ Score, scales = "free_y") +
  theme_classic() +
  labs(title = "Score distributions by Prediction",
       x = NULL, y = "Score value")

ggsave(filename = paste0(OUT_DIR, SIG_NAME, "_", METHOD_NAME, "_Score_Distributions_by_Prediction.png"),
         plot = p, width = 8, height = 6)
