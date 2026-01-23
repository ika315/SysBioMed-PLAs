# -------------------------------------------------------------------
# Benchmarking_Master_v2.R - FINAL VERSION
# -------------------------------------------------------------------
start_time <- Sys.time() 

library(Seurat)
library(AUCell)
library(UCell)
library(pROC)
library(dplyr)
library(ggplot2)
library(pheatmap)

args <- commandArgs(trailingOnly = TRUE)
METHOD_NAME    <- if(length(args) >= 1) args[1] else "AUCell"
SIG_NAME       <- if(length(args) >= 2) args[2] else "MANNE_UP"
SIG_FILE_BASE  <- if(length(args) >= 3) args[3] else "MANNE_COVID19_COMBINED_COHORT_VS_HEALTHY_DONOR_PLATELETS_UP.v2025.1.Hs"
USE_EXTENSION  <- if(length(args) >= 4) as.logical(args[4]) else TRUE
THRESH_MODE    <- if(length(args) >= 5) args[5] else "youden" 

# --- KONFIGURATION ---
GT_COLUMN    <- "pla.status"
POSITIVE_VAL <- "PLA"
PATH_DATA    <- "~/SysBioMed-PLAs/data/seu_sx_final.rds"

# --- ORDNERSTRUKTUR AUTOMATISCH ERSTELLEN ---
OUT_DIR <- paste0("plots/Platelet_Main/", SIG_NAME, "/", METHOD_NAME, "_Ext", USE_EXTENSION, "_", THRESH_MODE, "/")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create("results/metrics", recursive = TRUE, showWarnings = FALSE)
dir.create("results/celltype_data", recursive = TRUE, showWarnings = FALSE)
dir.create("results/extended_lists", recursive = TRUE, showWarnings = FALSE)

# --- DATEN LADEN ---
print(paste("Lade Daten für:", SIG_NAME, "mit", METHOD_NAME))
pbmc <- readRDS(PATH_DATA)
new_metadata <- read.csv("~/SysBioMed-PLAs/data/external_dataset_pla-status_metatable.csv", row.names = 1)
rownames(new_metadata) <- new_metadata$barcodes_clean
common_cells <- intersect(Cells(pbmc), rownames(new_metadata))
pbmc <- subset(pbmc, cells = common_cells)
pbmc <- AddMetaData(pbmc, metadata = new_metadata[common_cells, ])

# --- GENLISTE LADEN ---
base_dir <- getwd()
source(file.path(base_dir, "src", "read_and_extend_gene_list.R"))
PATH_SIG <- file.path(base_dir, "data", paste0(SIG_FILE_BASE, ".csv"))
genes <- read_gene_list(PATH_SIG)

# --- SCORING LOGIK ---
print(paste("--- Calculating Scores using", METHOD_NAME, "---"))

if (METHOD_NAME == "AUCell" || METHOD_NAME == "WeightedAUCell") {
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

# 2. Schritt: Gensequenz-Erweiterung (optional)
if (USE_EXTENSION) {
    EXT_FILE <- paste0("results/extended_lists/ext_", SIG_NAME, "_", METHOD_NAME, ".csv")
    if (file.exists(EXT_FILE)) {
        message("Lade existierende Liste...")
        extended_genes <- read.csv(EXT_FILE)$geneName
    } else {
        message("Berechne neue Extension...")
        res_ext <- extend_gene_set(pbmc, base_genes = genes, score_name = "Raw_Score_Original")
        extended_genes <- res_ext$extended_genes
        write.csv(data.frame(geneName = extended_genes), EXT_FILE, row.names = FALSE)
    }
    final_genes <- extended_genes
} else {
    final_genes <- genes
}

# 3. Schritt: Finales Scoring
if (METHOD_NAME == "AUCell" || METHOD_NAME == "WeightedAUCell") {
    pbmc$Raw_Score <- as.numeric(getAUC(AUCell_calcAUC(list(Platelet_Score = final_genes), rankings))[1, ])
} else if (METHOD_NAME == "UCell") {
    pbmc <- AddModuleScore_UCell(pbmc, features = list(Platelet_Score = final_genes), name = NULL)
    pbmc$Raw_Score <- pbmc$Platelet_Score
} else {
    pbmc <- AddModuleScore(pbmc, features = list(final_genes), name = "AMS")
    pbmc$Raw_Score <- pbmc$AMS1
}
pbmc$Z_Score <- as.vector(scale(pbmc$Raw_Score))

# --- THRESHOLDING & EVALUIERUNG ---
pbmc$GT_Response <- ifelse(pbmc[[GT_COLUMN]] == POSITIVE_VAL, 1, 0)
roc_obj <- roc(response = pbmc$GT_Response, predictor = pbmc$Z_Score, direction = "<", quiet = TRUE)

if (THRESH_MODE == "youden") {
    THRESHOLD_Z <- as.numeric(coords(roc_obj, x = "best", best.method = "youden")$threshold)
} else {
    # PLATZHALTER für die neue Logik
    THRESHOLD_Z <- 1.5 
}

pbmc$Prediction <- ifelse(pbmc$Z_Score > THRESHOLD_Z, "Positive", "Negative")
pbmc$Error_Type <- case_when(
    pbmc$Prediction == "Positive" & pbmc[[GT_COLUMN]] == POSITIVE_VAL ~ "TP",
    pbmc$Prediction == "Positive" & pbmc[[GT_COLUMN]] != POSITIVE_VAL ~ "FP",
    pbmc$Prediction == "Negative" & pbmc[[GT_COLUMN]] == POSITIVE_VAL ~ "FN",
    pbmc$Prediction == "Negative" & pbmc[[GT_COLUMN]] != POSITIVE_VAL ~ "TN"
)

# Zeit und Metriken
runtime_min <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))
tp <- sum(pbmc$Error_Type == "TP", na.rm = TRUE); fp <- sum(pbmc$Error_Type == "FP", na.rm = TRUE)
fn <- sum(pbmc$Error_Type == "FN", na.rm = TRUE); tn <- sum(pbmc$Error_Type == "TN", na.rm = TRUE)
prec <- if((tp + fp) > 0) tp / (tp + fp) else 0
rec  <- if((tp + fn) > 0) tp / (tp + fn) else 0
f1   <- if((prec + rec) > 0) 2 * prec * rec / (prec + rec) else 0

# Speichern
write.csv(data.frame(
    Method = METHOD_NAME, 
    Signature = SIG_NAME, 
    Ext = USE_EXTENSION, 
    Mode = THRESH_MODE, 
    AUC = as.numeric(auc(roc_obj)), 
    F1 = f1, 
    Prec = prec, 
    Rec = rec, 
    TP = tp, 
    FP = fp, 
    FN = fn, 
    TN = tn, 
    Threshold_Used = THRESH_MODE,
    Threshold_Value = THRESHOLD_Z,
    Runtime_Min = runtime_min
),
paste0("results/metrics/metrics_", SIG_NAME, "_", METHOD_NAME, "_Ext", USE_EXTENSION, "_", THRESH_MODE, ".csv"), row.names = FALSE)

ct_data <- pbmc@meta.data %>% 
    group_by(celltype_clean, celltype_l3) %>% 
    summarise(
        FP_Count = sum(Error_Type == "FP"), 
        TP_Count = sum(Error_Type == "TP"),
        Mean_Z = mean(Z_Score),
        Cell_Count = n(),
        .groups = "drop" 
    ) %>%
    mutate(
        Method = METHOD_NAME, 
        Signature = SIG_NAME, 
        Ext = USE_EXTENSION, 
        Mode = THRESH_MODE
    )

write.csv(ct_data, 
          paste0("results/celltype_data/ct_", SIG_NAME, "_", METHOD_NAME, "_Ext", USE_EXTENSION, "_", THRESH_MODE, ".csv"), 
          row.names = FALSE)

# --- PLOTS ---

# A. UMAPS
png(paste0(OUT_DIR, "1_Reference_UMAP.png"), 1200, 800)
print(DimPlot(pbmc, group.by="celltype_clean", label=T, repel=T) + labs(title="Reference Gating"))
dev.off()

png(paste0(OUT_DIR, "2_Score_UMAP.png"), 900, 700)
print(FeaturePlot(pbmc, features="Raw_Score") + scale_colour_viridis_c(option="magma") + labs(title="Raw Score Intensity"))
dev.off()

png(paste0(OUT_DIR, "3_Error_Mapping_UMAP.png"), 1000, 800)
print(DimPlot(pbmc, group.by="Error_Type") + scale_color_manual(values=c("TP"="#228B22","FP"="#FF4500","FN"="#1E90FF","TN"="#D3D3D3")) + 
      labs(title="Error Mapping", subtitle=paste("Threshold Z =", round(THRESHOLD_Z, 2))))
dev.off()

# B. VERTEILUNGEN 
# Density: PLA vs Platelet-free 
png(paste0(OUT_DIR, "4a_Density_ZScore_PLA_Status.png"), 1200, 800)
print(ggplot(pbmc@meta.data, aes(x=Z_Score, fill=!!sym(GT_COLUMN))) + 
      geom_density(alpha=0.5) + 
      facet_wrap(~celltype_clean) + 
      theme_minimal() + 
      labs(title="Z-Score Distribution: PLA vs Platelet-free", subtitle=paste("Signature:", SIG_NAME)))
dev.off()

png(paste0(OUT_DIR, "4b_Density_RawScore_PLA_Status.png"), 1200, 800)
print(ggplot(pbmc@meta.data, aes(x=Raw_Score, fill=!!sym(GT_COLUMN))) + 
      geom_density(alpha=0.5) + 
      facet_wrap(~celltype_clean) + 
      theme_minimal() + 
      labs(title="Raw Score Distribution: PLA vs Platelet-free", subtitle=paste("Method:", METHOD_NAME)))
dev.off()

# Density: Nur Platelets
png(paste0(OUT_DIR, "5a_Density_ZScore_Platelets_Only.png"), 1000, 600)
print(ggplot(subset(pbmc, celltype_l3 == "Platelets")@meta.data, aes(x=Z_Score)) + 
      geom_density(fill="red", alpha=0.4) + 
      theme_minimal() + 
      labs(title="Z-Score Distribution: Only Platelets (L3 Reference)"))
dev.off()

png(paste0(OUT_DIR, "5b_Density_RawScore_Platelets_Only.png"), 1000, 600)
print(ggplot(subset(pbmc, celltype_l3 == "Platelets")@meta.data, aes(x=Raw_Score)) + 
      geom_density(fill="darkred", alpha=0.4) + 
      theme_minimal() + 
      labs(title="Raw Score Distribution: Only Platelets (L3 Reference)"))
dev.off()

# Violin: Raw Score pro Zelltyp 
png(paste0(OUT_DIR, "6_Violin_RawScore_Celltypes.png"), 1200, 600)
print(VlnPlot(pbmc, features="Raw_Score", group.by="celltype_clean", pt.size=0) + labs(title="Raw Score Distribution"))
dev.off()

png(paste0(OUT_DIR, "6b_Split_Violin_PLA_Status.png"), 1400, 700)
print(VlnPlot(pbmc, features="Raw_Score", group.by="celltype_clean", split.by=GT_COLUMN, pt.size=0) +
      scale_fill_manual(values=c("PLA"="#FF4B4B", "platelet-free"="#4B8BFF")) + labs(title="Score Comparison: PLA vs. Free"))
dev.off()

# 7. Density nach Error-Klasse 
png(paste0(OUT_DIR, "7_Density_Error_Types_ZScore.png"), 1000, 600)
print(ggplot(pbmc@meta.data, aes(x=Z_Score)) + geom_density(fill="steelblue", alpha=0.6) + 
      facet_wrap(~Error_Type, scales="free_y") + geom_vline(xintercept=THRESHOLD_Z, linetype="dashed", color="red") + 
      theme_classic() + labs(title="Z-Score Density by Error Type"))
dev.off()

# 8. FP Count Barplot
fp_data <- pbmc@meta.data %>% filter(Error_Type == "FP") %>% group_by(celltype_clean) %>% tally() %>% arrange(desc(n))
if(nrow(fp_data) > 0) {
    png(paste0(OUT_DIR, "8_FP_Count_per_Celltype.png"), 1000, 700)
    print(ggplot(fp_data, aes(x=reorder(celltype_clean, -n), y=n, fill=celltype_clean)) + 
          geom_bar(stat="identity") + geom_text(aes(label=n), vjust=-0.5) + theme_minimal() + 
          theme(axis.text.x=element_text(angle=45, hjust=1), legend.position="none") + labs(title="False Positives per Celltype"))
    dev.off()
}

# C. PERFORMANCE KURVEN
png(paste0(OUT_DIR, "10_ROC_Curve.png"), 700, 700)
plot(roc_obj, col="#E41A1C", lwd=3, main=paste("ROC AUC:", round(auc(roc_obj),3)))
dev.off()

eval_df <- data.frame(score=pbmc$Z_Score, gt=pbmc$GT_Response) %>% arrange(desc(score)) %>%
           mutate(tp_c=cumsum(gt), fp_c=cumsum(1-gt), prec=tp_c/(tp_c+fp_c), rec=tp_c/sum(gt))
png(paste0(OUT_DIR, "11_PR_Curve.png"), 800, 700)
print(ggplot(eval_df, aes(x=rec, y=prec)) + geom_line(color="#E41A1C", linewidth=1.2) + 
      theme_minimal() + labs(title="Precision-Recall Curve", x="Recall", y="Precision"))
dev.off()

png(paste0(OUT_DIR, "12_Error_Composition_Bar.png"), 1200, 700)
err_perc <- pbmc@meta.data %>% group_by(celltype_clean, Error_Type) %>% tally() %>% 
            group_by(celltype_clean) %>% mutate(p=n/sum(n)*100)
print(ggplot(err_perc, aes(x=celltype_clean, y=p, fill=Error_Type)) + geom_bar(stat="identity") + 
      scale_fill_manual(values=c("TP"="#228B22","FP"="#FF4500","FN"="#1E90FF","TN"="#D3D3D3")) + 
      theme_minimal() + theme(axis.text.x=element_text(angle=45, hjust=1)) +
      labs(title="Error Composition per Cell Type", y="Percentage (%)", fill="Error Class"))
dev.off()

message("Done! Alle Plots in: ", OUT_DIR)
