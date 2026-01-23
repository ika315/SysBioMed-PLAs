# ============================================================
# Thresholding variants – AddModuleScore
# ============================================================

# ----------------------------
# Libraries
# ----------------------------
library(Seurat)
library(mclust)
library(pROC)
library(dplyr)
library(ggplot2)
library(Matrix)

# ----------------------------
# Paths & parameters
# ----------------------------
PATH_DATA <- "~/SysBioMed-PLAs/data/seu_sx_integration.rds"
GENE_LIST <- "~/SysBioMed-PLAs/data/updated_gene_list.csv"

plot_dir <- "plots/thresholding/AddModuleScore"
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

# ----------------------------
# Load helper functions
# ----------------------------
source("~/SysBioMed-PLAs/src/read_and_extend_gene_list.R")

# ----------------------------
# Load data (memory-safe)
# ----------------------------
message("Loading Seurat object...")

pbmc_full <- readRDS(PATH_DATA)
DefaultAssay(pbmc_full) <- "RNA"

pbmc <- DietSeurat(
  pbmc_full,
  assays = "RNA",
  reductions = c("pca", "umap"),
  graphs = c("RNA_nn", "RNA_snn"),
  misc = FALSE,
  images = FALSE
)

rm(pbmc_full)
gc()

# recompute PCA if needed
if (!"pca" %in% Reductions(pbmc)) {
  message("PCA not found – recomputing PCA")
  pbmc <- ScaleData(pbmc, verbose = FALSE)
  pbmc <- RunPCA(pbmc, npcs = 30, verbose = FALSE)
}

# recompute UMAP if needed
if (!"umap" %in% Reductions(pbmc)) {
  message("UMAP not found – recomputing UMAP")
  pbmc <- RunUMAP(pbmc, dims = 1:30)
}

# ----------------------------
# Ground truth (GT)
# ----------------------------
gt_col_name <- "celltype.l2"
pbmc$GT_Response <- ifelse(pbmc[[gt_col_name]] == "Platelet", 1, 0)
pbmc$GT_Class    <- ifelse(pbmc$GT_Response == 1, "Platelet", "Non-platelet")

pbmc$celltype.l1 <- factor(pbmc$celltype.l1)
pbmc$celltype.l2 <- factor(pbmc$celltype.l2)
pbmc$celltype.l3 <- factor(pbmc$celltype.l3)

# ----------------------------
# Load platelet gene set
# ----------------------------
genes <- read_gene_list(GENE_LIST)
sig_list <- list(Platelet_Signature = genes)

# ============================================================
# AddModuleScore scoring
# ============================================================

message("Computing AddModuleScore...")

pbmc <- AddModuleScore(
  pbmc,
  features = sig_list,
  name = "AMS"
)

pbmc$AMS_Score <- pbmc$AMS1

# ----------------------------
# Z-scoring
# ----------------------------
scores_raw <- pbmc$AMS_Score
mean_score <- mean(scores_raw)
sd_score   <- sd(scores_raw)

pbmc$AMS_ZScore <- (scores_raw - mean_score) / sd_score

# ============================================================
# Thresholding methods
# ============================================================

# ----------------------------
# 1) GMM (raw score)
# ----------------------------
gmm <- Mclust(scores_raw, G = 2)
pbmc$flag_gmm <- gmm$classification == which.max(gmm$parameters$mean)

gmm_means <- sort(gmm$parameters$mean)
thr_gmm <- mean(gmm_means)

# ----------------------------
# 2) MAD (raw score)
# ----------------------------
thr_mad <- median(scores_raw) + 3 * mad(scores_raw)
pbmc$flag_mad <- scores_raw > thr_mad

# ----------------------------
# 3) Quantile (raw score)
# ----------------------------
thr_q99 <- quantile(scores_raw, 0.99)
pbmc$flag_q99 <- scores_raw > thr_q99

# ----------------------------
# 4) Z-score elbow (unsupervised)
# ----------------------------
z <- pbmc$AMS_ZScore
z_sorted <- sort(z, decreasing = TRUE)

d1 <- diff(z_sorted)
d2 <- diff(d1)
knee <- which.min(d2) + 1

thr_z_elbow <- z_sorted[knee]
pbmc$flag_z_elbow <- z > thr_z_elbow

# ----------------------------
# 5) Youden index (supervised benchmark)
# ----------------------------
roc_obj <- roc(
  response  = pbmc$GT_Response,
  predictor = pbmc$AMS_ZScore,
  direction = "<",
  quiet = TRUE
)

youden_coords <- coords(roc_obj, x = "best", best.method = "youden")
thr_youden <- as.numeric(youden_coords$threshold)

pbmc$flag_youden <- pbmc$AMS_ZScore > thr_youden

# ============================================================
# Visualization
# ============================================================

flag_methods <- c(
  "flag_gmm",
  "flag_mad",
  "flag_q99",
  "flag_z_elbow",
  "flag_youden"
)

# ----------------------------
# UMAP: raw score
# ----------------------------
p_raw <- FeaturePlot(
  pbmc,
  features = "AMS_Score",
  reduction = "umap",
  pt.size = 0.4
) + ggtitle("AddModuleScore platelet signature (raw)")

ggsave(
  file.path(plot_dir, "UMAP_AMS_raw_score.png"),
  p_raw, width = 7, height = 6
)

# ----------------------------
# UMAPs: threshold flags
# ----------------------------
for (flag in flag_methods) {
  p_flag <- DimPlot(
    pbmc,
    group.by = flag,
    reduction = "umap",
    pt.size = 0.4
  ) + ggtitle(paste("Platelet-enriched cells:", flag))
  
  ggsave(
    file.path(plot_dir, paste0("UMAP_", flag, ".png")),
    p_flag, width = 7, height = 6
  )
}

# ----------------------------
# Violin: raw score + thresholds
# ----------------------------
thr_z_raw <- mean_score + thr_z_elbow * sd_score

p_vln_raw <- VlnPlot(pbmc, features = "AMS_Score", pt.size = 0) +
  geom_hline(yintercept = thr_mad,  linetype="dashed", color="blue") +
  geom_hline(yintercept = thr_q99,  linetype="dashed", color="purple") +
  geom_hline(yintercept = thr_gmm,  linetype="dashed", color="red") +
  geom_hline(yintercept = thr_z_raw, linetype="dashed", color="darkgreen") +
  labs(
    title = "AddModuleScore raw score distribution",
    subtitle = "Blue=MAD | Purple=Q99 | Red=GMM | Green=Z-elbow (projected)"
  )

ggsave(
  file.path(plot_dir, "VlnPlot_AMS_raw_score_thresholds.png"),
  p_vln_raw, width = 7, height = 6
)

# ----------------------------
# Violin: Z-score + thresholds
# ----------------------------
p_vln_z <- VlnPlot(pbmc, features = "AMS_ZScore", pt.size = 0) +
  geom_hline(yintercept = thr_z_elbow, linetype="dashed", color="darkgreen") +
  geom_hline(yintercept = thr_youden, linetype="dashed", color="black") +
  labs(
    title = "AddModuleScore Z-score distribution",
    subtitle = "Green=Z-elbow | Black=Youden (GT)"
  )

ggsave(
  file.path(plot_dir, "VlnPlot_AMS_ZScore_thresholds.png"),
  p_vln_z, width = 7, height = 6
)

# ----------------------------
# Error types (Youden)
# ----------------------------
pbmc$Error_Type_Youden <- case_when(
  pbmc$flag_youden & pbmc$GT_Response == 1 ~ "TP",
  pbmc$flag_youden & pbmc$GT_Response == 0 ~ "FP",
  !pbmc$flag_youden & pbmc$GT_Response == 1 ~ "FN",
  !pbmc$flag_youden & pbmc$GT_Response == 0 ~ "TN"
)

pbmc$Error_Type_Youden <- factor(
  pbmc$Error_Type_Youden,
  levels = c("TP", "FP", "FN", "TN")
)

# ----------------------------
# Error UMAP
# ----------------------------
p_err <- DimPlot(
  pbmc,
  group.by = "Error_Type_Youden",
  reduction = "umap"
) +
  scale_color_manual(values = c(
    "TP"="#228B22", "FP"="#FF4500",
    "FN"="#1E90FF", "TN"="#D3D3D3"
  )) +
  labs(
    title = "AddModuleScore error map (Youden)",
    subtitle = paste("Z threshold =", round(thr_youden, 2))
  )

ggsave(
  file.path(plot_dir, "UMAP_AMS_Youden_Errors.png"),
  p_err, width = 7, height = 6
)

# ----------------------------
# Density by error class
# ----------------------------
p_density <- ggplot(pbmc@meta.data, aes(x = AMS_ZScore)) +
  geom_density(fill="steelblue", alpha=0.6) +
  facet_wrap(~ Error_Type_Youden, scales="free_y") +
  geom_vline(xintercept = thr_youden, linetype="dashed") +
  theme_classic() +
  labs(title="AddModuleScore Z-score density by error class (Youden)")

ggsave(
  file.path(plot_dir, "AMS_ZScore_Density_Youden.png"),
  p_density, width = 10, height = 6
)

message("AddModuleScore thresholding & benchmarking complete.")