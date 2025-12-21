# ============================================================
# Diagnostics for VISION Naive B Cell Benchmarking
# Purpose:
#   Determine whether poor performance is caused by
#   (1) data / signature / setup issues, or
#   (2) intrinsic method–biology mismatch
# ============================================================

library(Seurat)
library(VISION)
library(dplyr)
library(ggplot2)
library(Matrix)

# ------------------------------------------------------------
# 0. Setup
# ------------------------------------------------------------

pbmc <- readRDS("results/pbmc_seu_sx_final_Vision_NaiveB.rds")
DefaultAssay(pbmc) <- "RNA_v3"

dir.create("diagnostics", showWarnings = FALSE, recursive = TRUE)

message("=== Loaded PBMC object ===")
print(pbmc)

# ============================================================
# 1. BASIC STRUCTURAL SANITY CHECKS
# ============================================================
# What this tests:
#   - Whether the Seurat object is intact
#   - Whether metadata needed for benchmarking exists
#
# Expected if OK:
#   - All required columns present
#   - Reasonable number of cells and genes
#
# If this fails:
#   - Pipeline or saving error (not a VISION issue)
# ============================================================

required_cols <- c(
  "Vision_Raw", "Vision_ZScore",
  "celltype.l1", "celltype.l2",
  "GT_Response", "Error_Type"
)

missing_cols <- setdiff(required_cols, colnames(pbmc@meta.data))
message("Missing metadata columns:")
print(missing_cols)

# ============================================================
# 2. GROUND TRUTH BALANCE
# ============================================================
# What this tests:
#   - Class imbalance of Naive B cells
#
# Expected if OK:
#   - Naive B cells ideally >5–10% of dataset
#
# If this fails:
#   - Precision–Recall curves and thresholds become unstable
#   - Method may appear worse than it is
# ============================================================

gt_tab <- table(pbmc$GT_Response)
print(gt_tab)
print(prop.table(gt_tab))

# ============================================================
# 3. SIGNATURE COVERAGE & SPARSITY (ROBUST VERSION)
# ============================================================
# What this tests:
#   - How many signature genes are present in the dataset
#   - Whether those genes are broadly expressed
#
# Expected if OK:
#   - ≥ 30–50 genes
#   - Most genes expressed in a reasonable fraction of cells
#
# If this fails:
#   - Signature is too sparse or mismatched (ID issue)
# ============================================================

# IMPORTANT: explicitly define which gene list we are testing
signature_genes <- genes   # or extended_genes if testing extended

expr_genes <- rownames(pbmc@assays$RNA_v3)

sig_genes <- intersect(expr_genes, signature_genes)

message("Number of signature genes found:")
print(length(sig_genes))

if (length(sig_genes) < 5) {
  warning(
    paste(
      "Too few signature genes matched the expression matrix (",
      length(sig_genes),
      "). Diagnostics beyond this point may be unreliable."
    )
  )
} else {
  
  expr_fraction <- rowMeans(
    pbmc@assays$RNA_v3@counts[sig_genes, , drop = FALSE] > 0
  )
  
  png("diagnostics/03_signature_expression_fraction.png", 800, 600)
  hist(
    expr_fraction,
    breaks = 30,
    main = "Fraction of cells expressing signature genes",
    xlab = "Fraction of cells"
  )
  dev.off()
  
  print(summary(expr_fraction))
}

# ============================================================
# 4. RAW VISION SCORE DISTRIBUTION
# ============================================================
# What this tests:
#   - Whether VISION scores have usable dynamic range
#
# Expected if OK:
#   - Smooth, unimodal-ish distribution
#
# Red flags:
#   - Extremely narrow range
#   - Spikes / discretization
# ============================================================

p_raw <- ggplot(pbmc@meta.data, aes(Vision_Raw)) +
  geom_density(fill = "grey70") +
  theme_minimal() +
  labs(title = "VISION raw score distribution")

ggsave(
  "diagnostics/04_vision_raw_distribution.png",
  p_raw, width = 7, height = 5
)

# ============================================================
# 5. Z-SCORE CALIBRATION
# ============================================================
# What this tests:
#   - Whether scaling behaves as expected
#
# Expected if OK:
#   - Mean ≈ 0, SD ≈ 1
#
# If SD << 1:
#   - No biological signal
# ============================================================

summary(pbmc$Vision_ZScore)
sd(pbmc$Vision_ZScore, na.rm = TRUE)

# ============================================================
# 6. BIOLOGICAL ALIGNMENT WITH CELL TYPES
# ============================================================
# What this tests:
#   - Whether Naive B cells actually score higher
#     than other B or immune populations
#
# Expected if OK:
#   - Clear upward shift for "B naive"
#
# If this fails:
#   - Signature not discriminative in this dataset
# ============================================================

p_ct <- ggplot(
  pbmc@meta.data,
  aes(x = celltype.l2, y = Vision_ZScore)
) +
  geom_boxplot(outlier.shape = NA) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "VISION Z-score by celltype.l2")

ggsave(
  "diagnostics/06_zscore_by_celltype_l2.png",
  p_ct, width = 10, height = 6
)

# ============================================================
# 7. ERROR-TYPE SEPARATION (KEY FAILURE DIAGNOSTIC)
# ============================================================
# What this tests:
#   - Whether TP > FP and TN > FN in score space
#
# Expected if OK:
#   - TP clearly higher than FP
#   - FN clearly lower than TN
#
# If FP ≈ TP and FN ≈ TN:
#   - Method captures generic signal, not target biology
# ============================================================

p_err <- ggplot(
  pbmc@meta.data,
  aes(Error_Type, Vision_ZScore, fill = Error_Type)
) +
  geom_boxplot(outlier.shape = NA) +
  theme_minimal() +
  labs(title = "VISION Z-score by error type")

ggsave(
  "diagnostics/07_zscore_by_error_type.png",
  p_err, width = 8, height = 6
)

# ============================================================
# 8. TECHNICAL COVARIATE CORRELATIONS
# ============================================================
# What this tests:
#   - Whether VISION scores are driven by sequencing depth
#
# Expected if OK:
#   - Weak or no correlation
#
# If strong correlation:
#   - Pooling amplifies technical effects
# ============================================================

tech_vars <- c("nCount_RNA", "nFeature_RNA", "percent.mt")

for (v in tech_vars) {
  p <- ggplot(pbmc@meta.data, aes(.data[[v]], Vision_Raw)) +
    geom_point(alpha = 0.1) +
    geom_smooth(method = "lm", se = FALSE) +
    theme_minimal() +
    labs(title = paste("VISION raw vs", v))
  
  ggsave(
    paste0("diagnostics/08_vision_vs_", v, ".png"),
    p, width = 6, height = 5
  )
}

# ============================================================
# 9. BASELINE COMPARISON: MEAN EXPRESSION SCORE
# ============================================================
# What this tests:
#   - Whether VISION adds information beyond simple averaging
#
# Expected if OK:
#   - Moderate correlation (VISION ≠ trivial mean)
#
# If correlation ≈ 1:
#   - VISION adds no value here
# ============================================================

pbmc$MeanExprScore <- Matrix::colMeans(
  pbmc@assays$RNA_v3@data[sig_genes, ]
)

cor_vision_mean <- cor(
  pbmc$Vision_Raw,
  pbmc$MeanExprScore,
  use = "complete.obs"
)

message("Correlation Vision vs mean expression:")
print(cor_vision_mean)

# ============================================================
# 10. OPTIONAL: B-CELL-ONLY DIAGNOSTIC
# ============================================================
# What this tests:
#   - Whether signal exists *within* B cells only
#
# Expected if OK:
#   - Improved separation inside B-cell compartment
# ============================================================

pbmc_B <- subset(pbmc, celltype.l1 == "B")

p_B <- ggplot(
  pbmc_B@meta.data,
  aes(celltype.l2, Vision_ZScore)
) +
  geom_boxplot(outlier.shape = NA) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "VISION Z-score within B cells only")

ggsave(
  "diagnostics/10_zscore_Bcells_only.png",
  p_B, width = 8, height = 5
)

message("=== VISION diagnostics completed ===")