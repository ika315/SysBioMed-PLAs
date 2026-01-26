# -------------------------------------------------------------------
# Comparison_Generator_v3.R - FINAL Presentation Version
# -------------------------------------------------------------------
library(ggplot2)
library(dplyr)
library(tidyr)
library(purrr)
library(pheatmap)
library(RColorBrewer)
library(ggrepel)

# Pfade
METRICS_DIR <- "results/metrics/"
CT_DIR      <- "results/celltype_data/"
BASE_OUT    <- "plots/Global_Comparison/"

dir.create(BASE_OUT, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(BASE_OUT, "Heatmaps"), showWarnings = FALSE)

# --- 1. DATEN LADEN & ROBUSTES PARSING ---

parse_filename <- function(fname, type="metrics") {
  clean_name <- gsub(paste0(type, "_|.csv"), "", fname)
  parts <- strsplit(clean_name, "_")[[1]]
  # Format: Signature_Method_ExtStatus_Mode
  return(data.frame(
    Signature = parts[1],
    Method    = parts[2],
    Ext       = grepl("ExtTRUE", clean_name),
    Mode      = parts[length(parts)]
  ))
}

metric_files <- list.files(METRICS_DIR, pattern = "*.csv", full.names = TRUE)
all_metrics <- map_df(metric_files, function(f) {
  df <- read.csv(f)
  meta <- parse_filename(basename(f), "metrics")

  if(!"Method" %in% colnames(df)) df$Method <- meta$Method
  if(!"Signature" %in% colnames(df)) df$Signature <- meta$Signature
  if(!"Ext" %in% colnames(df)) df$Ext <- meta$Ext
  if(!"Mode" %in% colnames(df)) df$Mode <- meta$Mode
  return(df)
}) %>% mutate(
  Experiment = paste(Method, Signature, Mode, sep = " | "),
  Approach   = paste(Method, Mode, sep = " + ")
)

# Zelltyp-Daten laden
ct_files <- list.files(CT_DIR, pattern = "*.csv", full.names = TRUE)
all_ct <- map_df(ct_files, function(f) {
  df <- read.csv(f)
  meta <- parse_filename(basename(f), "ct")
  if(!"Method" %in% colnames(df)) df$Method <- meta$Method
  if(!"Signature" %in% colnames(df)) df$Signature <- meta$Signature
  if(!"Ext" %in% colnames(df)) df$Ext <- meta$Ext
  if(!"Mode" %in% colnames(df)) df$Mode <- meta$Mode
  df$Experiment <- paste(df$Method, df$Signature, df$Mode, sep = " | ")
  return(df)
})

# --- 2. PLOT: VERBESSERUNG DURCH EXTENSION ---
extension_impact <- all_metrics %>%
  select(Method, Signature, Mode, Ext, F1) %>%
  pivot_wider(names_from = Ext, values_from = F1, names_prefix = "Ext_") %>%
  filter(!is.na(Ext_TRUE) & !is.na(Ext_FALSE)) %>%
  mutate(Improvement = Ext_TRUE - Ext_FALSE)

png(file.path(BASE_OUT, "Roadmap_Extension_Benefit.png"), 1000, 700, res = 120)
print(ggplot(extension_impact, aes(x = Signature, y = Improvement, fill = Method)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  facet_wrap(~Mode) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Benefit of Gene Set Extension",
       subtitle = "Positive values indicate F1-Score increase with Extension",
       y = "Delta F1-Score (Ext - NoExt)", x = "Gene Signature"))
dev.off()

# --- 3. Computational Comparison ---

# 1. Normaler Boxplot
png(file.path(BASE_OUT, "Method_Runtime_Benchmarking.png"), 800, 600, res = 120)
p1 <- ggplot(filter(all_metrics, Ext == FALSE), aes(x = Method, y = Runtime_Min, fill = Method)) +
  geom_boxplot(alpha = 0.6) +
  geom_jitter(width = 0.2, alpha = 0.5) +
  theme_minimal() +
  labs(title = "Computational Efficiency per Method (Baseline)",
       y = "Runtime in Minutes", x = "Enrichment Method")
print(p1)
dev.off()

# 2. Log-Scale Boxplot
png(file.path(BASE_OUT, "Method_Runtime_Benchmarking_Log.png"), 800, 600, res = 120)
p2 <- ggplot(filter(all_metrics, Ext == FALSE), aes(x = Method, y = Runtime_Min, fill = Method)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA) + 
  geom_jitter(width = 0.2, alpha = 0.5) +
  scale_y_log10() + 
  theme_minimal() + 
  labs(title = "Computational Efficiency (Baseline Lists)", 
       subtitle = "Y-Axis on Log10 Scale",
       y = "Runtime in Minutes (Log)", x = "Enrichment Method")
print(p2)
dev.off()

# 3. Cleveland Dot Plot (FIXED)
png(file.path(BASE_OUT, "Method_Runtime_Benchmarking_Cleveland.png"), 800, 600, res = 120)
p3 <- ggplot(filter(all_metrics, Ext == FALSE), aes(x = Runtime_Min, y = Method, color = Method)) +
  geom_point(size = 3) +
  theme_light() +
  labs(title = "Method Efficiency Profile", x = "Runtime (Minutes)", y = NULL)
print(p3) 
dev.off()  

# 4. Bar Plot mit Error Bars (FIXED)
runtime_summary <- all_metrics %>%
  filter(Ext == FALSE) %>%
  group_by(Method) %>%
  summarise(mean_run = mean(Runtime_Min, na.rm = TRUE), 
            sd_run = sd(Runtime_Min, na.rm = TRUE), .groups = "drop")

png(file.path(BASE_OUT, "Method_Runtime_Benchmarking_Bar_Plot.png"), 800, 600, res = 120)
p4 <- ggplot(runtime_summary, aes(x = Method, y = mean_run, fill = Method)) +
  geom_bar(stat = "identity", alpha = 0.8, width = 0.6) +
  geom_errorbar(aes(ymin = mean_run - sd_run, ymax = mean_run + sd_run), width = 0.2) +
  theme_classic() +
  labs(title = "Average Computational Efficiency", y = "Mean Runtime (Min)")
print(p4) 
dev.off()  

# --- 4. HEATMAPS: ZELLTYP PERFORMANCE ---
df_heatmap <- all_ct %>% filter(Ext == TRUE)

create_presentation_heatmap <- function(data, ct_col, metric, filename) {
  if(!ct_col %in% colnames(data)) {
    message("Spalte ", ct_col, " nicht gefunden. Überspringe Heatmap.")
    return(NULL)
  }
  
  plot_mat <- data %>%
    group_by(!!sym(ct_col), Experiment) %>%
    summarise(val = mean(!!sym(metric), na.rm=T), .groups = "drop") %>%
    pivot_wider(names_from = Experiment, values_from = val, values_fill = 0)

  mat <- as.matrix(plot_mat[,-1])
  rownames(mat) <- plot_mat[[1]]

  mat[is.na(mat)] <- 0
  cluster_rows <- if(nrow(mat) > 1 && any(apply(mat, 1, sd) > 0)) TRUE else FALSE
  cluster_cols <- if(ncol(mat) > 1 && any(apply(mat, 2, sd) > 0)) TRUE else FALSE

  png(file.path(BASE_OUT, "Heatmaps", filename), 1800, 1100, res = 120)
  tryCatch({
    pheatmap(mat,
             color = colorRampPalette(c("white", "orange", "red"))(100),
             display_numbers = TRUE, 
             number_format = "%.2f",
             main = paste("Cell-type Specific:", metric, "for", ct_col),
             cluster_rows = cluster_rows,
             cluster_cols = cluster_cols,
             clustering_distance_cols = "euclidean",
             angle_col = 45)
  }, error = function(e) {
    pheatmap(mat, cluster_rows = FALSE, cluster_cols = FALSE, 
             display_numbers = TRUE, main = paste(metric, "- No Clustering"))
  })
  dev.off()
}

target_cell_columns <- c("celltype_clean", "celltype.l3")
target_modes <- c("youden", "gmm_dist_dual", "gmm_dist_platelet")
target_metrics <- c("FP", "TP", "Balanced_Accuracy", "Mean_Z")

for(ct_col in target_cell_columns) {
  for(m in target_modes) {
    for(met in target_metrics) {
      
      file_name <- paste0("Heatmap_", met, "_", m, "_", ct_col, ".png")
      
      create_presentation_heatmap(
        data = df_heatmap %>% filter(Mode == m),
        ct_col = ct_col,
        metric = met,
        filename = file_name
      )
    }
  }
}

# --- 5. Großer Vergleich: Precision vs Recall ---

#png(file.path(BASE_OUT, "Global_Approach_Comparison.png"), 1400, 900, res = 120)
#print(ggplot(all_metrics, aes(x = Prec, y = Rec, color = Mode, shape = Method)) +
#  geom_point(aes(size = F1), alpha = 0.7) +
#  facet_grid(Ext ~ Signature) +
#  theme_bw() +
#  scale_color_brewer(palette = "Set1") +
#  labs(title = "Universal Performance Comparison",
#       subtitle = "Precision vs Recall (Faceted by Signature and Extension)",
#       x = "Precision", y = "Recall (Sensitivity)",
#       size = "F1-Score", color = "Threshold Mode", shape = "Method"))
#dev.off()

png(file.path(BASE_OUT, "PR_Comparison_Incl_Youden.png"), 1400, 900, res = 120)
print(ggplot(filter(all_metrics, Mode %in% c("youden", "gmm_dist_dual", "gmm_dist_platelet")), 
             aes(x = Prec, y = Rec, color = Mode, shape = Method)) +
  geom_point(aes(size = F1), alpha = 0.7) +
  facet_grid(Ext ~ Signature) + theme_bw() + scale_color_brewer(palette = "Set1") +
  labs(title = "Performance: Supervised vs. Unsupervised"))
dev.off()

png(file.path(BASE_OUT, "PR_Comparison_GMM_Only.png"), 1400, 900, res = 120)
print(ggplot(filter(all_metrics, Mode %in% c("gmm_dist_dual", "gmm_dist_platelet")), 
             aes(x = Prec, y = Rec, color = Mode, shape = Method)) +
  geom_point(aes(size = F1), alpha = 0.7) +
  facet_grid(Ext ~ Signature) + theme_bw() + scale_color_manual(values=c("#377EB8", "#4DAF4A")) +
  labs(title = "Unsupervised Discovery: Platelet vs. Dual-Gate"))
dev.off()

# --- 6. Perfromance Scatterplot ---
dir.create(file.path(BASE_OUT, "Methods_Detail"), showWarnings = FALSE)

for(meth in unique(all_metrics$Method)) {
  meth_data <- all_metrics %>% filter(Method == meth, Ext == TRUE)
  if(nrow(meth_data) == 0) next
  
  png(file.path(BASE_OUT, "Methods_Detail", paste0("Profile_", meth, ".png")), 1200, 800, res = 120)
  p <- ggplot(meth_data, aes(x = Rec, y = Prec, color = Signature)) +
    geom_point(aes(size = F1), alpha = 0.8) + # Punkte etwas transparenter für Überlagerungen
    facet_wrap(~Mode) + 
    theme_bw() + 
    scale_size_continuous(range = c(2, 8)) + # Größere Punkte für bessere Sichtbarkeit
    labs(title = paste("Method Performance Profile:", meth), 
         subtitle = "Extended Lists Only | Faceted by Threshold Mode",
         x = "Recall", y = "Precision",
         size = "F1-Score", color = "Gene Signature") +
    theme(legend.position = "right")
  
  print(p)
  dev.off()
}

# --- 7. SIGNATURE RANKING (Lollipop) ---
png(file.path(BASE_OUT, "Signature_Performance_Ranking.png"), 1200, 800, res = 120)
summary_perf <- all_metrics %>%
  filter(Ext == TRUE) %>%
  group_by(Signature, Method) %>%
  summarise(mean_F1 = mean(F1, na.rm=T), .groups = "drop")

print(ggplot(summary_perf, aes(x = reorder(Signature, mean_F1), y = mean_F1, color = Method)) +
  geom_segment(aes(xend = Signature, yend = 0), size = 1, alpha = 0.3) +
  geom_point(size = 4) +
  coord_flip() +
  theme_minimal() +
  labs(title = "Average F1-Score per Signature (Extended Lists)",
       x = "Signature", y = "Mean F1-Score"))
dev.off()

print("Vergleichs-Plots wurden erfolgreich erstellt!")
