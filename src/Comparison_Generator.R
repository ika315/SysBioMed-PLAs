# -------------------------------------------------------------------
# Comparison_Generator.R
# Zentrales Skript für alle Vergleiche
# -------------------------------------------------------------------

library(ggplot2)
library(dplyr)
library(tidyr)
library(ggrepel)
library(pheatmap)
library(RColorBrewer)
library(purrr)

# Pfade
METRICS_DIR <- "results/metrics/"
CT_DIR      <- "results/celltype_data/"
BASE_OUT    <- "plots/Global_Comparison/"

# Ordner erstellen
dir.create(file.path(BASE_OUT, "Performance"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(BASE_OUT, "Heatmaps"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(BASE_OUT, "Runtime"), recursive = TRUE, showWarnings = FALSE)

# --- 1. DATEN EINLESEN ---

# A. Globale Metriken (AUC, F1, Runtime)
metric_files <- list.files(METRICS_DIR, pattern = "*.csv", full.names = TRUE)
all_metrics <- map_df(metric_files, read.csv) %>%
  mutate(Condition = paste(Method, ifelse(Ext, "Ext", "NoExt"), sep = "_"))

# B. Zelltyp-Daten (FP-Counts, Z-Scores)
ct_files <- list.files(CT_DIR, pattern = "*.csv", full.names = TRUE)
all_ct <- map_df(ct_files, read.csv) %>%
  mutate(Experiment = paste(Method, Signature, ifelse(Ext, "Ext", "NoExt"), Mode, sep = " | "))

plot_performance <- function(df, group_val, type = "Signature") {
  plot_df <- df %>%
    pivot_longer(cols = c(AUC, F1, Prec, Rec), names_to = "Metric", values_to = "Value")
  
  p <- ggplot(plot_df, aes(x = Condition, y = Value, color = Metric)) +
    geom_segment(aes(x = Condition, xend = Condition, y = 0, yend = Value), color = "lightgrey") +
    geom_point(size = 3) +
    facet_wrap(~Metric, scales = "free_y") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = paste("Performance for", group_val), x = NULL, y = "Score")
  
  ggsave(file.path(BASE_OUT, "Performance", paste0("Perf_", group_val, ".png")), p, width = 10, height = 6)
}

# Plots für jede Signatur erstellen
unique(all_metrics$Signature) %>% walk(~plot_performance(filter(all_metrics, Signature == .x), .x))

#Runtime Plot

png(file.path(BASE_OUT, "Runtime", "Method_Runtime_Comparison.png"), width = 800, height = 600, res = 120)
ggplot(all_metrics, aes(x = Method, y = Runtime, fill = Method)) +
  geom_boxplot(alpha = 0.7) +
  geom_jitter(width = 0.2) +
  theme_minimal() +
  labs(title = "Computational Efficiency", y = "Runtime (Minutes)", x = "Method")
dev.off()

#Heatmaps 

create_ct_heatmap <- function(df, ct_column, value_column, title_suffix) {
  plot_mat <- df %>%
    group_by(!!sym(ct_column), Experiment) %>%
    summarise(val = sum(!!sym(value_column)), .groups = "drop") %>%
    pivot_wider(names_from = Experiment, values_from = val, values_fill = 0)
  
  mat <- as.matrix(plot_mat[,-1])
  rownames(mat) <- plot_mat[[ct_column]]
  
  plot_val <- if(value_column == "FP_Count") log1p(mat) else mat
  
  png(file.path(BASE_OUT, "Heatmaps", paste0("Heatmap_", value_column, "_", title_suffix, ".png")), 
      width = 1500, height = 1000, res = 120)
  pheatmap(plot_val, 
           display_numbers = mat, 
           number_format = "%.0f",
           color = colorRampPalette(c("white", "orange", "red"))(100),
           main = paste("Global Comparison:", value_column, "(", title_suffix, ")"),
           clustering_distance_cols = "euclidean")
  dev.off()
}

# Heatmaps für FP_Count (Clean & L3)
create_ct_heatmap(all_ct, "celltype_clean", "FP_Count", "Gated_Celltypes")
create_ct_heatmap(all_ct, "celltype_l3", "FP_Count", "L3_Celltypes")

# Heatmaps für Mean_Z (Zelltyp-spezifische Score-Intensität)
create_ct_heatmap(all_ct, "celltype_clean", "Mean_Z", "ZScore_Gated")

print("Alle Vergleiche wurden erfolgreich erstellt!")
