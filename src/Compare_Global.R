# -------------------------------------------------------------------
# Compare_Global.R - Vergleich über alle Signaturen 
# -------------------------------------------------------------------

library(ggplot2)
library(dplyr)
library(tidyr)
library(ggrepel)

RESULTS_DIR <- "results/Global_Comparison/"
BASE_PLOT_DIR <- "plots/Global_Comparison/"
#if (!dir.exists(PLOT_DIR)) dir.create(PLOT_DIR, recursive = TRUE)

#all_files <- list.files(RESULTS_DIR, pattern = "^metrics_.*\\.csv$", full.names = TRUE)

#all_data <- lapply(all_files, function(f) {
#  df <- read.csv(f, stringsAsFactors = FALSE, check.names = FALSE)
#    if (!"Signature" %in% colnames(df)) {
#      fname <- basename(f)
#      clean_name <- gsub("^metrics_|.csv$", "", fname)
#      sig_part <- gsub("_[^_]+$", "", clean_name)
#      df$Signature <- sig_part
#  }
#  
#  colnames(df) <- trimws(colnames(df))
#  
#  return(df)
#})

all_files <- list.files(RESULTS_DIR, pattern = "^metrics_.*\\.csv$", full.names = TRUE)
all_data <- lapply(all_files, function(f) {
  df <- read.csv(f, stringsAsFactors = FALSE, check.names = FALSE)
  fname <- basename(f)
  clean_name <- gsub("^metrics_|.csv$", "", fname)
  parts <- strsplit(clean_name, "_")[[1]]
  df$Method <- parts[length(parts)]
  df$Signature <- paste(parts[1:(length(parts)-1)], collapse = "_")
  return(df)
})

global_df <- bind_rows(all_data)

create_plots <- function(data, group_name, sub_dir, x_axis_var) {
  current_dir <- file.path(BASE_PLOT_DIR, sub_dir, group_name)
  if (!dir.exists(current_dir)) dir.create(current_dir, recursive = TRUE)
  
  plot_df <- data %>%
    select(!!sym(x_axis_var), AUC, F1_Score, Precision, Recall) %>%
    pivot_longer(cols = -!!sym(x_axis_var), names_to = "Metric", values_to = "Value")

  png(filename = file.path(current_dir, "Performance_Lollipop.png"), width = 1000, height = 700, res = 120)
  p1 <- ggplot(plot_df, aes(x = !!sym(x_axis_var), y = Value, color = !!sym(x_axis_var))) +
    geom_segment(aes(x = !!sym(x_axis_var), xend = !!sym(x_axis_var), y = 0, yend = Value), color = "grey") +
    geom_point(size = 4) +
    facet_wrap(~Metric) +
    ylim(0, 1) + theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none") +
    labs(title = paste("Performance Profile:", group_name), y = "Score Value", x = NULL)
  print(p1)
  dev.off()

  png(filename = file.path(current_dir, "PR_Scatter.png"), width = 800, height = 700, res = 120)
  p2 <- ggplot(data, aes(x = Recall, y = Precision, color = !!sym(x_axis_var), label = !!sym(x_axis_var))) +
    geom_point(size = 5) + geom_text_repel() +
    geom_abline(intercept = 0, slope = 1, linetype = "dotted") +
    xlim(0, 1) + ylim(0, 1) + theme_minimal() +
    labs(title = paste("PR Trade-off:", group_name))
  print(p2)
  dev.off()
}

unique_sigs <- unique(global_df$Signature)
for (sig in unique_sigs) {
  sig_data <- filter(global_df, Signature == sig)
  create_plots(sig_data, sig, "By_Signature", "Method")
}

unique_methods <- unique(global_df$Method)
for (meth in unique_methods) {
  meth_data <- filter(global_df, Method == meth)
  create_plots(meth_data, meth, "By_Method", "Signature")
}

print("Alle Detail-Plots wurden in plots/Global_Comparisons/ erstellt!")



#png(filename = paste0(PLOT_DIR, "Global_Method_Stability_F1.png"), width = 1000, height = 700)

#ggplot(global_df, aes(x = Method, y = F1_Score, fill = Method)) +
#  geom_boxplot(alpha = 0.6, outlier.shape = NA) +
#  geom_jitter(aes(color = Signature), width = 0.1, size = 3) +
#  theme_minimal() +
#  ylim(0, 1) +
#  labs(title = "Method Performance Stability",
#       subtitle = "F1-Score across all tested Gene Signatures",
#       y = "F1-Score (F-Measure)", x = "Scoring Method") +
#  scale_fill_brewer(palette = "Set1")

#dev.off()

#png(filename = paste0(PLOT_DIR, "Global_Performance_Heatmap.png"), width = 1200, height = 800)

#ggplot(global_df, aes(x = Method, y = Signature, fill = F1_Score)) +
#  geom_tile(color = "white") +
#  scale_fill_gradientn(colors = c("blue", "white", "red"), limits = c(0, 1)) +
#  geom_text(aes(label = round(F1_Score, 2)), size = 5) +
#  theme_minimal() +
#  labs(title = "Global Performance Heatmap (F1-Score)",
#       subtitle = "Rows: Signatures | Columns: Methods",
#       fill = "F1-Score")

#dev.off()

#ranking <- global_df %>%
#  arrange(desc(F1_Score)) %>%
#  select(Signature, Method, F1_Score, AUC, Precision, Recall)

#print("--- GLOBAL RANKING (Top 10) ---")
#print(head(ranking, 10))

#write.csv(ranking, paste0(RESULTS_DIR, "GLOBAL_FINAL_RANKING.csv"), row.names = FALSE)
