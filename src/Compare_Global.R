# -------------------------------------------------------------------
# Compare_Global.R - Vergleich über alle Signaturen 
# -------------------------------------------------------------------

library(ggplot2)
library(dplyr)
library(tidyr)
library(ggrepel)

RESULTS_DIR <- "results/Global_Comparison/"
PLOT_DIR <- "plots/Global_Comparison/"
if (!dir.exists(PLOT_DIR)) dir.create(PLOT_DIR, recursive = TRUE)

all_files <- list.files(RESULTS_DIR, pattern = "^metrics_.*\\.csv$", full.names = TRUE)

all_data <- lapply(all_files, function(f) {
  df <- read.csv(f, stringsAsFactors = FALSE, check.names = FALSE)
  
  if (!"Signature" %in% colnames(df)) {
      fname <- basename(f)
      clean_name <- gsub("^metrics_|.csv$", "", fname)
      sig_part <- gsub("_[^_]+$", "", clean_name)
      df$Signature <- sig_part
  }
  
  colnames(df) <- trimws(colnames(df))
  
  return(df)
})

global_df <- bind_rows(all_data)

if (!"Signature" %in% colnames(global_df)) {
    stop("Kritischer Fehler: Spalte 'Signature' konnte nicht erzeugt werden!")
}

png(filename = paste0(PLOT_DIR, "Global_Method_Stability_F1.png"), width = 1000, height = 700)

ggplot(global_df, aes(x = Method, y = F1_Score, fill = Method)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA) +
  geom_jitter(aes(color = Signature), width = 0.1, size = 3) +
  theme_minimal() +
  ylim(0, 1) +
  labs(title = "Method Performance Stability",
       subtitle = "F1-Score across all tested Gene Signatures",
       y = "F1-Score (F-Measure)", x = "Scoring Method") +
  scale_fill_brewer(palette = "Set1")

dev.off()

png(filename = paste0(PLOT_DIR, "Global_Performance_Heatmap.png"), width = 1200, height = 800)

ggplot(global_df, aes(x = Method, y = Signature, fill = F1_Score)) +
  geom_tile(color = "white") +
  scale_fill_gradientn(colors = c("blue", "white", "red"), limits = c(0, 1)) +
  geom_text(aes(label = round(F1_Score, 2)), size = 5) +
  theme_minimal() +
  labs(title = "Global Performance Heatmap (F1-Score)",
       subtitle = "Rows: Signatures | Columns: Methods",
       fill = "F1-Score")

dev.off()

ranking <- global_df %>%
  arrange(desc(F1_Score)) %>%
  select(Signature, Method, F1_Score, AUC, Precision, Recall)

print("--- GLOBAL RANKING (Top 10) ---")
print(head(ranking, 10))

write.csv(ranking, paste0(RESULTS_DIR, "GLOBAL_FINAL_RANKING.csv"), row.names = FALSE)
