#-------------------------------------------------------------------
# Compare_Methods.R
# Compare performance metrics of different PLA scoring variants
# -------------------------------------------------------------------

library(ggplot2)
library(dplyr)
library(tidyr)

# --- CONFIGURATION ---
METHODS_TO_COMPARE <- c("AUCell", "UCell", "AddModuleScore")

# Path to results
RESULTS_DIR <- "results/"
PLOT_DIR <- "plots/Comparison/"
if (!dir.exists(PLOT_DIR)) dir.create(PLOT_DIR, recursive = TRUE)

# Read data
all_metrics <- list()

for (m in METHODS_TO_COMPARE) {
    file_path <- paste0(RESULTS_DIR, "metrics_", m, ".csv")
    if (file.exists(file_path)) {
        all_metrics[[m]] <- read.csv(file_path)
    } else {
        warning(paste("File not found for method:", m))
    }
}

comparison_df <- bind_rows(all_metrics)

plot_df <- comparison_df %>%
    select(Method, AUC, F1_Score, Precision, Recall) %>%
    pivot_longer(cols = -Method, names_to = "Metric", values_to = "Value")

# Plotting
png(filename = paste0(PLOT_DIR, "Method_Comparison_Barplot.png"), width = 1000, height = 700)

ggplot(plot_df, aes(x = Method, y = Value, fill = Method)) +
    geom_bar(stat = "identity", position = "dodge", alpha = 0.8) +
    facet_wrap(~Metric, scales = "fixed") +
    ylim(0, 1) +
    labs(title = "Performance Comparison of PLA Scoring Methods",
         subtitle = "Metrics based on Ground Truth (PLA_Gating)",
         y = "Value (0 to 1)", x = "") +
    theme_minimal() +
    theme(strip.text = element_text(face = "bold", size = 12), axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none") +
    scale_fill_brewer(palette = "Set1")

dev.off()

print(paste("Comparison plot generated at:", PLOT_DIR))
