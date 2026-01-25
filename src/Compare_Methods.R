#-------------------------------------------------------------------
# Compare_Methods.R
# Compare performance metrics of different PLA scoring variants
# -------------------------------------------------------------------

library(ggplot2)
library(dplyr)
library(tidyr)
library(fmsb)
library(ggrepel)

# --- CONFIGURATION ---
METHODS_TO_COMPARE <- c("AUCell", "UCell", "AddModuleScore")
CURRENT_SIG <- "HP_ABNORMAL"

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

png(filename = paste0(PLOT_DIR, "PLA_Comparison.png"), width = 1000, height = 750, res = 120)

ggplot(plot_df, aes(x = Method, y = Value, color = Method)) +
  geom_segment(aes(x = Method, xend = Method, y = 0, yend = Value), color = "grey") +
  geom_point(size = 5) +
  facet_wrap(~Metric, scales = "fixed") +
  ylim(0, 1) +
  theme_bw() + 
  labs(title = "Benchmark Performance Profile",
       subtitle = "Comparison of scoring algorithms (Ground Truth: PLA_Gating)",
       y = "Score Value", x = NULL) +
  scale_color_brewer(palette = "Dark2") +
  theme(strip.background = element_rect(fill = "grey90"),
        strip.text = element_text(face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")

dev.off()

# Radar Chart
radar_df <- comparison_df %>% 
  select(AUC, Precision, Recall, F1_Score)
rownames(radar_df) <- comparison_df$Method

radar_df <- rbind(rep(1, 4), rep(0, 4), radar_df)

png(filename = paste0(PLOT_DIR, "PLA_Performance_Radar.png"), width = 800, height = 800, res = 120)
colors_border <- c("#1B9E77", "#D95F02", "#7570B3")
colors_in     <- c(rgb(0.1, 0.6, 0.5, 0.2), rgb(0.8, 0.4, 0.1, 0.2), rgb(0.4, 0.4, 0.7, 0.2))

radarchart(radar_df, pcol=colors_border, pfcol=colors_in, plwd=3, plty=1,
           cglcol="grey", cglty=1, axislabcol="grey", caxislabels=seq(0,1,0.25), cglwd=0.8,
           vlcex=1.2, title="Multi-Metric Method Profile")
legend(x=0.7, y=1.2, legend = rownames(radar_df)[-c(1,2)], bty = "n", pch=20, col=colors_border, cex=1.2, pt.cex=2)
dev.off()

# PRECISION-RECALL SCATTERPLOT
png(filename = paste0(PLOT_DIR, "PLA_Precision_vs_Recall_Scatter.png"), width = 900, height = 800, res = 120)

p_scatter <- ggplot(comparison_df, aes(x = Recall, y = Precision, color = Method, label = Method)) +
  geom_point(size = 6, alpha = 0.8) +
  geom_text_repel(size = 5, fontface = "bold") +
  geom_abline(intercept = 0, slope = 1, linetype = "dotted", color = "grey") + 
  xlim(0.5, 1) + ylim(0.5, 1) + # Fokus auf den relevanten Bereich
  labs(title = "Method Performance Trade-off",
       subtitle = "Top-Right corner represents ideal performance",
       x = "Recall (Sensitivity)", y = "Precision (Reliability)") +
  theme_bw() +
  scale_color_brewer(palette = "Set1")

print(p_scatter)
dev.off()


