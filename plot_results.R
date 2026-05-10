#!/usr/bin/env Rscript

# Bio-Network Discovery Plotting Script (V2)
# Usage: Rscript plot_results.R results/Aspirin.csv

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  stop("Please provide the path to a results CSV file (e.g., results/Aspirin.csv)")
}

input_file <- args[1]
output_dir <- dirname(input_file)
drug_name <- gsub(".csv", "", basename(input_file))

# Load libraries
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(scales)
})

# Load data
df <- read.csv(input_file)

# 1. Discovery Volcano Plot (Discovery Score vs Z-Score)
# We use Discovery Score (0-1) and Z-Score (Statistical significance)
p1 <- ggplot(df, aes(x = z_score, y = discovery_score, color = category)) +
  geom_point(aes(size = surprise_score), alpha = 0.7) +
  geom_hline(yintercept = 0.8, linetype = "dashed", color = "red") +
  geom_vline(xintercept = 2.0, linetype = "dashed", color = "blue") +
  theme_minimal() +
  labs(
    title = paste("Discovery Landscape:", drug_name),
    subtitle = "Red line: Discovery Threshold (0.8) | Blue line: Significance Threshold (Z=2.0)",
    x = "Statistical Significance (Z-Score)",
    y = "Composite Discovery Score (0-1)",
    size = "Biological Novelty",
    color = "System Category"
  ) +
  theme(legend.position = "bottom") +
  # Add labels for top hits
  geom_text(data = subset(df, discovery_score > 0.9 | z_score > 5), 
            aes(label = pathway_name), vjust = -1, size = 3, check_overlap = TRUE)

volcano_path <- file.path(output_dir, paste0(drug_name, "_r_volcano.png"))
ggsave(volcano_path, p1, width = 10, height = 7, dpi = 300)

# 2. System Category Distribution
p2 <- ggplot(df %>% filter(is_significant == "True"), aes(x = reorder(category, discovery_score), y = discovery_score)) +
  geom_boxplot(outlier.shape = NA, fill = "lightblue", alpha = 0.5) +
  geom_jitter(aes(color = polarity), width = 0.2) +
  coord_flip() +
  theme_minimal() +
  labs(
    title = paste("Significant Discoveries by Category:", drug_name),
    x = "Biological System",
    y = "Discovery Score",
    color = "Predicted Effect"
  )

dist_path <- file.path(output_dir, paste0(drug_name, "_r_category_dist.png"))
ggsave(dist_path, p2, width = 10, height = 6, dpi = 300)

cat(paste0("[+] R-based analysis complete for ", drug_name, "\n"))
cat(paste0("[+] Plots saved to:\n    - ", volcano_path, "\n    - ", dist_path, "\n"))
