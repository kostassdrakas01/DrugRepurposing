#!/usr/bin/env Rscript

# Analytical Discovery Plotting (V3)
# Goal: Actionable analysis of Target-Pathway dependencies

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
  library(tidyr)
  library(scales)
})

# Load data
df <- read.csv(input_file)

# 1. Target-Pathway Intersection Analysis (Heatmap)
# This shows which proteins are actually driving each discovery
df_long <- df %>%
  filter(discovery_score >= 0.70) %>% # Focus on relevant hits
  separate_rows(targets, sep = ";") %>%
  filter(targets != "")

p1 <- ggplot(df_long, aes(x = targets, y = reorder(pathway_name, discovery_score))) +
  geom_tile(aes(fill = discovery_score), color = "white") +
  scale_fill_gradientn(colors = c("#e0f3f8", "#abd9e9", "#4575b4", "#313695"), limits = c(0.7, 1.0)) +
  theme_minimal() +
  labs(
    title = paste("Target-Pathway Dependency Map:", drug_name),
    subtitle = "Analysis of which targets drive specific biological discoveries",
    x = "Molecular Target (Protein)",
    y = "Biological Pathway",
    fill = "Confidence"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()
  )

heatmap_path <- file.path(output_dir, paste0(drug_name, "_target_interaction.png"))
ggsave(heatmap_path, p1, width = 12, height = 8, dpi = 300)

# 2. Quadrant Analysis: Statistical vs Biological Novelty
# High Z-Score = Statistically unlikely to be random
# High Surprise = Cross-category (repurposing potential)
p2 <- ggplot(df, aes(x = surprise_score, y = z_score)) +
  geom_vline(xintercept = 0.5, linetype = "dotted", color = "grey") +
  geom_hline(yintercept = 2.0, linetype = "dotted", color = "grey") +
  geom_point(aes(color = category, size = discovery_score), alpha = 0.6) +
  geom_text(data = subset(df, z_score > 4 | surprise_score > 0.9), 
            aes(label = pathway_name), vjust = -1, size = 2.5, check_overlap = TRUE) +
  theme_minimal() +
  annotate("text", x = 0.25, y = 10, label = "Clinical Indication (Known)", color = "grey", alpha = 0.5) +
  annotate("text", x = 0.75, y = 10, label = "High-Priority Discovery", color = "#4575b4", fontface = "bold") +
  annotate("text", x = 0.25, y = 0.5, label = "Low-Confidence / Noise", color = "grey", alpha = 0.5) +
  annotate("text", x = 0.75, y = 0.5, label = "Stochastic Surprise", color = "grey", alpha = 0.5) +
  labs(
    title = paste("Discovery Quadrant Analysis:", drug_name),
    subtitle = "Mapping Statistical Significance against Biological Surprise",
    x = "Biological Surprise (Novelty)",
    y = "Statistical Significance (Z-Score)",
    color = "System Category",
    size = "Composite Score"
  ) +
  theme(legend.position = "right")

quadrant_path <- file.path(output_dir, paste0(drug_name, "_quadrant_analysis.png"))
ggsave(quadrant_path, p2, width = 11, height = 7, dpi = 300)

cat(paste0("[+] Analytical plots complete for ", drug_name, "\n"))
cat(paste0("[+] Actionable Analysis Visuals:\n    - ", heatmap_path, " (Target Dependency Map)\n    - ", quadrant_path, " (Quadrant Strategy Map)\n"))
