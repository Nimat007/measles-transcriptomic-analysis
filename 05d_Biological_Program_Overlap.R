# ==============================================================================
# SCRIPT 05d: FINAL BIOLOGICAL-PROGRAM OVERLAP DIAGNOSTIC
# ==============================================================================

library(here)
library(Seurat)
library(tidyverse)
library(yaml)
library(patchwork)

# ------------------------------------------------------------------------------
# 1. Setup & Load Existing Data
# ------------------------------------------------------------------------------
message("--> Loading configuration, existing markers, and provisional object...")
config <- yaml::read_yaml(here::here("config.yml"))

dir_processed <- here::here(config$directories$data_processed)
dir_results   <- here::here("results", "02_single_cell")
dir_plots     <- here::here("results", "plots")

seu <- readRDS(file.path(dir_processed, "GSE168775_Clustered_Provisional.rds"))
emp_markers <- read.csv(file.path(dir_results, "Positive_Cluster_Markers.csv"))
target_programs <- config$scRNA_params$annotation$marker_panels

# ------------------------------------------------------------------------------
# 2. Empirical Marker Overlap Matrix
# ------------------------------------------------------------------------------
message("--> Computing discrete gene overlap between empirical markers and canonical programs...")

overlap_results <- list()

for (cluster_id in unique(emp_markers$cluster)) {
  cluster_genes <- emp_markers %>% filter(cluster == cluster_id) %>% pull(gene)
  
  for (prog_name in names(target_programs)) {
    prog_genes <- target_programs[[prog_name]]
    shared_genes <- intersect(cluster_genes, prog_genes)
    
    overlap_results[[paste(cluster_id, prog_name, sep = "_")]] <- data.frame(
      Cluster = cluster_id,
      Program = prog_name,
      Overlap_Count = length(shared_genes),
      Shared_Genes = paste(shared_genes, collapse = ", ")
    )
  }
}

overlap_df <- bind_rows(overlap_results) %>% filter(Overlap_Count > 0)
write.csv(overlap_df, file.path(dir_results, "Biological_Program_Marker_Overlap.csv"), row.names = FALSE)

# ------------------------------------------------------------------------------
# 3. Continuous Module Scoring (AddModuleScore)
# ------------------------------------------------------------------------------
message("--> Computing continuous Module Scores for targeted programs...")

# Ensure assay is SCT for accurate scoring
DefaultAssay(seu) <- "SCT"
Idents(seu) <- "seurat_clusters"

# Add module scores safely, accounting for genes present in the dataset
prog_names_clean <- c()
for (prog_name in names(target_programs)) {
  genes <- intersect(target_programs[[prog_name]], rownames(seu))
  if (length(genes) > 0) {
    seu <- AddModuleScore(seu, features = list(genes), name = prog_name, assay = "SCT", ctrl = 5)
    # Seurat automatically appends '1' to the name given in AddModuleScore
    prog_names_clean <- c(prog_names_clean, paste0(prog_name, "1"))
  }
}

# ------------------------------------------------------------------------------
# 4. Generate Aggregated Averages & DotPlot
# ------------------------------------------------------------------------------
message("--> Exporting Module Score summaries and visual diagnostics...")

# Extract mean module scores per cluster
mod_scores <- seu@meta.data %>%
  group_by(seurat_clusters) %>%
  summarise(across(all_of(prog_names_clean), mean), .groups = "drop")

# Clean column names for export
colnames(mod_scores) <- gsub("1$", "", colnames(mod_scores))
write.csv(mod_scores, file.path(dir_results, "Biological_Program_Mean_ModuleScores.csv"), row.names = FALSE)

# Generate a unified DotPlot of the module scores
p_dot <- DotPlot(seu, features = prog_names_clean, cluster.idents = TRUE) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10)) +
  labs(
    title = "Biological Program Overlap (Module Scores)",
    x = "Canonical Program",
    y = "0.6 Cluster ID"
  ) +
  scale_x_discrete(labels = function(x) gsub("1$", "", x))

ggsave(file.path(dir_plots, "Biological_Program_Overlap_DotPlot.png"), p_dot, width = 12, height = 7, dpi = 300)

message("SUCCESS: Biological-Program Overlap Diagnostic complete. Review CSVs and DotPlot.")