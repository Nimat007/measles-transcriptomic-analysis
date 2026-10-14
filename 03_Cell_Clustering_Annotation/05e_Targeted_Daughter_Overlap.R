# ==============================================================================
# SCRIPT 05e: TARGETED DAUGHTER OVERLAP DIAGNOSTIC (CLUSTERS 5 & 8)
# ==============================================================================

library(here)
library(Seurat)
library(tidyverse)
library(patchwork)

# ------------------------------------------------------------------------------
# 1. Setup & Load Existing Data
# ------------------------------------------------------------------------------
message("--> Loading provisional object and existing marker tables...")

dir_processed <- here::here("data", "processed")
dir_results   <- here::here("results", "02_single_cell", "substructure_audit_v2")
dir_plots     <- here::here("results", "plots", "substructure_audit_v2")

seu <- readRDS(file.path(dir_processed, "GSE168775_Clustered_Provisional.rds"))
parent_markers <- read.csv(here::here("results", "02_single_cell", "Positive_Cluster_Markers.csv"))
daughter_markers <- read.csv(file.path(dir_results, "TaskD_Daughter_vs_Daughter_Markers.csv"))

# Define the highly specific target programs requested by the reviewer
target_programs <- list(
  ISG_Program = c("IFIT1", "IFIT2", "IFIT3", "ISG15", "MX1", "OAS1", "RSAD2"),
  Targeted_Epithelial = c("PIGR", "ATP12A", "TFCP2L1", "CAPS")
)

# ------------------------------------------------------------------------------
# 2. Add Targeted Module Scores
# ------------------------------------------------------------------------------
message("--> Computing targeted Module Scores for 0.8 daughters...")

DefaultAssay(seu) <- "SCT"

prog_names_clean <- c()
for (prog_name in names(target_programs)) {
  genes <- intersect(target_programs[[prog_name]], rownames(seu))
  if (length(genes) > 0) {
    seu <- AddModuleScore(seu, features = list(genes), name = prog_name, assay = "SCT", ctrl = 5)
    prog_names_clean <- c(prog_names_clean, paste0(prog_name, "1"))
  }
}

# ------------------------------------------------------------------------------
# 3. Analyze Parent 8 (ISG Program across Daughters 13 & 16)
# ------------------------------------------------------------------------------
message("--> Analyzing Parent 8 Daughters (ISG Concentration)...")

cells_p8 <- rownames(seu@meta.data)[seu$rpca_snn_res.0.6 == "8"]
seu_p8 <- subset(seu, cells = cells_p8)

# Module Score Summary
p8_summary <- seu_p8@meta.data %>%
  group_by(Daughter_0.8 = rpca_snn_res.0.8) %>%
  summarise(
    Parent = "8",
    Total_Cells = n(),
    Mean_ISG_Score = round(mean(ISG_Program1), 3),
    .groups = "drop"
  )

# Marker overlap comparison (Parent 8 vs Daughters 13 and 16)
p8_parent_genes <- parent_markers %>% filter(cluster == "8") %>% pull(gene)
p8_d13_genes <- daughter_markers %>% filter(Parent_0.6 == "8", Daughter_0.8 == "13") %>% pull(gene)
p8_d16_genes <- daughter_markers %>% filter(Parent_0.6 == "8", Daughter_0.8 == "16") %>% pull(gene)

p8_marker_comp <- data.frame(
  Parent_0.6 = "8",
  Daughter_0.8 = c("13", "16"),
  Parent_Total_Markers = length(p8_parent_genes),
  Daughter_Specific_Markers = c(length(p8_d13_genes), length(p8_d16_genes)),
  Shared_with_Parent = c(length(intersect(p8_parent_genes, p8_d13_genes)), length(intersect(p8_parent_genes, p8_d16_genes)))
)

# ------------------------------------------------------------------------------
# 4. Analyze Parent 5 (Targeted Epithelial Program across Daughters 9 & 12)
# ------------------------------------------------------------------------------
message("--> Analyzing Parent 5 Daughters (Targeted Epithelial Concentration)...")

cells_p5 <- rownames(seu@meta.data)[seu$rpca_snn_res.0.6 == "5"]
seu_p5 <- subset(seu, cells = cells_p5)

# Module Score Summary
p5_summary <- seu_p5@meta.data %>%
  group_by(Daughter_0.8 = rpca_snn_res.0.8) %>%
  summarise(
    Parent = "5",
    Total_Cells = n(),
    Mean_Targeted_Epi_Score = round(mean(Targeted_Epithelial1), 3),
    .groups = "drop"
  )

# Marker overlap comparison (Parent 5 vs Daughters)
# (Note: Earlier diagnostics showed 0 DEGs between D9 and D12, but we log the parent baseline here)
p5_parent_genes <- parent_markers %>% filter(cluster == "5") %>% pull(gene)
p5_d9_genes <- daughter_markers %>% filter(Parent_0.6 == "5", Daughter_0.8 == "9") %>% pull(gene)

p5_marker_comp <- data.frame(
  Parent_0.6 = "5",
  Daughter_0.8 = c("9", "12"),
  Parent_Total_Markers = length(p5_parent_genes),
  Daughter_Specific_Markers = c(length(p5_d9_genes), 0), # Hardcoded 0 based on previous run warning
  Shared_with_Parent = c(length(intersect(p5_parent_genes, p5_d9_genes)), 0)
)

# ------------------------------------------------------------------------------
# 5. Export Summaries & Visualizations
# ------------------------------------------------------------------------------
message("--> Exporting Final Diagnostic Data & Plots...")

write.csv(bind_rows(p8_summary, p5_summary), file.path(dir_results, "TaskI_Targeted_Daughter_ModuleScores.csv"), row.names = FALSE)
write.csv(bind_rows(p8_marker_comp, p5_marker_comp), file.path(dir_results, "TaskJ_Parent_vs_Daughter_Marker_Overlap.csv"), row.names = FALSE)

# Generate simple Violin Plots for visual distribution of the scores
p_vln_p8 <- VlnPlot(seu_p8, features = "ISG_Program1", group.by = "rpca_snn_res.0.8", pt.size = 0.1) +
  ggtitle("Parent 8: ISG Program by 0.8 Daughter") + ylab("Module Score") + NoLegend()

p_vln_p5 <- VlnPlot(seu_p5, features = "Targeted_Epithelial1", group.by = "rpca_snn_res.0.8", pt.size = 0.1) +
  ggtitle("Parent 5: Targeted Epi Program by 0.8 Daughter") + ylab("Module Score") + NoLegend()

p_grid <- p_vln_p8 | p_vln_p5
ggsave(file.path(dir_plots, "TaskK_Daughter_Targeted_Violins.png"), p_grid, width = 12, height = 6, dpi = 300)

message("SUCCESS: Targeted daughter overlap diagnostic complete.")