# ==============================================================================
# SCRIPT 06a: DUAL-TIER BIOLOGICAL ANNOTATION
# ==============================================================================

library(here)
library(Seurat)
library(tidyverse)
library(patchwork)

dir_processed <- here::here("data", "processed")
dir_results   <- here::here("results", "02_single_cell", "annotation")
dir_plots     <- here::here("results", "plots", "annotation")

if (!dir.exists(dir_results)) dir.create(dir_results, recursive = TRUE)
if (!dir.exists(dir_plots)) dir.create(dir_plots, recursive = TRUE)

message("--> Loading Provisional 0.6 Clustered Object...")
seu <- readRDS(file.path(dir_processed, "GSE168775_Clustered_Provisional.rds"))

# ------------------------------------------------------------------------------
# 1. Define Dual-Tier Annotation Map (Strict/Conservative)
# ------------------------------------------------------------------------------
anno_map <- data.frame(
  seurat_clusters = as.character(0:15),
  Broad_Lineage = c(
    "Basal",       # 0
    "Unresolved",  # 1
    "Secretory",   # 2
    "Unresolved",  # 3
    "Unresolved",  # 4
    "Secretory",   # 5
    "Secretory",   # 6
    "Unresolved",  # 7
    "Secretory",   # 8
    "Secretory",   # 9
    "Ciliated",    # 10
    "Secretory",   # 11
    "Unresolved",  # 12
    "Ciliated",    # 13
    "Ciliated",    # 14
    "Unresolved"   # 15
  ),
  High_Res_State = c(
    "Basal",                      # 0
    "Unresolved_1",               # 1
    "Secretory_Gradient_A",       # 2
    "Unresolved_3",               # 3
    "Remodeling_Associated",      # 4
    "Secretory",                  # 5
    "Secretory_Inflammatory",     # 6
    "ISG_High_Epithelial",        # 7
    "Secretory_Gradient_B",       # 8
    "Club_Secretory",             # 9
    "Ciliated",                   # 10
    "Secretory_Associated",       # 11
    "Infection_Associated_State", # 12
    "Ciliated_Related",           # 13
    "Deuterosomal_Precursor",     # 14
    "Cycling_Proliferating"       # 15
  ),
  stringsAsFactors = FALSE
)

# ------------------------------------------------------------------------------
# 2. Pre-Mapping Integrity Checks
# ------------------------------------------------------------------------------
message("--> Running Annotation Integrity Checks...")

if (any(duplicated(anno_map$seurat_clusters))) {
  stop("CRITICAL ERROR: anno_map contains duplicated seurat_clusters.")
}

existing_clusters <- unique(as.character(seu$seurat_clusters))
if (!all(existing_clusters %in% anno_map$seurat_clusters)) {
  stop("CRITICAL ERROR: One or more seurat_clusters in the object are missing from anno_map.")
}

write.csv(anno_map, file.path(dir_results, "Cluster_Annotation_Dictionary.csv"), row.names = FALSE)
unresolved_list <- anno_map %>% filter(Broad_Lineage == "Unresolved")
write.csv(unresolved_list, file.path(dir_results, "Unresolved_Clusters_Report.csv"), row.names = FALSE)

# ------------------------------------------------------------------------------
# 3. Apply Annotations (Barcode Safe)
# ------------------------------------------------------------------------------
message("--> Applying Annotations to Metadata...")

# Explicit index matching guarantees cell barcodes align perfectly
match_idx <- match(as.character(seu$seurat_clusters), anno_map$seurat_clusters)

seu$Broad_Lineage <- anno_map$Broad_Lineage[match_idx]
seu$High_Res_State <- anno_map$High_Res_State[match_idx]

# Post-Mapping Integrity Check
if (any(is.na(seu$Broad_Lineage)) || any(is.na(seu$High_Res_State))) {
  stop("CRITICAL ERROR: One or more cells received an NA annotation during mapping.")
}

# ------------------------------------------------------------------------------
# 4. Export Annotated Visualizations
# ------------------------------------------------------------------------------
message("--> Generating Annotated UMAPs...")

p_broad <- DimPlot(seu, group.by = "Broad_Lineage", label = TRUE, repel = TRUE) + 
  ggtitle("Broad Epithelial Lineages")
p_state <- DimPlot(seu, group.by = "High_Res_State", label = TRUE, repel = TRUE, label.size = 3) + 
  ggtitle("High Resolution Biological States")

ggsave(file.path(dir_plots, "UMAP_Broad_Lineage.png"), p_broad, width = 8, height = 6, dpi = 300)
ggsave(file.path(dir_plots, "UMAP_High_Res_State.png"), p_state, width = 10, height = 7, dpi = 300)

message("--> Saving Annotated Object...")
saveRDS(seu, file.path(dir_processed, "GSE168775_Annotated_Provisional.rds"))
message("SUCCESS: Script 06a Complete.")