# ==============================================================================
# SCRIPT 05c: TARGETED SUBSTRUCTURE DIAGNOSTIC v2 (CLUSTERS 2, 5, 8)
# ==============================================================================

library(here)
library(Seurat)
library(tidyverse)
library(patchwork)

# ------------------------------------------------------------------------------
# 1. Setup & Load Provisional Object
# ------------------------------------------------------------------------------
dir_processed <- here::here("data", "processed")
dir_results   <- here::here("results", "02_single_cell", "substructure_audit_v2")
dir_plots     <- here::here("results", "plots", "substructure_audit_v2")

if (!dir.exists(dir_results)) dir.create(dir_results, recursive = TRUE)
if (!dir.exists(dir_plots)) dir.create(dir_plots, recursive = TRUE)

message("--> Loading Provisional Clustered Object...")
seu <- readRDS(file.path(dir_processed, "GSE168775_Clustered_Provisional.rds"))

target_parents <- c("2", "5", "8")
meta_targets <- seu@meta.data %>% filter(rpca_snn_res.0.6 %in% target_parents)

# ------------------------------------------------------------------------------
# Task A & Summary: Identify daughters and apply size rule
# ------------------------------------------------------------------------------
message("--> Generating Daughter Mappings and Size Summaries...")
daughter_map <- meta_targets %>%
  group_by(Parent_0.6 = rpca_snn_res.0.6, Daughter_0.8 = rpca_snn_res.0.8) %>%
  summarise(Cell_Count = n(), .groups = "drop") %>%
  # Explicitly label the 50-cell rule as a practical interpretation screen, not a statistical validity cutoff
  mutate(Interpretation_Status = ifelse(Cell_Count >= 50, "Major (>=50 cells; Practical Screen)", "Minor Fragment (<50 cells; Practical Screen)")) %>%
  arrange(Parent_0.6, desc(Cell_Count))

write.csv(daughter_map, file.path(dir_results, "TaskA_Daughter_Mapping.csv"), row.names = FALSE)

# ------------------------------------------------------------------------------
# Task B: Parent-Specific Raw Library & Composition Summaries
# ------------------------------------------------------------------------------
message("--> Calculating Parent-Specific Library Compositions...")

comp_details <- list()
comp_summaries <- list()

for(p in target_parents) {
  p_meta <- meta_targets %>% filter(rpca_snn_res.0.6 == p)
  
  # Detail Table
  detail <- p_meta %>%
    group_by(Daughter_0.8 = rpca_snn_res.0.8, Library_ID) %>%
    summarise(Raw_Cell_Count = n(), .groups = "drop") %>%
    group_by(Daughter_0.8) %>%
    mutate(
      Parent_0.6 = p,
      Percent_of_Daughter = round((Raw_Cell_Count / sum(Raw_Cell_Count)) * 100, 2)
    ) %>%
    select(Parent_0.6, Daughter_0.8, Library_ID, Raw_Cell_Count, Percent_of_Daughter)
  
  comp_details[[p]] <- detail
  
  # Summary Table
  summary <- detail %>%
    group_by(Parent_0.6, Daughter_0.8) %>%
    summarise(
      Total_Cells = sum(Raw_Cell_Count),
      Number_of_Libraries = n_distinct(Library_ID),
      Largest_Library_Count = max(Raw_Cell_Count),
      Largest_Library_Percentage = max(Percent_of_Daughter),
      Number_of_Libraries_GE_10_cells = sum(Raw_Cell_Count >= 10),
      .groups = "drop"
    )
  
  comp_summaries[[p]] <- summary
}

write.csv(bind_rows(comp_details), file.path(dir_results, "TaskB_Parent_Specific_Composition.csv"), row.names = FALSE)
write.csv(bind_rows(comp_summaries), file.path(dir_results, "TaskB_Daughter_Size_Summary.csv"), row.names = FALSE)

# ------------------------------------------------------------------------------
# Task C: Day x Condition Composition (Raw Counts)
# ------------------------------------------------------------------------------
message("--> Calculating Day x Condition Compositions...")
day_cond_counts <- meta_targets %>%
  group_by(Parent_0.6 = rpca_snn_res.0.6, Daughter_0.8 = rpca_snn_res.0.8, Day, Condition) %>%
  summarise(Cell_Count = n(), .groups = "drop") %>%
  arrange(Parent_0.6, Daughter_0.8, Day, Condition)

write.csv(day_cond_counts, file.path(dir_results, "TaskC_Daughter_DayCondition_RawCounts.csv"), row.names = FALSE)

# ------------------------------------------------------------------------------
# Diagnostic Block: SCT Models & Cluster 5
# ------------------------------------------------------------------------------
message("\n==================================================")
message("   DIAGNOSTIC: SCT MODELS & CLUSTER 5")
message("==================================================")
message("Seurat version: ", packageVersion("Seurat"))

# Corrected SCTModel slot access for Seurat v5
message("SCT models before global PrepSCTFindMarkers: ", length(seu[["SCT"]]@SCTModel.list))

seu_marker <- seu
DefaultAssay(seu_marker) <- "SCT"
seu_marker <- PrepSCTFindMarkers(seu_marker, assay = "SCT", verbose = FALSE)

message("SCT models after global PrepSCTFindMarkers: ", length(seu_marker[["SCT"]]@SCTModel.list))

cl5_cells <- rownames(seu_marker@meta.data)[seu_marker$rpca_snn_res.0.6 == "5"]
cl5_meta <- seu_marker@meta.data[cl5_cells, ]

message("\nCluster 5 -> Daughter 9 library counts:")
print(table(cl5_meta$Library_ID[cl5_meta$rpca_snn_res.0.8 == "9"]))

message("\nCluster 5 -> Daughter 12 library counts:")
print(table(cl5_meta$Library_ID[cl5_meta$rpca_snn_res.0.8 == "12"]))

message(sprintf("\nMultiple libraries present in Daughter 9? : %s", length(unique(cl5_meta$Library_ID[cl5_meta$rpca_snn_res.0.8 == "9"])) > 1))
message(sprintf("Multiple libraries present in Daughter 12? : %s", length(unique(cl5_meta$Library_ID[cl5_meta$rpca_snn_res.0.8 == "12"])) > 1))
message("==================================================\n")

# ------------------------------------------------------------------------------
# Task D: Daughter-vs-Daughter Marker Analysis
# ------------------------------------------------------------------------------
message("--> Executing Daughter-vs-Daughter Marker Analysis...")
all_daughter_markers <- list()

for (p in target_parents) {
  message(sprintf("    Finding distinct markers for daughters of 0.6 Cluster %s...", p))
  
  cells_parent <- rownames(meta_targets)[meta_targets$rpca_snn_res.0.6 == p]
  seu_sub <- subset(seu_marker, cells = cells_parent)
  Idents(seu_sub) <- "rpca_snn_res.0.8"
  
  # Protect against minor fragments (< 3 cells causing strict mathematical failures)
  valid_idents <- names(which(table(Idents(seu_sub)) >= 3))
  seu_sub <- subset(seu_sub, idents = valid_idents)
  
  if(length(unique(Idents(seu_sub))) > 1) {
    # recorrect_umi = FALSE explicitly prevents the unequal library warning 
    # since data is already globally recorrected in seu_marker.
    m <- FindAllMarkers(seu_sub, assay = "SCT", recorrect_umi = FALSE, 
                        only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, verbose = FALSE)
    
    if(nrow(m) > 0) {
      m$Parent_0.6 <- p
      all_daughter_markers[[p]] <- m
    } else {
      message(sprintf("      NOTE: Valid DE executed, but 0 significant genes found between daughters of Cluster %s.", p))
    }
  }
}

if(length(all_daughter_markers) > 0) {
  df_daughter_markers <- bind_rows(all_daughter_markers) %>%
    select(Parent_0.6, cluster, gene, avg_log2FC, pct.1, pct.2, p_val_adj) %>%
    rename(Daughter_0.8 = cluster)
  
  write.csv(df_daughter_markers, file.path(dir_results, "TaskD_Daughter_vs_Daughter_Markers.csv"), row.names = FALSE)
} else {
  # Defensive empty DF creation if ALL parents fail to yield markers
  df_empty <- data.frame(Parent_0.6=character(), Daughter_0.8=character(), gene=character(), avg_log2FC=numeric(), pct.1=numeric(), pct.2=numeric(), p_val_adj=numeric())
  write.csv(df_empty, file.path(dir_results, "TaskD_Daughter_vs_Daughter_Markers.csv"), row.names = FALSE)
}

# ------------------------------------------------------------------------------
# Task E: Marker-Program Interpretation Panel
# ------------------------------------------------------------------------------
message("--> Generating Canonical Program Audits for Daughters...")
target_programs <- c(
  "KRT5", "TP63", "KRT14", "KRT15", "KRT17",          # Basal
  "FOXJ1", "TPPP3", "PIFO", "CAPS", "TUBB4B", "DNAH5",# Ciliated
  "DEUP1", "FOXN4", "CDC20B", "CCNO", "MCIDAS",       # Deuterosomal
  "SCGB1A1", "SCGB3A1", "BPIFA1", "BPIFB1",           # Secretory
  "MUC5B", "MUC5AC", "SPDEF", "AGR2", "CLCA1",        # Mucous
  "MKI67", "TOP2A", "PBK", "ASPM", "CENPF",           # Cycling
  "IFIT1", "IFIT2", "IFIT3", "ISG15", "MX1", "OAS1", "RSAD2" # ISG
)
avail_programs <- intersect(target_programs, rownames(seu))

seu_target_only <- subset(seu, cells = rownames(meta_targets))
Idents(seu_target_only) <- "rpca_snn_res.0.8"

p_dot <- DotPlot(seu_target_only, features = avail_programs) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8)) +
  labs(title = "Canonical Marker Programs across 0.8 Daughter Clusters", y = "0.8 Daughter ID")

ggsave(file.path(dir_plots, "TaskE_Daughter_Program_DotPlot.png"), p_dot, width = 16, height = 6, dpi = 300)

# ------------------------------------------------------------------------------
# Task F: Focused UMAP Visualization
# ------------------------------------------------------------------------------
message("--> Generating Focused UMAPs...")
umap_list <- list()

for (p in target_parents) {
  cells_parent <- rownames(meta_targets)[meta_targets$rpca_snn_res.0.6 == p]
  seu_sub <- subset(seu, cells = cells_parent)
  
  umap_list[[p]] <- DimPlot(seu_sub, reduction = "umap_rpca", group.by = "rpca_snn_res.0.8", pt.size = 1.2) +
    ggtitle(sprintf("0.6 Parent: %s", p)) +
    theme(legend.position = "bottom")
}

p_umap_grid <- wrap_plots(umap_list, ncol = 3)
ggsave(file.path(dir_plots, "TaskF_Focused_Daughter_UMAPs.png"), p_umap_grid, width = 15, height = 6, dpi = 300)

message("SUCCESS: v2 Substructure diagnostics generated. All results written to substructure_audit_v2.")