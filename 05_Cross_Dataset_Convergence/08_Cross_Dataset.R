# ==============================================================================
# SCRIPT 08: POST-PSEUDOBULK SUMMARY & CROSS-DATASET CONVERGENCE
# ==============================================================================
# 1. Install BiocManager if you don't already have it
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
# Force binary package installation to bypass Rtools compilation
options(pkgType = "binary")
BiocManager::install("h5mread", type = "binary")
# 2. Install ComplexHeatmap via Bioconductor
BiocManager::install("ComplexHeatmap")
library(here)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)
library(patchwork)
library(ComplexHeatmap)
library(RColorBrewer)

# ------------------------------------------------------------------------------
# DIRECTORIES & EXACT FILE PATHS
# ------------------------------------------------------------------------------
dir_pb_de   <- here::here("results", "04_pseudobulk_de")
dir_harm    <- here::here("results", "05_identifier_harmonization")
dir_results <- here::here("results", "05_cross_dataset_convergence")

if (!dir.exists(dir_results)) dir.create(dir_results, recursive = TRUE)

file_deg_prim <- file.path(dir_pb_de, "DEGs_Primary_Broad_Lineages.csv")
file_deg_sec  <- file.path(dir_pb_de, "DEGs_Secondary_Clusters.csv")
file_aud_prim <- file.path(dir_pb_de, "Audit_Primary_Model_Usage.csv")
file_aud_sec  <- file.path(dir_pb_de, "Audit_Secondary_Model_Usage.csv")
file_bulk_met <- file.path(dir_harm, "Bulk_Mapped_Gene_Metrics.csv")

# ------------------------------------------------------------------------------
# INPUT VALIDATION
# ------------------------------------------------------------------------------
message("--> Validating exact input files...")
req_files <- c(file_deg_prim, file_deg_sec, file_aud_prim, file_aud_sec, file_bulk_met)
if (any(!file.exists(req_files))) {
  stop("CRITICAL ERROR: Missing required inputs. Ensure Script 07 and Harmonization outputs are present.")
}

# Explicitly check Bulk_Mapped_Gene_Metrics.csv columns before loading completely
temp_bulk <- read.csv(file_bulk_met, nrows = 1)
required_bulk_cols <- c(
  "Mapped_Gene_Symbol",
  "Mean_Bulk_LFC",
  "Median_Bulk_LFC",
  "Direction_Consistency",
  "Number_Positive_Timepoints",
  "Number_Negative_Timepoints"
)

missing_bulk_cols <- setdiff(required_bulk_cols, colnames(temp_bulk))
if (length(missing_bulk_cols) > 0) {
  stop(
    paste(
      "CRITICAL ERROR: Missing required columns in Bulk_Mapped_Gene_Metrics.csv:",
      paste(missing_bulk_cols, collapse = ", ")
    )
  )
}

deg_prim <- read.csv(file_deg_prim)
deg_sec  <- read.csv(file_deg_sec)
aud_prim <- read.csv(file_aud_prim)
aud_sec  <- read.csv(file_aud_sec)

# Read authoritative mapped bulk metrics
# Read authoritative mapped bulk metrics
bulk_metrics <- read.csv(file_bulk_met) %>%
  dplyr::rename(Gene = "Mapped_Gene_Symbol") %>%
  dplyr::mutate(
    Bulk_Direction = case_when(
      Mean_Bulk_LFC > 0 ~ "Upregulated",
      Mean_Bulk_LFC < 0 ~ "Downregulated",
      TRUE ~ "No_Net_Change"
    )
  )

if (anyDuplicated(bulk_metrics$Gene) > 0) {
  stop("CRITICAL ERROR: Duplicate Gene entries detected in Bulk_Mapped_Gene_Metrics.csv.")
}

# ------------------------------------------------------------------------------
# COMBINE & MAP CONTRASTS
# ------------------------------------------------------------------------------
message("--> Mapping Contrasts and Joining Audit Data...")

deg_all <- bind_rows(
  deg_prim %>% mutate(Resolution = "Primary_Lineage", Target_Group = as.character(Target_Group)),
  deg_sec %>% mutate(Resolution = "Secondary_Cluster", Target_Group = as.character(Target_Group))
) %>%
  mutate(Mapped_Audit_Contrast = case_when(
    Contrast == "Bystander_vs_Mock" ~ "Bys_vs_Mock",
    Contrast == "Infected_vs_Bystander" ~ "Inf_vs_Bys",
    TRUE ~ "UNMAPPED_ERROR"
  ))

# Exclude viral features explicitly to prevent inflated host DEG counts
viral_features <- c("GFP", "N", "P", "M", "F", "H", "L")
deg_all <- deg_all %>%
  filter(!Gene %in% viral_features)

aud_all <- bind_rows(
  aud_prim %>% mutate(Resolution = "Primary_Lineage", Group = as.character(Group)),
  aud_sec %>% mutate(Resolution = "Secondary_Cluster", Group = as.character(Group))
)

deg_all <- deg_all %>%
  left_join(aud_all, by = c("Target_Group" = "Group", "Mapped_Audit_Contrast" = "Contrast", "Resolution"))

if (any(is.na(deg_all$Status))) {
  stop("CRITICAL ERROR: Contrast mapping failed to find matching Audit rows (NA generated).")
}

# ==============================================================================
# SECTION 1: COMPLETE PSEUDOBULK RESULT AUDIT
# ==============================================================================
message("--> Generating Post-Pseudobulk Result Audit...")

pb_audit_summary <- deg_all %>%
  group_by(Target_Group, Contrast, Resolution, Mapped_Audit_Contrast) %>%
  summarise(
    Total_Tested_Genes = n(),
    Non_NA_pvalue = sum(!is.na(pvalue)),
    Non_NA_padj = sum(!is.na(padj)),
    Pct_NA_padj = round(sum(is.na(padj)) / n() * 100, 2),
    FDR_Sig_Genes = sum(padj < 0.05, na.rm = TRUE),
    FDR_Sig_Pos_LFC = sum(padj < 0.05 & log2FoldChange > 0, na.rm = TRUE),
    FDR_Sig_Neg_LFC = sum(padj < 0.05 & log2FoldChange < 0, na.rm = TRUE),
    FDR_Sig_Large_Effect = sum(padj < 0.05 & abs(log2FoldChange) >= 1, na.rm = TRUE),
    Median_Abs_LFC_Sig = median(abs(log2FoldChange[which(padj < 0.05)]), na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(aud_all, by = c("Target_Group" = "Group", "Mapped_Audit_Contrast" = "Contrast", "Resolution")) %>%
  select(Resolution, Target_Group, Contrast, Status, Residual_DF, Support_Class, 
         Shared_Days, N_Shared_Days, N_Cond1, N_Cond2,
         Total_Tested_Genes:Median_Abs_LFC_Sig)

write.csv(pb_audit_summary, file.path(dir_results, "Pseudobulk_Result_Summary.csv"), row.names = FALSE)

# ==============================================================================
# SECTION 2: DEG CLASSIFICATION
# ==============================================================================
message("--> Separating FDR-Significant Findings...")

sig_degs <- deg_all %>%
  filter(padj < 0.05) %>%
  mutate(DEG_Category = ifelse(abs(log2FoldChange) >= 1, 
                               "Large_Effect_FDR_Significant", 
                               "FDR_Significant")) %>%
  select(Gene, Target_Group, Contrast, Resolution, log2FoldChange, lfcSE, stat, pvalue, padj, DEG_Category, Residual_DF, Support_Class)

sig_prim <- sig_degs %>% filter(Resolution == "Primary_Lineage")
sig_sec  <- sig_degs %>% filter(Resolution == "Secondary_Cluster")

write.csv(sig_sec, file.path(dir_results, "Exploratory_Cluster_Summary.csv"), row.names = FALSE)
write.csv(sig_prim, file.path(dir_results, "Primary_Lineage_FDR_Significant_Genes.csv"), row.names = FALSE)
write.csv(sig_prim %>% filter(DEG_Category == "Large_Effect_FDR_Significant"), 
          file.path(dir_results, "Primary_Lineage_LargeEffect_FDR_Genes.csv"), row.names = FALSE)

# ==============================================================================
# SECTION 3: PRIMARY DISCOVERY SET (INFECTED VS BYSTANDER)
# ==============================================================================
message("--> Summarizing Primary Infected vs Bystander Discovery Set...")

inf_prim <- sig_prim %>% filter(Contrast == "Infected_vs_Bystander")

if (nrow(inf_prim) > 0) {
  gene_overlap <- inf_prim %>%
    group_by(Gene) %>%
    summarise(
      N_Supporting_Lineages = n_distinct(Target_Group),
      Lineages_Detected = paste(sort(unique(Target_Group)), collapse = "|"),
      Direction_Consistent = length(unique(sign(log2FoldChange))) == 1,
      Max_Abs_LFC = max(abs(log2FoldChange)),
      Min_padj = min(padj),
      .groups = "drop"
    ) %>%
    arrange(desc(N_Supporting_Lineages), desc(Max_Abs_LFC))
  write.csv(gene_overlap, file.path(dir_results, "Primary_Lineage_Gene_Overlap_Summary.csv"), row.names = FALSE)
} else {
  write.csv(data.frame(Message="No FDR-significant primary genes available."), file.path(dir_results, "Primary_Lineage_Gene_Overlap_Summary.csv"), row.names = FALSE)
}

# ==============================================================================
# SECTION 4: BULK-SINGLE-CELL CROSS-DATASET CONVERGENCE
# ==============================================================================
message("--> Performing Exact Bulk/Single-Cell Convergence...")

overlap_full <- inf_prim %>%
  inner_join(bulk_metrics, by = "Gene") %>%
  mutate(
    Overlap_Type = case_when(
      sign(log2FoldChange) == sign(Mean_Bulk_LFC) ~ "Concordant",
      sign(log2FoldChange) != sign(Mean_Bulk_LFC) ~ "Discordant",
      TRUE ~ "Unknown"
    )
  ) %>%
  select(Gene, Broad_Lineage = Target_Group, Overlap_Type,
         SingleCell_LFC = log2FoldChange, SingleCell_padj = padj, SingleCell_Category = DEG_Category,
         Residual_DF, Support_Class,
         Mean_Bulk_LFC, Median_Bulk_LFC, Bulk_Direction, Direction_Consistency,
         Number_Positive_Timepoints, Number_Negative_Timepoints)

write.csv(overlap_full, file.path(dir_results, "Bulk_SingleCell_Convergence_Full.csv"), row.names = FALSE)

# Isolate True Concordant Candidates and create ranked inspection tables
concordant_genes <- overlap_full %>% filter(Overlap_Type == "Concordant")

if (nrow(concordant_genes) > 0) {
  gene_level_concordance <- concordant_genes %>%
    group_by(Gene) %>%
    summarise(
      Supporting_SC_Lineages = n_distinct(Broad_Lineage),
      Lineages = paste(sort(unique(Broad_Lineage)), collapse = "|"),
      Max_Abs_SC_LFC = max(abs(SingleCell_LFC)),
      Min_SC_padj = min(SingleCell_padj),
      Mean_Bulk_LFC = first(Mean_Bulk_LFC),
      Median_Bulk_LFC = first(Median_Bulk_LFC),
      Bulk_Direction = first(Bulk_Direction),
      Direction_Consistency = first(Direction_Consistency),
      Number_Positive_Timepoints = first(Number_Positive_Timepoints),
      Number_Negative_Timepoints = first(Number_Negative_Timepoints),
      .groups = "drop"
    ) %>%
    arrange(desc(Supporting_SC_Lineages), Min_SC_padj, desc(Max_Abs_SC_LFC))
} else {
  gene_level_concordance <- data.frame(
    Gene = character(),
    Supporting_SC_Lineages = integer(),
    Lineages = character(),
    Max_Abs_SC_LFC = numeric(),
    Min_SC_padj = numeric(),
    Mean_Bulk_LFC = numeric(),
    Median_Bulk_LFC = numeric(),
    Bulk_Direction = character(),
    Direction_Consistency = character(),
    Number_Positive_Timepoints = integer(),
    Number_Negative_Timepoints = integer()
  )
}

write.csv(gene_level_concordance, file.path(dir_results, "Bulk_SingleCell_Convergence_GeneLevel.csv"), row.names = FALSE)

# Explicitly create multi-lineage table for upcoming network input inspection
multi_lineage_concordance <- gene_level_concordance %>% filter(Supporting_SC_Lineages >= 2)
write.csv(multi_lineage_concordance, file.path(dir_results, "Bulk_SingleCell_Convergence_MultiLineage.csv"), row.names = FALSE)

# Generate Lineage-Support Sensitivity Report
sensitivity_report <- data.frame(
  Support_Threshold = c(">= 1 Lineage", ">= 2 Lineages", "3 Lineages"),
  Candidate_Count = c(
    nrow(gene_level_concordance),
    nrow(multi_lineage_concordance),
    sum(gene_level_concordance$Supporting_SC_Lineages == 3)
  )
)

write.csv(sensitivity_report, file.path(dir_results, "Lineage_Support_Sensitivity_Report.csv"), row.names = FALSE)

message("\n[ SENSITIVITY REPORT ]")
print(sensitivity_report)

# ==============================================================================
# SECTION 5: VISUALIZATIONS
# ==============================================================================
message("--> Generating Data-Driven Visualizations...")

# 1. Summary Bar Plot (FDR Significant counts)
valid_summary <- pb_audit_summary %>% filter(Resolution == "Primary_Lineage" & Status == "Success")
if(nrow(valid_summary) > 0) {
  p_summary <- valid_summary %>%
    select(Target_Group, Contrast, Upregulated = FDR_Sig_Pos_LFC, Downregulated = FDR_Sig_Neg_LFC) %>%
    pivot_longer(cols = c(Upregulated, Downregulated), names_to = "Direction", values_to = "Gene_Count") %>%
    ggplot(aes(x = Target_Group, y = Gene_Count, fill = Direction)) +
    geom_bar(stat = "identity", position = "dodge") +
    facet_wrap(~ Contrast) +
    theme_minimal() +
    scale_fill_manual(values = c("Upregulated" = "firebrick", "Downregulated" = "dodgerblue")) +
    labs(title = "Primary Lineage FDR-Significant Genes",
         subtitle = "Caution: Models subject to structural limitations (e.g., Reduced-support DF = 2)",
         x = "Primary Lineage", y = "Count of FDR-Significant Genes (padj < 0.05)")
  
  ggsave(file.path(dir_results, "Pseudobulk_DEG_Counts.png"), p_summary, width = 10, height = 6, dpi = 300)
}

# 2. Volcano Plots (Primary Lineage Infected vs Bystander)
if(nrow(inf_prim) > 0) {
  volcano_list <- list()
  for (lin in unique(deg_prim$Target_Group[deg_prim$Contrast == "Infected_vs_Bystander"])) {
    lin_raw <- deg_prim %>% filter(Target_Group == lin & Contrast == "Infected_vs_Bystander" & !is.na(padj))
    lin_audit <- aud_prim %>% filter(Group == lin & Contrast == "Inf_vs_Bys")
    
    if(nrow(lin_raw) == 0 || nrow(lin_audit) == 0) next
    
    df_stat <- lin_audit$Residual_DF[1]
    df_class <- lin_audit$Support_Class[1]
    
    lin_raw <- lin_raw %>% mutate(Sig = case_when(
      padj < 0.05 & log2FoldChange > 0 ~ "Up",
      padj < 0.05 & log2FoldChange < 0 ~ "Down",
      TRUE ~ "Not Sig"
    ))
    
    top_genes <- lin_raw %>% filter(padj < 0.05) %>% slice_max(abs(log2FoldChange), n = 10) %>% pull(Gene)
    
    p_volc <- ggplot(lin_raw, aes(x = log2FoldChange, y = -log10(padj), color = Sig)) +
      geom_point(alpha = 0.6, size = 1.5) +
      scale_color_manual(values = c("Up" = "firebrick", "Down" = "dodgerblue", "Not Sig" = "grey70")) +
      geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
      geom_text_repel(data = filter(lin_raw, Gene %in% top_genes), aes(label = Gene), size = 3, color = "black", max.overlaps = 20) +
      theme_minimal() +
      theme(legend.position = "none") +
      labs(title = paste(lin, "- Infected vs Bystander"),
           subtitle = sprintf("Support Class: %s (Residual DF = %s)", df_class, df_stat),
           x = "Log2 Fold Change", y = "-Log10(padj)")
    
    volcano_list[[lin]] <- p_volc
  }
  
  if (length(volcano_list) > 0) {
    p_grid <- wrap_plots(volcano_list, ncol = length(volcano_list))
    ggsave(file.path(dir_results, "Primary_Lineage_Volcano_Grid.png"), p_grid, width = 5 * length(volcano_list), height = 5, dpi = 300)
  }
}

# 3. Cross-Dataset Convergence Scatter & Heatmap

if (exists("concordant_genes") && nrow(concordant_genes) >= 2) {
  
  hm_data <- concordant_genes %>% 
    select(Gene, Broad_Lineage, SingleCell_LFC) %>%
    pivot_wider(names_from = Broad_Lineage, values_from = SingleCell_LFC)
  
  hm_data <- tibble::column_to_rownames(hm_data, "Gene") %>% as.matrix()
  
  hm_height <- max(1200, min(6000, nrow(hm_data) * 20))
  
  png(file.path(dir_results, "Convergent_Gene_Heatmap.png"), width = 1600, height = hm_height, res = 300)
  ht <- Heatmap(hm_data,
                name = "scRNA-seq\nLog2FC",
                cluster_rows = FALSE,     # FIX: Disables row hclust on NA values
                cluster_columns = FALSE,
                na_col = "gray85", 
                row_title = "Directionally Concordant Candidates",
                column_title = "Single-Cell Infection Effect (Primary Lineages)",
                col = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100),
                show_row_names = ifelse(nrow(hm_data) <= 150, TRUE, FALSE))
  draw(ht)
  dev.off()
}

# ==============================================================================
# SECTION 6: PROVENANCE
# ==============================================================================
writeLines(capture.output(sessionInfo()), file.path(dir_results, "SessionInfo_08.txt"))
message("==========================================================================")
message("   POST-PSEUDOBULK SUMMARY & CONVERGENCE COMPLETE")
message("==========================================================================")