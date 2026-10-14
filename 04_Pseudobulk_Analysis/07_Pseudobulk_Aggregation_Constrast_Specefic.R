# ==============================================================================
# SCRIPT 07: PSEUDOBULK AGGREGATION & CONTRAST-SPECIFIC DESEQ2 (STRICT)
# ==============================================================================

library(here)
library(Seurat)
library(DESeq2)
library(tidyverse)

dir_processed <- here::here("data", "processed")
dir_results   <- here::here("results", "04_pseudobulk_de")
if (!dir.exists(dir_results)) dir.create(dir_results, recursive = TRUE)

message("--> Loading Annotated Object...")
seu <- readRDS(file.path(dir_processed, "GSE168775_Annotated_Final.rds"))

# Set working parameters
MIN_CELLS <- 10
PRIMARY_LINEAGES <- c("Basal", "Secretory", "Ciliated") # Explicitly excluding 'Unresolved'

# ==============================================================================
# 1. CORE AGGREGATION FUNCTION (EXPLICIT EXPERIMENTAL UNIT IDs)
# ==============================================================================
create_pseudobulk <- function(seu_obj, group_var, min_cells) {
  
  # 1. Build Metadata Map and Filter
  meta <- seu_obj@meta.data %>%
    group_by(!!sym(group_var), Library_ID, Condition, Day) %>%
    summarise(Cell_Count = n(), .groups = "drop") %>%
    filter(Cell_Count >= min_cells) %>%
    # Use dashes instead of underscores for the separator
    mutate(PB_ID = paste("PB", !!sym(group_var), Library_ID, Condition, Day, sep = "---")) %>%
    # Force spaces, slashes, AND underscores to become dashes
    mutate(PB_ID = gsub("[ _/\\\\]", "-", PB_ID))
  
  # Absolute safeguard against duplicated PB_IDs
  if(any(duplicated(meta$PB_ID))) {
    stop("CRITICAL ERROR: Duplicate Pseudobulk IDs generated. Check grouping variables.")
  }
  
  # 2. Assign PB_ID to cells and subset Seurat object to valid cells
  seu_obj$PB_ID <- paste("PB", seu_obj@meta.data[[group_var]], seu_obj$Library_ID, seu_obj$Condition, seu_obj$Day, sep = "---")
  seu_obj$PB_ID <- gsub("[ _/\\\\]", "-", seu_obj$PB_ID) 
  
  valid_cells <- rownames(seu_obj@meta.data)[seu_obj$PB_ID %in% meta$PB_ID]
  seu_sub <- subset(seu_obj, cells = valid_cells)
  
  # 3. Aggregate sparse matrix
  pb_counts <- AggregateExpression(seu_sub, group.by = "PB_ID", assays = "RNA", slot = "counts")$RNA
  
  # 4. Strict Alignment Check & Reordering
  meta_final <- as.data.frame(meta) %>% column_to_rownames("PB_ID")
  
  if (!setequal(rownames(meta_final), colnames(pb_counts))) {
    stop("CRITICAL ERROR: Pseudobulk metadata IDs and count matrix IDs do not match exactly.")
  }
  
  # Explicitly reorder counts to match metadata row order
  pb_counts <- pb_counts[, rownames(meta_final), drop = FALSE]
  
  # Return structured list
  return(list(
    counts = pb_counts,
    meta = meta_final
  ))
}

# ==============================================================================
# 2. DESEQ2 WRAPPER & STRICT MATRIX VALIDATION FUNCTION
# ==============================================================================
run_contrast_de <- function(pb_data, target_group, contrast_name, cond_filter, day_filter = NULL, required_shared_days = NULL) {
  
  group_col <- names(pb_data$meta)[1]
  
  # Subset metadata carefully preserving rownames
  sub_meta <- pb_data$meta %>%
    rownames_to_column("PB_ID") %>%
    filter(!!sym(group_col) == target_group) %>%
    filter(Condition %in% cond_filter)
  
  if(!is.null(day_filter)) {
    sub_meta <- sub_meta %>% filter(Day %in% day_filter)
  }
  
  sub_meta <- sub_meta %>% column_to_rownames("PB_ID")
  
  # Extract condition counts safely for reporting
  cond_counts <- table(factor(sub_meta$Condition, levels = cond_filter))
  n_cond1 <- unname(cond_counts[1])
  n_cond2 <- unname(cond_counts[2])
  
  # Base return template for exclusions
  ret_template <- list(results = NULL, model = NA, res_df = NA, support_class = NA, shared_days = NA, n_shared_days = NA, n_cond1 = n_cond1, n_cond2 = n_cond2)
  
  # ----------------------------------------------------------------------------
  # PRE-MODEL VALIDATION GATES
  # ----------------------------------------------------------------------------
  
  # Gate 1: Require at least 2 pseudobulk samples per condition
  if (any(cond_counts < 2)) {
    ret_template$status <- "Excluded: Insufficient Replication (<2 libraries per condition)"
    return(ret_template)
  }
  
  # Extract shared days for Gate 2 and reporting
  day_cond1 <- unique(as.character(sub_meta$Day[sub_meta$Condition == cond_filter[1]]))
  day_cond2 <- unique(as.character(sub_meta$Day[sub_meta$Condition == cond_filter[2]]))
  shared_days <- intersect(day_cond1, day_cond2)
  
  ret_template$shared_days <- paste(shared_days, collapse = "|")
  ret_template$n_shared_days <- length(shared_days)
  
  # Gate 2: Check for shared temporal overlap
  if (length(shared_days) == 0) {
    ret_template$status <- "Excluded: No Shared Timepoints"
    return(ret_template)
  }
  
  # Gate 2.5: Explicitly check required shared timepoints if specified
  if (!is.null(required_shared_days)) {
    missing_shared_days <- setdiff(required_shared_days, shared_days)
    if (length(missing_shared_days) > 0) {
      ret_template$status <- paste0("Excluded: Missing Required Shared Timepoint(s): ", paste(missing_shared_days, collapse = ", "))
      return(ret_template)
    }
  }
  
  # Gate 3: Check minimal temporal structure for longitudinal model
  if (length(unique(sub_meta$Day)) < 2) {
    ret_template$status <- "Excluded: Insufficient Timepoints for Longitudinal Model"
    return(ret_template)
  }
  
  # Apply factors
  sub_meta$Condition <- factor(sub_meta$Condition, levels = cond_filter)
  sub_meta$Day <- factor(
    sub_meta$Day,
    levels = c("Day3", "Day7", "Day14")
  )
  
  # Gate 4 & 5: Strict Matrix Math Evaluator for ~ Day + Condition
  X_day <- model.matrix(~ Day + Condition, data = sub_meta)
  rank_X <- qr(X_day)$rank
  cols_X <- ncol(X_day)
  residual_df <- nrow(X_day) - rank_X
  
  # Immediately calculate support class so it is logged even if the model is excluded
  support_class <- case_when(
    residual_df >= 3 ~ "Well-supported",
    residual_df == 2 ~ "Reduced-support",
    residual_df == 1 ~ "Fragile",
    TRUE ~ "Unknown"
  )
  
  ret_template$res_df <- residual_df
  ret_template$support_class <- support_class
  
  if (rank_X < cols_X) {
    ret_template$status <- "Excluded: Matrix Rank Deficient"
    return(ret_template)
  }
  if (residual_df <= 0) {
    ret_template$status <- "Excluded: Zero/Negative Residual DF"
    return(ret_template)
  }
  
  # ----------------------------------------------------------------------------
  # DESEQ2 EXECUTION (STRICTLY ~ Day + Condition)
  # ----------------------------------------------------------------------------
  model_formula <- "~ Day + Condition"
  sub_counts <- pb_data$counts[, rownames(sub_meta), drop = FALSE]
  
  dds <- DESeqDataSetFromMatrix(countData = sub_counts, 
                                colData = sub_meta, 
                                design = as.formula(model_formula))
  
  # Pre-filter low count genes
  keep <- rowSums(counts(dds) >= 10) >= min(table(sub_meta$Condition))
  dds <- dds[keep,]
  
  dds <- DESeq(dds, quiet = TRUE)
  
  # Extract contrast
  res <- results(dds, contrast = c("Condition", cond_filter[2], cond_filter[1]), alpha = 0.05)
  res_df <- as.data.frame(res) %>% 
    rownames_to_column("Gene") %>%
    mutate(Target_Group = target_group, Contrast = contrast_name)
  
  return(list(
    status = "Success",
    results = res_df,
    model = model_formula,
    res_df = residual_df,
    support_class = support_class,
    shared_days = paste(shared_days, collapse = "|"),
    n_shared_days = length(shared_days),
    n_cond1 = n_cond1,
    n_cond2 = n_cond2
  ))
}

# ==============================================================================
# 3. EXECUTE PRIMARY ANALYSIS: BROAD LINEAGES
# ==============================================================================
message("--> Generating Primary Lineage Pseudobulks (N >= 10)...")
pb_lineage <- create_pseudobulk(seu, "Broad_Lineage", MIN_CELLS)

all_lineage_de <- list()
lineage_audit <- list()

for(lin in PRIMARY_LINEAGES) {
  message(sprintf("    Testing Lineage: %s", lin))
  
  # 3A. Bystander vs Mock (Day 3 & Day 14 ONLY, REQUIRED)
  res_bys_mock <- run_contrast_de(pb_lineage, lin, "Bystander_vs_Mock", 
                                  cond_filter = c("Mock", "Bystander"), 
                                  day_filter = c("Day3", "Day14"),
                                  required_shared_days = c("Day3", "Day14"))
  
  if(!is.null(res_bys_mock$results)) all_lineage_de[[paste(lin, "Bys_Mock")]] <- res_bys_mock$results
  lineage_audit[[paste(lin, "Bys_Mock")]] <- data.frame(
    Group = lin, Contrast = "Bys_vs_Mock", Status = res_bys_mock$status, Model = res_bys_mock$model,
    Residual_DF = res_bys_mock$res_df, Support_Class = res_bys_mock$support_class,
    Shared_Days = res_bys_mock$shared_days, N_Shared_Days = res_bys_mock$n_shared_days,
    N_Cond1 = res_bys_mock$n_cond1, N_Cond2 = res_bys_mock$n_cond2
  )
  
  # 3B. Infected vs Bystander (All Days REQUIRED)
  res_inf_bys <- run_contrast_de(pb_lineage, lin, "Infected_vs_Bystander", 
                                 cond_filter = c("Bystander", "Infected"), 
                                 day_filter = NULL, 
                                 required_shared_days = c("Day3", "Day7", "Day14"))
  
  if(!is.null(res_inf_bys$results)) all_lineage_de[[paste(lin, "Inf_Bys")]] <- res_inf_bys$results
  lineage_audit[[paste(lin, "Inf_Bys")]] <- data.frame(
    Group = lin, Contrast = "Inf_vs_Bys", Status = res_inf_bys$status, Model = res_inf_bys$model,
    Residual_DF = res_inf_bys$res_df, Support_Class = res_inf_bys$support_class,
    Shared_Days = res_inf_bys$shared_days, N_Shared_Days = res_inf_bys$n_shared_days,
    N_Cond1 = res_inf_bys$n_cond1, N_Cond2 = res_inf_bys$n_cond2
  )
}

# Export Primary Results
write.csv(bind_rows(all_lineage_de), file.path(dir_results, "DEGs_Primary_Broad_Lineages.csv"), row.names = FALSE)
write.csv(bind_rows(lineage_audit), file.path(dir_results, "Audit_Primary_Model_Usage.csv"), row.names = FALSE)

# ==============================================================================
# 4. EXECUTE SECONDARY ANALYSIS: CLUSTERS (EXPLORATORY)
# ==============================================================================
message("--> Generating Secondary Cluster Pseudobulks (N >= 10)...")
pb_cluster <- create_pseudobulk(seu, "seurat_clusters", MIN_CELLS)
all_clusters <- unique(pb_cluster$meta$seurat_clusters)

all_cluster_de <- list()
cluster_audit <- list()

for(clu in all_clusters) {
  
  # 4A. Bystander vs Mock (Day 3 & Day 14 ONLY, REQUIRED)
  res_bys_mock <- run_contrast_de(pb_cluster, clu, "Bystander_vs_Mock", 
                                  cond_filter = c("Mock", "Bystander"), 
                                  day_filter = c("Day3", "Day14"),
                                  required_shared_days = c("Day3", "Day14"))
  
  if(!is.null(res_bys_mock$results)) all_cluster_de[[paste(clu, "Bys_Mock")]] <- res_bys_mock$results
  cluster_audit[[paste(clu, "Bys_Mock")]] <- data.frame(
    Group = clu, Contrast = "Bys_vs_Mock", Status = res_bys_mock$status, Model = res_bys_mock$model,
    Residual_DF = res_bys_mock$res_df, Support_Class = res_bys_mock$support_class,
    Shared_Days = res_bys_mock$shared_days, N_Shared_Days = res_bys_mock$n_shared_days,
    N_Cond1 = res_bys_mock$n_cond1, N_Cond2 = res_bys_mock$n_cond2
  )
  
  # 4B. Infected vs Bystander (All Days REQUIRED)
  res_inf_bys <- run_contrast_de(pb_cluster, clu, "Infected_vs_Bystander", 
                                 cond_filter = c("Bystander", "Infected"), 
                                 day_filter = NULL, 
                                 required_shared_days = c("Day3", "Day7", "Day14"))
  
  if(!is.null(res_inf_bys$results)) all_cluster_de[[paste(clu, "Inf_Bys")]] <- res_inf_bys$results
  cluster_audit[[paste(clu, "Inf_Bys")]] <- data.frame(
    Group = clu, Contrast = "Inf_vs_Bys", Status = res_inf_bys$status, Model = res_inf_bys$model,
    Residual_DF = res_inf_bys$res_df, Support_Class = res_inf_bys$support_class,
    Shared_Days = res_inf_bys$shared_days, N_Shared_Days = res_inf_bys$n_shared_days,
    N_Cond1 = res_inf_bys$n_cond1, N_Cond2 = res_inf_bys$n_cond2
  )
}

# Export Secondary Results
write.csv(bind_rows(all_cluster_de), file.path(dir_results, "DEGs_Secondary_Clusters.csv"), row.names = FALSE)
write.csv(bind_rows(cluster_audit), file.path(dir_results, "Audit_Secondary_Model_Usage.csv"), row.names = FALSE)

message("==========================================================================")
message("   PSEUDOBULK DESEQ2 ANALYSIS COMPLETE")
message("==========================================================================")