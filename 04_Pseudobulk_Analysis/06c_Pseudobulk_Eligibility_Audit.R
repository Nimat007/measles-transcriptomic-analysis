# ==============================================================================
# SCRIPT 06c: PSEUDOBULK ELIGIBILITY & MODEL MATRIX AUDIT
# ==============================================================================

library(here)
library(Seurat)
library(tidyverse)

dir_processed <- here::here("data", "processed")
dir_results   <- here::here("results", "03_pseudobulk_audit")
if (!dir.exists(dir_results)) dir.create(dir_results, recursive = TRUE)

message("--> Loading Annotated Object...")
# UPDATED: Now points to the corrected Final object
seu <- readRDS(file.path(dir_processed, "GSE168775_Annotated_Final.rds"))

# ------------------------------------------------------------------------------
# 1. Genuine GFP Leakage Guard Check
# ------------------------------------------------------------------------------
message("--> Auditing Existing FACS Leakage Filters...")

# Broad search for logging candidate metadata
keywords <- c("GFP", "viral", "virus", "measles")
candidate_cols <- grep(paste(keywords, collapse="|"), colnames(seu@meta.data), ignore.case = TRUE, value = TRUE)
message("    Candidate leakage/viral columns found: ", if(length(candidate_cols) > 0) paste(candidate_cols, collapse=", ") else "None")

# Strict explicit check for validated categorical/continuous columns
preferred_cols <- c("GFP", "Strict_GFP", "Viral_Status")
actual_gfp_col <- intersect(preferred_cols, colnames(seu@meta.data))[1]

if (!is.na(actual_gfp_col)) {
  message("    Executing leakage audit using explicit column: ", actual_gfp_col)
  gfp_audit <- seu@meta.data %>%
    dplyr::group_by(Condition, .data[[actual_gfp_col]]) %>%
    dplyr::summarise(Cell_Count = dplyr::n(), .groups = "drop")
  
  write.csv(gfp_audit, file.path(dir_results, "Pre_DE_GFP_Leakage_Audit.csv"), row.names = FALSE)
} else {
  message("    NOTE: No explicitly validated categorical GFP/viral metadata column found.")
  leakage_fail_msg <- data.frame(
    Diagnostic = "GFP Leakage Audit",
    Status = "Failed to Execute",
    Reason = "Independent leakage auditing is unavailable in the current object because no validated metadata column (e.g., GFP, Strict_GFP, Viral_Status) was found."
  )
  write.csv(leakage_fail_msg, file.path(dir_results, "Pre_DE_GFP_Leakage_Audit_UNAVAILABLE.csv"), row.names = FALSE)
}

# ------------------------------------------------------------------------------
# 2. Extract Base Composition Table
# ------------------------------------------------------------------------------
base_comp <- seu@meta.data %>%
  group_by(Library_ID, Condition, Day, Broad_Lineage, seurat_clusters) %>%
  summarise(Cell_Count = n(), .groups = "drop")

# ------------------------------------------------------------------------------
# 3. Dynamic Threshold & Matrix Estimability Function
# ------------------------------------------------------------------------------
audit_pseudobulk_eligibility <- function(df, group_col, thresholds = c(10, 20, 30)) {
  
  results_list <- list()
  expected_conditions <- c("Mock", "Bystander", "Infected")
  expected_days <- c("Day3", "Day7", "Day14")
  
  # Ensure all original groups are evaluated at every threshold so missingness is logged
  all_target_groups <- unique(df[[group_col]])
  
  for (n_thresh in thresholds) {
    retained_libs <- df %>%
      group_by(!!sym(group_col), Library_ID, Condition, Day) %>%
      summarise(Total_Cells = sum(Cell_Count), .groups = "drop") %>%
      filter(Total_Cells >= n_thresh)
    
    for (tg in all_target_groups) {
      tg_meta <- retained_libs %>% filter(!!sym(group_col) == tg) 
      
      # Explicit zero-retention logic
      if (nrow(tg_meta) == 0) {
        results_list[[paste(tg, n_thresh, sep="_")]] <- data.frame(
          Aggregation_Level = group_col,
          Target_Group = tg,
          N_Threshold = n_thresh,
          Missing_Condition = paste(expected_conditions, collapse="|"),
          Missing_Day = paste(expected_days, collapse="|"),
          Model_Cond_Cols = NA, Model_Cond_Rank = NA, Model_Cond_FullRank = FALSE, Model_Cond_ResDF = NA, Condition_Matrix_Estimable = FALSE,
          Model_DayCond_Cols = NA, Model_DayCond_Rank = NA, Model_DayCond_FullRank = FALSE, Model_DayCond_ResDF = NA, DayCond_Matrix_Estimable = FALSE,
          Bys_vs_Mock_Both_Conditions_Present = FALSE, Bys_vs_Mock_Obs_Mock = 0, Bys_vs_Mock_Obs_Bys = 0,
          Bys_vs_Mock_Shared_Days = "None", Bys_vs_Mock_Shared_Days_N = 0, Bys_vs_Mock_Contrast_Data_Status = "No libraries retained at this threshold",
          Inf_vs_Bys_Both_Conditions_Present = FALSE, Inf_vs_Bys_Obs_Bys = 0, Inf_vs_Bys_Obs_Inf = 0,
          Inf_vs_Bys_Shared_Days = "None", Inf_vs_Bys_Shared_Days_N = 0, Inf_vs_Bys_Contrast_Data_Status = "No libraries retained at this threshold"
        )
        next
      }
      
      # If libraries exist, proceed with full logic
      tg_meta <- tg_meta %>%
        mutate(Condition = factor(Condition, levels = expected_conditions),
               Day = factor(Day, levels = expected_days))
      
      n_mock <- sum(tg_meta$Condition == "Mock")
      n_bys  <- sum(tg_meta$Condition == "Bystander")
      n_inf  <- sum(tg_meta$Condition == "Infected")
      
      days_mock <- unique(as.character(tg_meta$Day[tg_meta$Condition == "Mock"]))
      days_bys  <- unique(as.character(tg_meta$Day[tg_meta$Condition == "Bystander"]))
      days_inf  <- unique(as.character(tg_meta$Day[tg_meta$Condition == "Infected"]))
      
      shared_bys_mock <- intersect(days_mock, days_bys)
      shared_inf_bys  <- intersect(days_bys, days_inf)
      
      missing_conds <- paste(setdiff(expected_conditions, unique(as.character(tg_meta$Condition))), collapse = "|")
      missing_days  <- paste(setdiff(expected_days, unique(as.character(tg_meta$Day))), collapse = "|")
      if(missing_conds == "") missing_conds <- "None"
      if(missing_days == "") missing_days <- "None"
      
      # 1. Model: ~ Condition
      cond_cols <- NA; cond_rank <- NA; cond_df <- NA
      cond_full_rank <- FALSE; cond_estimable <- FALSE
      
      if(length(unique(tg_meta$Condition)) > 1) {
        try({
          X1 <- model.matrix(~ Condition, data = tg_meta)
          cond_cols <- ncol(X1)
          cond_rank <- qr(X1)$rank
          cond_df <- nrow(X1) - cond_rank
          cond_full_rank <- identical(cond_rank, cond_cols)
          cond_estimable <- cond_full_rank && (cond_df > 0)
        }, silent = TRUE)
      }
      
      # 2. Model: ~ Day + Condition
      daycond_cols <- NA; daycond_rank <- NA; daycond_df <- NA
      daycond_full_rank <- FALSE; daycond_estimable <- FALSE
      
      if(length(unique(tg_meta$Condition)) > 1 && length(unique(tg_meta$Day)) > 1) {
        try({
          X2 <- model.matrix(~ Day + Condition, data = tg_meta)
          daycond_cols <- ncol(X2)
          daycond_rank <- qr(X2)$rank
          daycond_df <- nrow(X2) - daycond_rank
          daycond_full_rank <- identical(daycond_rank, daycond_cols)
          daycond_estimable <- daycond_full_rank && (daycond_df > 0)
        }, silent = TRUE)
      }
      
      # Contrast Data Structure & Interpretability Status
      bys_mock_present <- (n_mock > 0) && (n_bys > 0)
      if (!bys_mock_present) {
        bys_mock_status <- "Missing required condition"
      } else if (length(shared_bys_mock) == 0) {
        bys_mock_status <- "Conditions present but lack shared timepoints"
      } else {
        bys_mock_status <- "Shared timepoint(s) available (Review matrix estimability)"
      }
      
      inf_bys_present <- (n_bys > 0) && (n_inf > 0)
      if (!inf_bys_present) {
        inf_bys_status <- "Missing required condition"
      } else if (length(shared_inf_bys) == 0) {
        inf_bys_status <- "Conditions present but lack shared timepoints"
      } else {
        inf_bys_status <- "Shared timepoint(s) available (Review matrix estimability)"
      }
      
      results_list[[paste(tg, n_thresh, sep="_")]] <- data.frame(
        Aggregation_Level = group_col,
        Target_Group = tg,
        N_Threshold = n_thresh,
        
        Missing_Condition = missing_conds,
        Missing_Day = missing_days,
        
        Model_Cond_Cols = cond_cols,
        Model_Cond_Rank = cond_rank,
        Model_Cond_FullRank = cond_full_rank,
        Model_Cond_ResDF = cond_df,
        Condition_Matrix_Estimable = cond_estimable,
        
        Model_DayCond_Cols = daycond_cols,
        Model_DayCond_Rank = daycond_rank,
        Model_DayCond_FullRank = daycond_full_rank,
        Model_DayCond_ResDF = daycond_df,
        DayCond_Matrix_Estimable = daycond_estimable,
        
        Bys_vs_Mock_Both_Conditions_Present = bys_mock_present,
        Bys_vs_Mock_Obs_Mock = n_mock,
        Bys_vs_Mock_Obs_Bys = n_bys,
        Bys_vs_Mock_Shared_Days = ifelse(length(shared_bys_mock) > 0, paste(shared_bys_mock, collapse="|"), "None"),
        Bys_vs_Mock_Shared_Days_N = length(shared_bys_mock),
        Bys_vs_Mock_Contrast_Data_Status = bys_mock_status,
        
        Inf_vs_Bys_Both_Conditions_Present = inf_bys_present,
        Inf_vs_Bys_Obs_Bys = n_bys,
        Inf_vs_Bys_Obs_Inf = n_inf,
        Inf_vs_Bys_Shared_Days = ifelse(length(shared_inf_bys) > 0, paste(shared_inf_bys, collapse="|"), "None"),
        Inf_vs_Bys_Shared_Days_N = length(shared_inf_bys),
        Inf_vs_Bys_Contrast_Data_Status = inf_bys_status
      )
    }
  }
  return(bind_rows(results_list))
}

# ------------------------------------------------------------------------------
# 4. Execute Audits & Export
# ------------------------------------------------------------------------------
message("--> Evaluating thresholds (N=10, 20, 30) for individual Clusters...")
cluster_audit <- audit_pseudobulk_eligibility(base_comp, "seurat_clusters", c(10, 20, 30))
write.csv(cluster_audit, file.path(dir_results, "Threshold_Audit_Clusters.csv"), row.names = FALSE)

message("--> Evaluating thresholds (N=10, 20, 30) for Broad Lineages...")
lineage_audit <- audit_pseudobulk_eligibility(base_comp, "Broad_Lineage", c(10, 20, 30))
write.csv(lineage_audit, file.path(dir_results, "Threshold_Audit_Broad_Lineages.csv"), row.names = FALSE)

# Generate an explicit Library Exclusion table for N=10 (Baseline chosen for Script 07)
excluded_libs <- base_comp %>%
  group_by(Broad_Lineage, Library_ID, Condition, Day) %>%
  summarise(Total_Cells = sum(Cell_Count), .groups = "drop") %>%
  filter(Total_Cells < 10) %>%
  arrange(Broad_Lineage, Total_Cells)
write.csv(excluded_libs, file.path(dir_results, "Library_Exclusions_N10_Lineage.csv"), row.names = FALSE)

message("==========================================================================")
message("   PSEUDOBULK ELIGIBILITY & MODEL MATRIX AUDIT COMPLETE")
message("==========================================================================")