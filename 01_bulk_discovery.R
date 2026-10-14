# ==============================================================================
# SCRIPT 01: INDEPENDENT BULK LONGITUDINAL DISCOVERY (GSE229685)
# ==============================================================================

library(here)
library(yaml)
library(tidyverse)
library(DESeq2)
library(apeglm)

config <- yaml::read_yaml(here::here("config.yml"))
set.seed(config$random_seed)

dir_processed <- here::here(config$directories$data_processed)
dir_results   <- here::here(config$directories$results_bulk)

# ==============================================================================
# STEP 1: Metadata Formatting & Matrix Cleaning
# ==============================================================================
message("--> [1/6] Ingesting and Cleaning Data...")

counts   <- read.csv(file.path(dir_processed, "GSE229685_counts.csv"), row.names = 1)
metadata <- read.csv(file.path(dir_processed, "GSE229685_metadata.csv"), row.names = 1)

# Clean Metadata (Whitespace trimming & strict factoring)
metadata$subject <- as.factor(trimws(metadata$subject))
metadata$time    <- factor(trimws(metadata$time), levels = c("0", "12", "24", "48", "72", "120"))

# Clean Count Matrix (Handle NAs & Convert Floats to Integers for DESeq2)
counts[is.na(counts)] <- 0
counts <- round(as.matrix(counts))
mode(counts) <- "integer"

stopifnot("CRITICAL: Matrix columns and Metadata rows do not align!" = all(colnames(counts) == rownames(metadata)))

# ==============================================================================
# STEP 2: Exact Donor x Time Structure Audit
# ==============================================================================
actual_structure <- table(metadata$subject, metadata$time)

expected_structure <- matrix(
  c(1, 1, 1, 1, 1, 1,  # NCF1
    1, 1, 1, 1, 1, 1,  # NCF2
    1, 1, 0, 1, 1, 1,  # NCF3 (Missing 24h)
    1, 1, 1, 1, 1, 1), # NCF4
  nrow = 4, byrow = TRUE,
  dimnames = list(c("NCF1", "NCF2", "NCF3", "NCF4"), c("0", "12", "24", "48", "72", "120"))
)

# Align subject naming formats dynamically
if (all(rownames(actual_structure) %in% c("1", "2", "3", "4"))) {
  rownames(expected_structure) <- c("1", "2", "3", "4")
}

# Explicitly align dimnames before positional comparison
actual_structure <- actual_structure[rownames(expected_structure), colnames(expected_structure)]

if (!all(as.matrix(actual_structure) == expected_structure)) {
  stop("CRITICAL ERROR: Donor x time structure does not match the locked unbalanced design.")
}

# ==============================================================================
# STEP 3: DESeq2 Construction & Pre-filtering
# ==============================================================================
message("--> [2/6] Building DESeq2 Object...")

dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData   = metadata,
  design    = as.formula(config$bulk_params$model_design)
)

keep <- rowSums(counts(dds) >= config$bulk_params$min_count_floor) >= config$bulk_params$min_count_sample_threshold
dds  <- dds[keep, ]
message(sprintf("    Retained %d features passing explicit pre-filter.", nrow(dds)))

# ==============================================================================
# STEP 4: Global Temporal Dynamics (LRT)
# ==============================================================================
message("--> [3/6] Executing Likelihood Ratio Test (LRT)...")

dds_lrt <- DESeq(dds, test = "LRT", reduced = as.formula(config$bulk_params$reduced_design))
res_lrt <- results(dds_lrt, alpha = config$bulk_params$discovery_fdr)
lrt_sig_genes <- rownames(subset(res_lrt, padj < config$bulk_params$discovery_fdr))

# EXPORT: Save full LRT results for Script 04 (Reviewer Diagnostics)
write.csv(as.data.frame(res_lrt) %>% rownames_to_column("Gene"), 
          file.path(dir_results, "GSE229685_LRT_full_results.csv"), row.names = FALSE)

# ==============================================================================
# STEP 5: Wald Time-vs-Baseline Contrasts & apeglm Shrinkage
# ==============================================================================
message("--> [4/6] Computing Wald Contrasts & apeglm Shrinkage...")

dds_wald <- DESeq(dds, test = "Wald")
available_coefs <- resultsNames(dds_wald) # RESTORED: Validate coefficients
timepoints <- c("12", "24", "48", "72", "120")
expected_coefs <- paste0("time_", timepoints, "_vs_0")

if (!all(expected_coefs %in% available_coefs)) {
  stop("CRITICAL ERROR: Expected time coefficients missing. Available: ", paste(available_coefs, collapse = ", "))
}

contrast_results <- list()

for (tp in timepoints) {
  coef_name <- paste0("time_", tp, "_vs_0")
  res_w <- results(dds_wald, name = coef_name, alpha = config$bulk_params$discovery_fdr)
  res_shrunk <- lfcShrink(dds_wald, coef = coef_name, type = "apeglm", res = res_w, quiet = TRUE)
  
  contrast_results[[tp]] <- as.data.frame(res_shrunk) %>%
    rownames_to_column("Gene") %>%
    select(Gene, padj, log2FoldChange) %>%
    rename_with(~paste0(.x, "_", tp, "h"), c(padj, log2FoldChange))
}

master_res <- purrr::reduce(contrast_results, full_join, by = "Gene")

# ==============================================================================
# STEP 6: Parameterized Persistence & Directional Filtering
# ==============================================================================
message("--> [5/6] Applying Operational Persistence Rules...")

final_candidates <- master_res %>%
  filter(Gene %in% lrt_sig_genes) %>%
  rowwise() %>%
  mutate(
    is_persistent = {
      p_vals <- c_across(starts_with("padj_"))
      lfcs   <- c_across(starts_with("log2FoldChange_"))
      sig_indices <- which(!is.na(p_vals) & p_vals < config$bulk_params$discovery_fdr)
      
      # Parameterized sensitivity check
      if (length(sig_indices) < config$bulk_params$persistence_required) { 
        FALSE 
      } else {
        sig_lfcs <- lfcs[sig_indices]
        all(sig_lfcs > 0) || all(sig_lfcs < 0) 
      }
    }
  ) %>% ungroup() %>% filter(is_persistent == TRUE) %>% select(-is_persistent)

# Dynamic file naming based on config parameter
out_filename <- sprintf("Persistent_Bulk_Gene_Set_%dof5.csv", config$bulk_params$persistence_required)
write.csv(final_candidates, file.path(dir_results, out_filename), row.names = FALSE)

message(sprintf("SUCCESS: Discovered %d target-blind persistent genes.", nrow(final_candidates)))