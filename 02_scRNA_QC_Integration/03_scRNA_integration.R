# ==============================================================================
# SCRIPT 03: REPRESENTATION BENCHMARKING (Seurat V5 Optimized)
# ==============================================================================

library(here)
library(yaml)
library(tidyverse)
library(Seurat)
library(harmony)

config <- yaml::read_yaml(here::here("config.yml"))
set.seed(config$random_seed)

dir_processed <- here::here(config$directories$data_processed)
dir_logs      <- here::here(config$directories$logs)

seu <- readRDS(file.path(dir_processed, "GSE168775_QC_Cleaned.rds"))

# ==============================================================================
# STEP 1: Unintegrated Baseline (SCT)
# ==============================================================================
message("--> [1/4] Executing SCTransform (Unintegrated Baseline)...")

# Seurat v5 Memory Fix: Split the RNA layer virtually rather than splitting the whole object.
# This prevents RAM bloat and layer misalignment.
seu[["RNA"]] <- split(seu[["RNA"]], f = seu$Library_ID)

# SCTransform automatically detects the split layers and normalizes them independently.
seu <- SCTransform(seu, vst.flavor = "v2", verbose = FALSE)

# Compute standard PCA on the unified variable features
seu <- RunPCA(seu, assay = "SCT", npcs = 30, verbose = FALSE, reduction.name = "pca_unintegrated")

# Dimensionality Assertions
stopifnot("CRITICAL: PCA embedding row count mismatch!" = nrow(Embeddings(seu, "pca_unintegrated")) == ncol(seu))
stopifnot("CRITICAL: PCA embeddings contain NA values!" = !anyNA(Embeddings(seu, "pca_unintegrated")))

seu <- RunUMAP(seu, dims = 1:30, reduction = "pca_unintegrated", reduction.name = "umap_unintegrated", verbose = FALSE)

# ==============================================================================
# STEP 2: RPCA Representation
# ==============================================================================
message("--> [2/4] Executing RPCA Representation...")

# Seurat v5 Native Integration: Bypasses FindIntegrationAnchors to prevent std::bad_alloc (RAM Crash)
seu <- IntegrateLayers(
  object = seu,
  method = RPCAIntegration,
  normalization.method = "SCT",
  orig.reduction = "pca_unintegrated",
  new.reduction = "pca_rpca",
  verbose = FALSE
)

seu <- RunUMAP(seu, dims = 1:30, reduction = "pca_rpca", reduction.name = "umap_rpca", verbose = FALSE)

# ==============================================================================
# STEP 3: Harmony Representation
# ==============================================================================
message("--> [3/4] Executing Harmony Representation...")

seu <- IntegrateLayers(
  object = seu,
  method = HarmonyIntegration,
  normalization.method = "SCT",
  orig.reduction = "pca_unintegrated",
  new.reduction = "harmony",
  assay = "SCT",
  theta = 2, # Explicit Harmony parameter
  verbose = FALSE
)

seu <- RunUMAP(seu, dims = 1:30, reduction = "harmony", reduction.name = "umap_harmony", verbose = FALSE)

# ==============================================================================
# STEP 4: Assemble Object & Cleanup
# ==============================================================================
message("--> [4/4] Assembling Master Object...")

# Neutral terminology
seu$Condition_Class <- ifelse(seu$Day == "Day3" & seu$Condition == "Mock", "Baseline", "Exposed")

# SEURAT V5 FIX: Rejoin ONLY the RNA layers. 
# SCT layers do not require/support JoinLayers() in v5.
seu[["RNA"]] <- JoinLayers(seu[["RNA"]])
DefaultAssay(seu) <- "SCT"

saveRDS(seu, file.path(dir_processed, "GSE168775_Integrated_Benchmarked.rds"))

time_stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
writeLines(capture.output(sessionInfo()), file.path(dir_logs, paste0("SessionInfo_03_", time_stamp, ".txt")))

message("SUCCESS: Representation benchmarking complete. RAM safely preserved. Ready for Script 04 (Reviewer Diagnostics).")