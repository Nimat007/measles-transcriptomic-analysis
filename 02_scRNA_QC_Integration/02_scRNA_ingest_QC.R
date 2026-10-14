# ==============================================================================
# SCRIPT 02: SINGLE-CELL INGESTION, VIRAL AUDIT & ADAPTIVE QC (GSE168775)
# ==============================================================================

library(here)
library(yaml)
library(tidyverse)
library(Seurat)
library(scDblFinder)
library(SingleCellExperiment)

config <- yaml::read_yaml(here::here("config.yml"))
set.seed(config$random_seed)

dir_raw       <- here::here(config$directories$data_raw, config$scRNA_params$dataset)
dir_processed <- here::here(config$directories$data_processed)
dir_results   <- here::here(config$directories$results_sc)
dir_logs      <- here::here(config$directories$logs)

if (!dir.exists(dir_results)) dir.create(dir_results, recursive = TRUE)

# ==============================================================================
# STEP 1: Strict Metadata Locking & Ingestion
# ==============================================================================
message("--> [1/6] Ingesting independent sequencing libraries...")

# GSM -> Condition mapping
locked_metadata <- data.frame(
  Library_ID = c("GSM5169176", "GSM5169175", "GSM5169177", 
                 "GSM5169181", "GSM5169179", 
                 "GSM5169182", "GSM5169180", "GSM5169178"),
  Day        = c("Day3", "Day3", "Day3", "Day7", "Day7", "Day14", "Day14", "Day14"),
  Condition  = c("Mock", "Bystander", "Infected", "Bystander", "Infected", "Mock", "Bystander", "Infected"),
  stringsAsFactors = FALSE
)

h5_files <- list.files(dir_raw, pattern = "\\.h5$", full.names = TRUE)
if (length(h5_files) != 8) stop("CRITICAL ERROR: Expected exactly 8 independent library .h5 files.")

library_list <- list()
for (f in h5_files) {
  lib_id <- sub("_.*", "", basename(f))
  if (!lib_id %in% locked_metadata$Library_ID) stop(paste("CRITICAL ERROR: Unexpected Library ID:", lib_id))
  
  meta <- locked_metadata %>% filter(Library_ID == lib_id)
  mat  <- Read10X_h5(f)
  if (is.list(mat)) mat <- mat$`Gene Expression` 
  
  seu <- CreateSeuratObject(counts = mat, project = lib_id)
  seu$Library_ID <- lib_id
  seu$Day        <- meta$Day
  seu$Condition  <- meta$Condition
  library_list[[lib_id]] <- seu
}

# ==============================================================================
# STEP 2: Viral Dictionary & Organelle Audit
# ==============================================================================
message("--> [2/6] Auditing Viral and Organelle Features...")

base_features <- rownames(library_list[[1]])
viral_dict <- config$scRNA_params$viral_feature_dictionary

if (!all(viral_dict %in% base_features)) stop("CRITICAL: Configured viral features are missing from matrix.")

library_list <- lapply(library_list, function(seu) {
  # Viral
  seu$Viral_UMIs      <- colSums(GetAssayData(seu, layer = "counts")[viral_dict, , drop = FALSE])
  seu$Viral_Burden_Xc <- log1p(seu$Viral_UMIs)
  
  # Organelles
  mt_genes <- grep("^MT-", rownames(seu), value = TRUE)
  seu[["percent.mt"]] <- if(length(mt_genes) > 0) PercentageFeatureSet(seu, features = mt_genes) else NA
  
  ribo_genes <- grep("^RP[SL]", rownames(seu), value = TRUE)
  seu[["percent.ribo"]] <- if(length(ribo_genes) > 0) PercentageFeatureSet(seu, features = ribo_genes) else NA
  
  return(seu)
})

# ==============================================================================
# STEP 3: Zero-Safeguarded Adaptive QC & Dynamic Complexity Flag
# ==============================================================================
message("--> [3/6] Applying Zero-Safeguarded Adaptive QC...")

calc_mad_floor <- function(x, abs_floor) {
  mad_val <- mad(x, na.rm = TRUE)
  if (mad_val <= 1e-5) return(abs_floor) 
  return(max(abs_floor, median(x, na.rm = TRUE) - (config$scRNA_params$qc$mad_threshold * mad_val)))
}

library_list <- lapply(library_list, function(seu) {
  f_floor <- calc_mad_floor(seu$nFeature_RNA, config$scRNA_params$qc$absolute_min_features)
  c_floor <- calc_mad_floor(seu$nCount_RNA, config$scRNA_params$qc$absolute_min_counts)
  
  qc_pass <- seu$nFeature_RNA >= f_floor & seu$nCount_RNA >= c_floor
  
  # Symmetrical Ribosomal and Mitochondrial upper-tail filtering
  if (!all(is.na(seu$percent.ribo))) {
    r_ceil  <- median(seu$percent.ribo, na.rm=T) + (config$scRNA_params$qc$mad_threshold * max(1e-5, mad(seu$percent.ribo, na.rm=T)))
    qc_pass <- qc_pass & (seu$percent.ribo <= r_ceil)
  }
  if (!all(is.na(seu$percent.mt))) {
    mt_ceil <- median(seu$percent.mt, na.rm=T) + (config$scRNA_params$qc$mad_threshold * max(1e-5, mad(seu$percent.mt, na.rm=T)))
    qc_pass <- qc_pass & (seu$percent.mt <= mt_ceil)
  }
  
  seu$QC_Drop <- !qc_pass
  
  # Dynamic High-Complexity Flag (Replaces hardcoded 8000 threshold)
  mad_val_f <- mad(seu$nFeature_RNA, na.rm = TRUE)
  adaptive_high_cutoff <- median(seu$nFeature_RNA, na.rm = TRUE) + (3 * max(1e-5, mad_val_f))
  seu$High_Complexity_Flag <- seu$nFeature_RNA > adaptive_high_cutoff
  
  return(seu)
})

# ==============================================================================
# STEP 4: Doublet Screening (tryCatch safeguarded)
# ==============================================================================
message("--> [4/6] Running Doublet Screening...")

library_list <- lapply(library_list, function(seu) {
  if (ncol(seu) < 50) {
    seu$doublet_class <- "untested"
    return(seu)
  }
  
  sce <- as.SingleCellExperiment(seu)
  sce <- tryCatch({
    scDblFinder(sce)
  }, error = function(e) {
    sce$scDblFinder.class <- "failed"
    return(sce)
  })
  
  seu$doublet_class <- sce$scDblFinder.class
  return(seu)
})

# ==============================================================================
# STEP 5: Merge, Drop Debris/Doublets, & Export Metadata
# ==============================================================================
message("--> [5/6] Exporting Audit Metadata and Cleaning Object...")

# Merge list into a single un-integrated Seurat Object
merged_unfiltered <- merge(library_list[[1]], y = library_list[2:8], project = "MeV_scRNA")

# Seurat v5 Fix: Join layers to prevent PCA silent dropping warnings later
merged_unfiltered <- JoinLayers(merged_unfiltered)

# Export Full Audit Trail for Reviewer Diagnostics
write.csv(merged_unfiltered@meta.data, file.path(dir_results, "scRNA_PreQC_Audit_Metadata.csv"), row.names = FALSE)

# Clean object: Retain singlets and safely handle untested/failed doublets
merged_clean <- subset(merged_unfiltered, QC_Drop == FALSE & doublet_class %in% c("singlet", "untested", "failed"))

saveRDS(merged_clean, file.path(dir_processed, "GSE168775_QC_Cleaned.rds"))

# ==============================================================================
# STEP 6: Provenance
# ==============================================================================
time_stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
writeLines(capture.output(sessionInfo()), file.path(dir_logs, paste0("SessionInfo_02_", time_stamp, ".txt")))

message(sprintf("SUCCESS: Retained %d clean cells across %d libraries.", ncol(merged_clean), length(unique(merged_clean$Library_ID))))