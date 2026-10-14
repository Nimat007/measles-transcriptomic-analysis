# ==============================================================================
# SCRIPT 06b: STRICT EVIDENCE-BASED ANNOTATION CORRECTION
# ==============================================================================

library(here)
library(Seurat)
library(dplyr)

message("=== ANNOTATION CORRECTION START ===")

# 1. Load the provisional object
file_in <- here::here("data", "processed", "GSE168775_Annotated_Provisional.rds")
if (!file.exists(file_in)) stop("CRITICAL ERROR: Provisional Seurat object not found.")

message("--> Loading provisional Seurat object...")
seu <- readRDS(file_in)

# 2. Apply the Evidence-Based Dictionary
message("--> Reassigning Broad_Lineage based on gene-level marker audit...")
seu@meta.data <- seu@meta.data %>%
  mutate(Broad_Lineage = case_when(
    # Basal: Includes Classic (0) and Matrix-Remodeling Basal (5)
    seurat_clusters %in% c("0", "5") ~ "Basal",
    
    # Secretory: Includes Classic and Inflammatory Secretory
    seurat_clusters %in% c("2", "4", "6", "7", "8", "9", "11") ~ "Secretory",
    
    # Ciliated: Strictly FOXJ1/PIFO/CAPS positive clusters
    seurat_clusters %in% c("1", "10", "14") ~ "Ciliated",
    
    # Unresolved: LncRNA stress, Mitochondrial stress, and Mitosis (Cycling)
    seurat_clusters %in% c("3", "12", "13", "15") ~ "Unresolved",
    
    # Fallback (Should theoretically be empty if all 0-15 are covered)
    TRUE ~ "Unresolved"
  ))

# 3. Print the new distribution to console
message("\n[ NEW LINEAGE COMPOSITION ]")
print(table(seu$seurat_clusters, seu$Broad_Lineage, useNA = "ifany"))

# 4. Save the Final Object
file_out <- here::here("data", "processed", "GSE168775_Annotated_Final.rds")
message("\n--> Saving corrected object to: ", file_out)
saveRDS(seu, file_out)

message("=== ANNOTATION CORRECTION COMPLETE ===")