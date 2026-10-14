# ==============================================================================
# SCRIPT 05: UNSUPERVISED CLUSTER MARKER CHARACTERIZATION & RESOLUTION AUDIT
# ==============================================================================

library(here)
library(yaml)
library(tidyverse)
library(Seurat)
library(patchwork)

# ------------------------------------------------------------------------------
# 1. Environment Setup & Configuration Validation
# ------------------------------------------------------------------------------
config <- yaml::read_yaml(here::here("config.yml"))

# Enforce explicit config parameters (fixed parameters and random seed for reproducibility)
req_keys <- c("reduction", "dims", "k_param", "resolution", "sensitivity_resolutions", 
              "algorithm", "n_start", "n_iter", "random_seed")
stopifnot("CRITICAL: Missing required clustering keys in config.yml!" = all(req_keys %in% names(config$scRNA_params$clustering)))

set.seed(config$scRNA_params$clustering$random_seed)

# Environment Guard: Enforce Seurat >= 5.5.1 as the validated pipeline environment for PrepSCTFindMarkers
stopifnot("CRITICAL: Seurat version < 5.5.1 detected. Please use the validated pipeline environment >= 5.5.1." = packageVersion("Seurat") >= "5.5.1")

dir_processed <- here::here(config$directories$data_processed)
dir_results   <- here::here(config$directories$results_sc)
dir_plots     <- here::here(config$directories$plots)
dir_logs      <- here::here(config$directories$logs)

for (d in c(dir_processed, dir_results, dir_plots, dir_logs)) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
}

# ------------------------------------------------------------------------------
# 2. Load Benchmarked Object & Validate Input Geometry
# ------------------------------------------------------------------------------
message("--> [1/9] Loading Benchmarked Object & Validating Coordinates...")
seu <- readRDS(file.path(dir_processed, "GSE168775_Integrated_Benchmarked.rds"))

target_reduction <- config$scRNA_params$clustering$reduction

# Strict Geometric Assertions
stopifnot("CRITICAL: Target reduction missing from object!" = target_reduction %in% Reductions(seu))
stopifnot("CRITICAL: Reduction embedding row count does not match cell count!" = nrow(Embeddings(seu, target_reduction)) == ncol(seu))
stopifnot("CRITICAL: Target reduction contains NA values!" = !anyNA(Embeddings(seu, target_reduction)))

# ------------------------------------------------------------------------------
# 3. Explicit Graph Construction & Multi-Resolution Clustering
# ------------------------------------------------------------------------------
message(sprintf("--> [2/9] Constructing Graph (%s, k=%d) & Running Multi-Resolution Clustering...", 
                target_reduction, config$scRNA_params$clustering$k_param))

# Explicitly name NN and SNN graphs to prevent implicit default graph selection
seu <- FindNeighbors(
  object = seu,
  reduction = target_reduction,
  dims = 1:config$scRNA_params$clustering$dims,
  k.param = config$scRNA_params$clustering$k_param,
  graph.name = c("rpca_nn", "rpca_snn"),
  verbose = FALSE
)

stopifnot("CRITICAL: rpca_snn graph missing after FindNeighbors!" = "rpca_snn" %in% names(seu@graphs))

# Execute clustering across sensitivity resolutions using fixed parameters and random seed
sens_res <- config$scRNA_params$clustering$sensitivity_resolutions

for (r in sens_res) {
  seu <- FindClusters(
    object = seu,
    graph.name = "rpca_snn",
    resolution = r,
    algorithm = config$scRNA_params$clustering$algorithm,
    n.start = config$scRNA_params$clustering$n_start,
    n.iter = config$scRNA_params$clustering$n_iter,
    random.seed = config$scRNA_params$clustering$random_seed,
    verbose = FALSE
  )
}

# ASSERTION: Verify all configured sensitivity resolution metadata columns successfully generated
expected_res_cols <- paste0("rpca_snn_res.", sens_res)
stopifnot("CRITICAL: Not all sensitivity resolution columns were generated in metadata!" = all(expected_res_cols %in% colnames(seu@meta.data)))

# Lock primary resolution into active identity and metadata column
primary_res <- config$scRNA_params$clustering$resolution
primary_col <- paste0("rpca_snn_res.", primary_res)

seu$seurat_clusters <- seu[[primary_col]][, 1]
Idents(seu) <- "seurat_clusters"

# ------------------------------------------------------------------------------
# 4. Composition, Raw-Count, & Proportion Audits
# ------------------------------------------------------------------------------
message("--> [3/9] Exporting Composition, Raw-Count, & Proportion Audits...")

# A. Basic Cluster Composition Audit
audit_table <- seu@meta.data %>%
  dplyr::group_by(seurat_clusters) %>%
  dplyr::summarise(
    Total_Cells = dplyr::n(),
    Percent_Total = round((dplyr::n() / ncol(seu)) * 100, 2),
    Unique_Libraries = dplyr::n_distinct(Library_ID),
    Mock_Cells = sum(Condition == "Mock", na.rm = TRUE),
    Bystander_Cells = sum(Condition == "Bystander", na.rm = TRUE),
    Infected_Cells = sum(Condition == "Infected", na.rm = TRUE),
    .groups = "drop"
  )
write.csv(audit_table, file.path(dir_results, "Cluster_Composition_Audit.csv"), row.names = FALSE)

# B. Cluster-by-Library Proportions (Rows = Clusters, Cols = Libraries)
cluster_by_lib_prop <- as.data.frame(unclass(round(prop.table(table(seu$seurat_clusters, seu$Library_ID), margin = 1) * 100, 2)))
write.csv(cluster_by_lib_prop, file.path(dir_results, "Cluster_by_Library_Proportions.csv"))

# C. Library-by-Cluster Proportions (Rows = Libraries, Cols = Clusters)
lib_by_cluster_prop <- as.data.frame(unclass(round(prop.table(table(seu$Library_ID, seu$seurat_clusters), margin = 1) * 100, 2)))
write.csv(lib_by_cluster_prop, file.path(dir_results, "Library_by_Cluster_Proportions.csv"))

# D. Library-by-Cluster RAW Cell Counts (Crucial for Pseudobulk DE Sparsity Gates)
lib_by_cluster_counts <- as.data.frame(unclass(table(seu$Library_ID, seu$seurat_clusters)))
write.csv(lib_by_cluster_counts, file.path(dir_results, "Library_by_Cluster_Cell_Counts.csv"))

# E. Cluster x Day x Condition RAW Cell Counts
day_cond_counts <- as.data.frame(table(Cluster = seu$seurat_clusters, Day = seu$Day, Condition = seu$Condition))
write.csv(day_cond_counts, file.path(dir_results, "Cluster_Day_Condition_Cell_Counts.csv"), row.names = FALSE)

# ------------------------------------------------------------------------------
# 5. Cluster Correspondence & Overlap Diagnostics (0.6 vs 0.4 and 0.6 vs 0.8)
# ------------------------------------------------------------------------------
message("--> [4/9] Exporting Resolution Correspondence / Overlap Diagnostics...")

tab_0.6_vs_0.4 <- table(Res_0.6 = seu$rpca_snn_res.0.6, Res_0.4 = seu$rpca_snn_res.0.4)
write.csv(as.data.frame(unclass(tab_0.6_vs_0.4)), file.path(dir_results, "Cluster_Correspondence_0.6_vs_0.4.csv"))

tab_0.6_vs_0.8 <- table(Res_0.6 = seu$rpca_snn_res.0.6, Res_0.8 = seu$rpca_snn_res.0.8)
write.csv(as.data.frame(unclass(tab_0.6_vs_0.8)), file.path(dir_results, "Cluster_Correspondence_0.6_vs_0.8.csv"))

max_overlap_0.4 <- apply(tab_0.6_vs_0.4, 1, max) / rowSums(tab_0.6_vs_0.4) * 100
max_overlap_0.8 <- apply(tab_0.6_vs_0.8, 1, max) / rowSums(tab_0.6_vs_0.8) * 100

overlap_df <- data.frame(
  Cluster_0.6 = rownames(tab_0.6_vs_0.4),
  Max_Overlap_Pct_in_0.4 = round(max_overlap_0.4, 2),
  Max_Overlap_Pct_in_0.8 = round(max_overlap_0.8, 2)
)
write.csv(overlap_df, file.path(dir_results, "Resolution_Correspondence_Max_Overlap.csv"), row.names = FALSE)

# ------------------------------------------------------------------------------
# 6. Save Pure Clustered Object (Prior to SCT DE Recorrection)
# ------------------------------------------------------------------------------
message("--> [5/9] Saving Pure Clustered Geometry Object (GSE168775_Clustered_Provisional.rds)...")
saveRDS(seu, file.path(dir_processed, "GSE168775_Clustered_Provisional.rds"))

# ------------------------------------------------------------------------------
# 7. SCT Model Recorrection & Granularity Audits
# ------------------------------------------------------------------------------
message("--> [6/9] Preparing SCT Models for Exploratory Diagnostics & DE...")

# Create a dedicated copy to protect the pure geometry object stored on disk
seu_marker <- seu

DefaultAssay(seu_marker) <- "SCT"
stopifnot("CRITICAL: Default assay is not SCT!" = DefaultAssay(seu_marker) == "SCT")

Idents(seu_marker) <- "seurat_clusters"

stopifnot("CRITICAL: seurat_clusters missing from metadata!" = "seurat_clusters" %in% colnames(seu_marker@meta.data))
stopifnot("CRITICAL: Cluster identities contain NA!" = !anyNA(Idents(seu_marker)))
stopifnot("CRITICAL: Cell-wise identity mismatch!" = all(as.character(Idents(seu_marker)) == as.character(seu_marker$seurat_clusters)))

# Strictly technical recorrection step: prepares the multiple independent SCT models for comparable downstream DE.
seu_marker <- PrepSCTFindMarkers(seu_marker, assay = "SCT", verbose = FALSE)

message("--> [7/9] Auditing Marker Program Granularity vs. Disappearance...")

expr_mat <- GetAssayData(seu_marker, assay = "SCT", layer = "data")
stopifnot("CRITICAL: Expression matrix dimension mismatch!" = ncol(expr_mat) == ncol(seu_marker))

# NOTE: >0.5 log-SCT is explicitly an exploratory high-expression threshold.
# It tracks whether a marker program persists structurally when clusters merge.
target_programs <- list(
  Deuterosomal = c("CDC20B", "DEUP1", "FOXN4"),
  Mucous       = c("MUC5B"),
  Cycling      = c("MKI67", "TOP2A"),
  Ciliated     = c("CAPS", "FOXJ1"),
  Basal        = c("KRT5", "TP63", "KRT15"),
  Secretory    = c("BPIFA1", "SCGB1A1")
)

granularity_results <- list()
for (prog_name in names(target_programs)) {
  genes <- intersect(target_programs[[prog_name]], rownames(expr_mat))
  if (length(genes) > 0) {
    
    if (length(genes) > 1) {
      high_exp_mask <- colSums(expr_mat[genes, , drop = FALSE] > 0.5) >= 2
    } else {
      high_exp_mask <- colSums(expr_mat[genes, , drop = FALSE] > 0.5) >= 1
    }
    
    for (r in sens_res) {
      col_name <- paste0("rpca_snn_res.", r)
      df <- seu_marker@meta.data %>%
        mutate(High_Expressing = high_exp_mask) %>%
        dplyr::group_by(Cluster = .data[[col_name]]) %>%
        dplyr::summarise(
          Program = prog_name,
          Resolution = r,
          Total_Cells = dplyr::n(),
          High_Expressing_Cells = sum(High_Expressing),
          Pct_High_Expressing = round((High_Expressing_Cells / Total_Cells) * 100, 2),
          .groups = "drop"
        )
      granularity_results[[paste(prog_name, r, sep = "_")]] <- df
    }
  }
}
write.csv(bind_rows(granularity_results), file.path(dir_results, "Granularity_vs_Disappearance_Audit.csv"), row.names = FALSE)

dist_genes <- intersect(c("CDC20B", "DEUP1", "FOXN4", "MUC5B"), rownames(seu_marker))
if(length(dist_genes) > 0) {
  p_dist <- VlnPlot(seu_marker, features = dist_genes, group.by = "rpca_snn_res.0.4", assay = "SCT", pt.size = 0.1, stack = TRUE, flip = TRUE) +
    ggtitle("Granularity vs Disappearance: Distribution Audit at Res 0.4")
  ggsave(file.path(dir_plots, "Granularity_vs_Disappearance_Distributions_0.4.png"), p_dist, width = 10, height = 8, dpi = 300)
}

ambig_clusters <- c("1", "3", "8", "10", "11", "14", "15")
avail_ambig <- intersect(ambig_clusters, levels(seu_marker$rpca_snn_res.0.6))

if (length(avail_ambig) > 0) {
  ambig_tracker <- seu_marker@meta.data %>%
    filter(rpca_snn_res.0.6 %in% avail_ambig) %>%
    dplyr::group_by(Res_0.6_Cluster = rpca_snn_res.0.6, Res_0.4_Cluster = rpca_snn_res.0.4, Res_0.8_Cluster = rpca_snn_res.0.8) %>%
    dplyr::summarise(Cell_Count = dplyr::n(), .groups = "drop")
  write.csv(ambig_tracker, file.path(dir_results, "Ambiguous_Cluster_Resolution_Tracker.csv"), row.names = FALSE)
}

key_markers <- unlist(config$scRNA_params$annotation$marker_panels) %>% unname() %>% unique()
avail_key_markers <- intersect(key_markers, rownames(seu_marker))

p_res_list <- list()
for (r in sens_res) {
  col_name <- paste0("rpca_snn_res.", r)
  p_res_list[[as.character(r)]] <- DotPlot(seu_marker, features = avail_key_markers, group.by = col_name) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8)) +
    labs(title = sprintf("Resolution %s", r), x = NULL, y = "Cluster")
}
p_res_grid <- wrap_plots(p_res_list, ncol = 1)
ggsave(file.path(dir_plots, "Resolution_Sensitivity_DotPlot.png"), p_res_grid, width = 14, height = 16, dpi = 300)

# ------------------------------------------------------------------------------
# 8. Unsupervised Cluster Marker Characterization
# ------------------------------------------------------------------------------
message("--> [8/9] Extracting Positive Cluster Markers (Non-integrated expression)...")

cluster_markers <- FindAllMarkers(
  object = seu_marker,
  assay = "SCT",
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25,
  verbose = FALSE
)

stopifnot("CRITICAL: FindAllMarkers returned an empty table!" = nrow(cluster_markers) > 0)
stopifnot("CRITICAL: Reported clusters in marker table do not match active clusters!" = all(unique(cluster_markers$cluster) %in% levels(seu_marker$seurat_clusters)))
stopifnot("CRITICAL: Gene symbols contain NA or blank strings!" = all(!is.na(cluster_markers$gene) & cluster_markers$gene != ""))
stopifnot("CRITICAL: Missing expected statistic columns in marker table!" = all(c("p_val_adj", "avg_log2FC", "gene", "cluster") %in% colnames(cluster_markers)))

write.csv(cluster_markers, file.path(dir_results, "Positive_Cluster_Markers.csv"), row.names = FALSE)

# ------------------------------------------------------------------------------
# 9. Summary Diagnostics & Provenance Export
# ------------------------------------------------------------------------------
message("--> [9/9] Generating Final Diagnostics and Provenance Logs...")

clusters_detected <- length(levels(seu_marker$seurat_clusters))
clusters_with_markers <- length(unique(cluster_markers$cluster))
total_markers <- nrow(cluster_markers)

top_markers <- cluster_markers %>%
  dplyr::group_by(cluster) %>%
  dplyr::slice_max(n = 1, order_by = avg_log2FC, with_ties = FALSE) %>%
  dplyr::select(cluster, gene, avg_log2FC, p_val_adj)

message("\n==========================================================================")
message("   UNSUPERVISED CLUSTER MARKER CHARACTERIZATION REPORT")
message("==========================================================================")
message(sprintf("   Clusters Detected                  : %d", clusters_detected))
message(sprintf("   Clusters with >=1 Positive Marker  : %d", clusters_with_markers))
message(sprintf("   Total Positive Markers Identified  : %d", total_markers))
message("   Top Marker Per Cluster:")
print(as.data.frame(top_markers))
message("==========================================================================\n")

probe_genes <- c("CDC20B", "DEUP1", "FOXN4", "MKI67", "TOP2A", "PBK", "ASPM", "FOXI1", "POU2F3")
avail_probes <- intersect(probe_genes, rownames(seu_marker))
if (length(avail_probes) > 0) {
  p_probe <- FeaturePlot(seu_marker, features = avail_probes, reduction = "umap_rpca", pt.size = 0.4, order = TRUE, ncol = 3)
  ggsave(file.path(dir_plots, "Rare_Population_Probe.png"), p_probe, width = 12, height = 10, dpi = 300)
}

p_umap <- DimPlot(seu_marker, reduction = "umap_rpca", label = TRUE, repel = TRUE) +
  ggtitle(sprintf("Unsupervised Clusters (RPCA Space, Res: %s)", primary_res))
ggsave(file.path(dir_plots, "Cluster_UMAP.png"), p_umap, width = 8, height = 6, dpi = 300)

time_stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
prov_log <- c(
  paste("Date:", time_stamp),
  paste("Seurat Version:", packageVersion("Seurat")),
  paste("Reduction:", target_reduction),
  paste("Graph Name:", "rpca_snn"),
  paste("Dimensions:", config$scRNA_params$clustering$dims),
  paste("k_param:", config$scRNA_params$clustering$k_param),
  paste("Primary Resolution:", primary_res),
  paste("Sensitivity Resolutions:", paste(sens_res, collapse = ", ")),
  paste("Algorithm:", config$scRNA_params$clustering$algorithm),
  paste("n_start:", config$scRNA_params$clustering$n_start),
  paste("n_iter:", config$scRNA_params$clustering$n_iter),
  paste("Random Seed:", config$scRNA_params$clustering$random_seed),
  paste("Total Input Cells:", ncol(seu)),
  paste("Clusters Identified:", clusters_detected)
)

writeLines(c(prov_log, "", capture.output(sessionInfo())), con = file.path(dir_logs, paste0("SessionInfo_05_", time_stamp, ".txt")))

message("SUCCESS: Script 05 refactoring complete. All requested audits rendered.")
