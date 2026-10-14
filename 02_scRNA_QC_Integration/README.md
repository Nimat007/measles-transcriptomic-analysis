# 02 — Single-Cell RNA-seq Quality Control and Integration

## 1. Overview

This work package processes the single-cell RNA-sequencing (scRNA-seq) dataset used in the Measles Virus (MeV) host-response research project. It covers library-level data ingestion, quality-control assessment, viral transcript quantification, doublet screening, normalization, dimensionality reduction, and integration benchmarking.

The objective is to establish a quality-controlled single-cell dataset and compare alternative representations of the data before downstream clustering, cell-type annotation, and biological analysis.

## 2. Dataset

- **Dataset:** GSE168775
- **Data type:** Single-cell RNA sequencing (scRNA-seq)
- **Biological system:** Human airway epithelial cultures exposed to measles virus
- **Experimental groups:** Mock, Bystander, and Infected
- **Sampling time points:** Day 3, Day 7, and Day 14
- **Input format:** 10x Genomics HDF5 (`.h5`) count matrices
- **Library structure:** Eight sequencing libraries, each assigned an explicit library ID, time point, and experimental condition.

The dataset is used to investigate cellular transcriptional responses and distinguish patterns across mock, bystander, and infected populations.

## 3. Scripts

| Script | Purpose |
|---|---|
| `02_scRNA_ingest_QC.R` | Library ingestion, metadata validation, viral burden calculation, adaptive QC, doublet screening, and export of the cleaned object. |
| `03_scRNA_integration.R` | SCTransform normalization, PCA, RPCA integration, Harmony integration, and generation of UMAP representations. |

## 4. Workflow

### Step 1 — Configuration and metadata validation

The pipeline reads project settings from `config.yml` and uses a fixed mapping between the eight GEO library identifiers, experimental conditions, and sampling days.

The ingestion script checks that exactly eight `.h5` files are present and rejects unexpected library identifiers. Each library is initially processed as a separate Seurat object.

### Step 2 — Viral transcript quantification

A predefined viral-feature dictionary is loaded from the project configuration. The script checks that all configured viral features are present in the count matrix.

For each cell, it calculates:

- `Viral_UMIs`: total UMI counts across the configured viral features.
- `Viral_Burden_Xc`: the natural logarithm of one plus the viral UMI count, calculated as `log1p(Viral_UMIs)`.

These variables provide a quantitative description of detected viral transcript burden. Their interpretation depends on the validity of the configured feature dictionary and the limitations of the deposited count matrices.

### Step 3 — Cell-level quality-control metrics

Quality-control metrics are calculated separately for each library, including:

- `nFeature_RNA`: number of detected features per cell.
- `nCount_RNA`: total RNA UMI counts per cell.
- `percent.ribo`: percentage of counts associated with annotated ribosomal genes, where available.
- `percent.mt`: percentage of counts associated with annotated mitochondrial genes, where available.

The script applies adaptive, library-specific thresholds based on medians and median absolute deviations (MAD), together with configured absolute minimums for detected features and total counts.

Upper-tail filtering is applied to ribosomal and mitochondrial percentages when those metrics are available. A separate adaptive high-complexity flag is recorded for cells with unusually high feature counts.

**Important:** Mitochondrial QC is conditional on mitochondrial genes being present in the feature annotation. A missing mitochondrial annotation is not equivalent to a measured mitochondrial percentage of zero.

### Step 4 — Doublet screening

Doublet screening is performed independently for each library using `scDblFinder`.

Libraries with fewer than 50 cells are marked as untested. Errors during doublet screening are recorded as failed rather than immediately terminating the pipeline.

The script exports a pre-QC metadata audit and constructs a cleaned Seurat object using the configured QC flags and doublet-status selection rule. In the current implementation, cells classified as `singlet`, `untested`, or `failed` can be retained if they pass the QC filters. Consequently, failed or untested doublet screening should not be interpreted as confirmed singlet status.

### Step 5 — SCTransform normalization and unintegrated baseline

The cleaned object is loaded into the integration workflow. RNA assay layers are split by sequencing library, and SCTransform normalization is performed using the v2 flavor.

Principal component analysis (PCA) is then performed to create an unintegrated representation. The script checks that the PCA embedding has the expected number of cell rows and contains no missing values.

A UMAP representation is generated from the
