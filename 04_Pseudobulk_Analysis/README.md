# 04 — Pseudobulk Eligibility, Aggregation, and Differential Expression

## 1. Overview

This work package develops a sample-aware differential-expression workflow for the single-cell RNA-sequencing dataset used in the Measles Virus (MeV) host-response thesis.

The analysis has two main components:

1. **Pseudobulk eligibility auditing:** Determine whether cell populations have sufficient representation across libraries and whether the intended experimental contrasts can be evaluated using the available sample structure.
2. **Pseudobulk aggregation and contrast-specific analysis:** Aggregate single-cell counts at the appropriate sample or library level and conduct differential-expression analyses for the specified biological comparisons.

The central methodological motivation is to avoid treating individual cells from the same biological sample as independent biological replicates. Cell-level expression measurements provide resolution, but biological replication must be represented appropriately in statistical testing.

## 2. Dataset and Experimental Structure

- **Dataset:** GSE168775
- **Data type:** Single-cell RNA sequencing
- **Biological system:** Human airway epithelial cultures exposed to measles virus
- **Experimental conditions:** Mock, Bystander, and Infected
- **Sampling time points:** Day 3, Day 7, and Day 14
- **Primary input:** The annotated single-cell Seurat object produced by the preceding clustering and annotation workflow.

The analysis uses cell annotations, library identifiers, experimental condition, and sampling-day metadata to evaluate eligibility and define the intended comparisons.

## 3. Scripts

| Script | Purpose |
|---|---|
| `06c_Pseudobulk_Eligibility_Audit.R` | Audits cell-count eligibility, condition and time-point coverage, and model-matrix estimability at the cluster and broad-lineage levels. |
| `07_Pseudobulk_Aggregation_Constrast_Specefic.R` | Performs pseudobulk aggregation and contrast-specific downstream analysis using the eligibility decisions and sample structure established in the preceding audit. |

The original script filename is retained for traceability.

## 4. Workflow

### Step 1 — Prepare the cell-composition table

The eligibility-audit workflow summarizes cell counts using the available cluster, broad-lineage, library, condition, and sampling-day metadata.

These summaries establish how many cells are available for each population in each library and experimental group.

### Step 2 — Audit minimum cell-count thresholds

Eligibility is evaluated using three candidate minimum cell-count thresholds:

- 10 cells
- 20 cells
- 30 cells

The audit is performed for individual clusters and broad lineages. This comparison documents how the choice of a minimum cell-count requirement affects which population–library combinations remain eligible for downstream analysis.

The script also exports an explicit list of broad-lineage/library combinations with fewer than 10 cells, providing a record of the exclusions associated with the threshold selected for the downstream workflow.

A minimum cell count is an operational eligibility rule; it does not, by itself, establish adequate statistical power.

### Step 3 — Evaluate condition and time-point coverage

For each target population and threshold, the audit examines whether the required experimental conditions are present and whether the conditions share sampling time points.

The specified comparisons include:

- Bystander versus Mock.
- Infected versus Bystander.

The audit records missing conditions, missing days, shared time points, and the number of observations contributing to each comparison.

This step helps distinguish a comparison that lacks the required experimental groups from one in which both groups exist but do not share a suitable sampling time point.

### Step 4 — Assess model-matrix estimability

The audit evaluates condition-only and day-by-condition model matrices, recording their dimensions, rank, residual degrees of freedom, and full-rank status.

These diagnostics help identify cases where the available experimental structure may not support the intended model or contrast.

The presence of both conditions is not sufficient to guarantee an estimable contrast. Model-matrix diagnostics are therefore used alongside the cell-count and time-point audits.

### Step 5 — Aggregate counts for pseudobulk analysis

The subsequent aggregation workflow uses the annotated single-cell data and sample metadata to construct pseudobulk expression profiles for the eligible populations.

The key principle is that cells are aggregated within the appropriate biological sample or library and target population, rather than treated as independent biological replicates in the differential-expression test.

The aggregation unit and metadata must remain consistent with the experimental design. Multiple libraries should not automatically be assumed to represent independent biological donors unless the study design supports that interpretation.

### Step 6 — Conduct contrast-specific differential-expression analysis

The downstream workflow evaluates the specified biological contrasts using the eligible pseudobulk profiles and the corresponding experimental metadata.

The interpretation of each contrast depends on its design matrix, the samples retained after eligibility filtering, and the availability of the relevant conditions and time points.

Results from Bystander-versus-Mock and Infected-versus-Bystander comparisons address different biological questions and should be interpreted separately.

### Step 7 — Document eligibility and exclusions

The audit outputs preserve information about population-level cell counts, threshold decisions, condition coverage, and model-matrix diagnostics.

These records provide a transparent basis for documenting why a population or comparison was included, excluded, or flagged for review before downstream statistical interpretation.

## 5. Main Audit Outputs

The eligibility-audit script explicitly writes the following tables:

| Output | Description |
|---|---|
| `Threshold_Audit_Clusters.csv` | Cluster-level eligibility and model-matrix diagnostics at the 10-, 20-, and 30-cell thresholds. |
| `Threshold_Audit_Broad_Lineages.csv` | Broad-lineage eligibility and model-matrix diagnostics at the same thresholds. |
| `Library_Exclusions_N10_Lineage.csv` | Broad-lineage/library combinations with fewer than 10 cells, documented for the threshold used by the downstream workflow. |

The contrast-specific aggregation script may generate additional pseudobulk matrices, differential-expression tables, and diagnostic outputs. Their exact filenames should be documented after checking the script and confirming which files were generated.

## 6. Methodological Considerations

- **Biological replication:** The sample or biological replicate, rather than the individual cell, must define the appropriate unit of replication for the differential-expression model.
- **Eligibility versus power:** Meeting a cell-count threshold does not guarantee sufficient biological replication or adequate statistical power.
- **Model estimability:** Condition and time-point coverage must be assessed before interpreting a contrast.
- **Contrast specificity:** Bystander-versus-Mock and Infected-versus-Bystander comparisons address different questions and require separate interpretation.
- **Sample structure:** Library identifiers, donor identities, and biological replicates are not necessarily interchangeable. Their relationship must be established from the study metadata.
- **Exclusion transparency:** Population–library combinations removed by the eligibility rules should be documented rather than silently discarded.
- **Statistical interpretation:** Differential-expression results should be interpreted alongside effect sizes, uncertainty, multiple-testing correction, and the number of biological replicates supporting each contrast.

## 7. Relationship to the Overall Thesis

This work package connects single-cell clustering and annotation to sample-aware differential-expression analysis.

It establishes which populations and comparisons can be considered under the defined eligibility rules, then provides the framework for evaluating transcriptional differences across the experimental conditions.

The resulting evidence contributes to the broader thesis objective of investigating measles-associated host responses while distinguishing cell-level patterns from conclusions supported by sample-level replication.

Downstream cross-dataset analyses can then assess whether selected biological signals are consistent with the independent bulk transcriptomic evidence.

## 8. Reproducibility and Status

The workflow is implemented in R and depends on the annotated single-cell object, experimental metadata, configuration and package requirements, and the outputs of the preceding analysis stages.

Large intermediate expression objects do not need to be included in the public repository if their size or sharing restrictions make that impractical. Instead, document the source of each input and the steps required to generate it.

**Status:** Pseudobulk eligibility and contrast-specific analysis workflow. The exact retained populations, sample counts, differential-expression results, and final biological conclusions should be reported only after checking the generated audit tables and analysis outputs.
