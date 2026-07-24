# Methods-to-scripts map

| Manuscript analysis | Principal script(s) | Principal outputs |
|---|---|---|
| Fibroblast QC, Harmony, and doublet sensitivity | `fibroblast_QC_Harmony_sensitivity.R`; `fibroblast_doublet_sensitivity.R` | `results/QC_Harmony/` |
| Root selection and donor-level pseudotime summary | `trajectory_root_and_donor_analysis.R` | `results/discovery_donor_level_summary.csv`; root-sensitivity tables |
| Ordering-gene and cross-method sensitivity | `trajectory_ordering_gene_sensitivity.R`; `trajectory_slingshot_sensitivity.R` | `results/ordering_gene_sensitivity/`; Slingshot tables in `results/branch_trajectory/` |
| Dynamic genes and lineage-aware tradeSeq | `trajectory_branch_tradeSeq.R`; `trajectory_tradeSeq_knot_sensitivity.R` | `results/branch_trajectory/` |
| Clean D1/D2 filtering and enrichment | `define_clean_programs_and_enrichment.R` | `results/clean_programs/` |
| Donor-aware robustness | `trajectory_donor_aware_robustness.R` | `results/donor_aware_trajectory/` |
| Discovery WGCNA | `reanalyse_discovery_WGCNA.R`; `reanalyse_discovery_WGCNA_pca_metacell_sensitivity.R` | `results/WGCNA/`; `results/WGCNA_PCA_metacell_sensitivity/` |
| GSE191067 state localization | `reanalyse_GSE191067_standard_validation.R` | `results/GSE191067_standard_validation/` |
| Cross-cohort module preservation | `reanalyse_GSE191067_standard_module_preservation.R` | module-preservation tables in `results/GSE191067_standard_validation/` |
| PRJNA813172 acquisition and quantification | `bulk_01_prepare_run_manifest.py` through `bulk_04_tximport_gene_symbol.R` | local FASTQ, Salmon `quant.sf`, and gene-symbol matrices |
| Recurrence projection and statistics | `reanalyse_bulk_primary_translation.R` | `results/bulk_recurrence/` |
| Submission panels | `build_submission_figures.R` | generated figure panels |
| Independent verification | `tests/independent_numeric_audit.py`; `tests/audit_release.py` | audit JSON/report |

Scripts use repository-relative paths. Large or controlled inputs are supplied through the environment variables documented in `data/README.md`.

