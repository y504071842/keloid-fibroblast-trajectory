# Keloid fibroblast trajectory reproducibility package

This repository accompanies the manuscript **"A single-cell trajectory framework of fibroblast state remodeling in keloids and its association with postoperative recurrence."** It contains the analysis scripts, frozen gene definitions, final numeric outputs, supplementary data workbooks, input manifests, and independent audit code used for the reported results.

## Scope

The repository supports four linked analyses:

1. discovery single-cell fibroblast purification, trajectory reconstruction, dynamic-program definition, and donor-aware robustness analyses;
2. discovery-cohort WGCNA and co-expression alignment;
3. external single-cell state localization and module-preservation analysis in GSE191067;
4. recurrence-associated bulk RNA-seq projection in PRJNA813172.

Clean D1 and clean D2 are frozen at 42 and 975 genes, respectively. The public gene definitions are in `data/derived/clean_D1D2_frozen_genes.csv`. Final reference outputs are in `results/`; Supplementary Methods, the Supplementary Data S1-S8 workbooks, and their checksum index are in `reference/supplementary_data/`.

## Reproduction levels

### Numeric audit without restricted inputs

Install Python 3.11 or later with `numpy`, `pandas`, and `openpyxl`, then run:

```powershell
python tests/independent_numeric_audit.py
python tests/audit_release.py
```

These checks independently recompute the headline donor-level and recurrence statistics from the frozen tabular outputs and verify the release manifest.

### Analysis from frozen processed inputs

The discovery scripts require three frozen RDS objects listed in `data/manifests/input_objects.csv`. Place them in `data/frozen_inputs/`, or set `KFT_INPUT_ROOT` to their directory. External single-cell scripts require the processed GSE191067 Seurat object; place it in `data/frozen_inputs/` or set `KFT_GSE191067_OBJECT`.

```powershell
$env:KFT_INPUT_ROOT = "D:/path/to/discovery_inputs"
$env:KFT_GSE191067_OBJECT = "D:/path/to/GSE191067_standard_reclustered_scored.rds"
Rscript analysis/trajectory_root_and_donor_analysis.R
```

The bulk pipeline can be reconstructed from PRJNA813172 with the scripts prefixed `bulk_01` through `bulk_04`, followed by `reanalyse_bulk_primary_translation.R`. Set `PRJNA813172_LOG2TPM` and `PRJNA813172_METADATA` before the final projection.

## Repository map

- `analysis/`: executable R, Python, and PowerShell analysis scripts.
- `config/`: environment-variable contract and random seeds.
- `data/`: frozen gene definitions, input manifests, and data-access notes.
- `environment/`: R session information and software-version summary.
- `reference/supplementary_data/`: Supplementary Methods, final Supplementary Data S1-S8 workbooks, and their checksum index.
- `results/`: final machine-readable outputs used by the manuscript and figures.
- `tests/`: independent numeric and release-integrity audits.
- `docs/`: script-to-method mapping and draft availability statements.

## Data access

GSE191067 and PRJNA813172 are public. Discovery-cohort sequencing data and processed objects are available from the corresponding author upon reasonable request, subject to institutional ethics approval and applicable participant-consent restrictions; see `data/README.md`. No patient-derived RDS object is committed to this repository.

## Software

The final analysis used R 4.5.2 on Windows 11. Exact package versions are recorded in `environment/R_sessionInfo.txt`; principal software and command-line tool versions are summarized in `environment/software_versions.tsv`.

## License

Analysis code is released under the MIT License. Public source datasets remain subject to their original repository terms. Manuscript text, figures, and supplementary data are not relicensed by the code license.

## Citation

The accompanying authors are Yutong Yuan, Yuanbo Liu, Tianxiao Wei, Shan Zhu, and Nuo Si. Machine-readable citation metadata are provided in `CITATION.cff`. The public repository is available at <https://github.com/y504071842/keloid-fibroblast-trajectory>, and the latest archived release is available through the Zenodo concept DOI <https://doi.org/10.5281/zenodo.21503330>.
