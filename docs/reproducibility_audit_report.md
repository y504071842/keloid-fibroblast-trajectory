# Reproducibility audit report

Audit date: 2026-07-22 to 2026-07-23

## Scope

The audit used `manuscript 2.0.docx` and the final package under `initial draft 7_11v2` as the authoritative manuscript and result sources. It covered code portability, frozen input identity, dynamic-program membership, donor-aware trajectory analyses, WGCNA, GSE191067 state localization and module preservation, PRJNA813172 recurrence projection, and final machine-readable outputs.

## Passed checks

- All three discovery RDS inputs, the 1.33-GB GSE191067 processed object, the PRJNA813172 `log2(TPM + 1)` matrix and metadata, and the NCBI run report matched their recorded SHA-256 checksums.
- The NCBI run report independently yielded exactly 77 SRR accessions.
- Frozen membership was confirmed as 42 clean D1 genes and 975 clean D2 genes from 1,118 dynamic candidates.
- All 20 R scripts parsed successfully under R 4.5.2.
- The donor-aware trajectory pipeline reproduced every final CSV with zero numeric difference.
- The two-lineage tradeSeq pipeline reproduced all association, pattern, end-point, early-differentiation, estimability, lineage-classification, and cell-metadata tables with zero numeric difference.
- Ordering-gene sensitivity and Slingshot reproduced the final summaries and cell-level tables with zero numeric difference after correcting the public script's input-column lookup.
- Discovery WGCNA reproduced module assignments, aligned modules, donor-equal correlations, bootstrap summaries, gene overlap, and network parameters with zero numeric difference.
- GSE191067 reclustering, state mapping, program coverage, donor-state scores, and state composition reproduced the final values. Files with different byte-level hashes were semantically identical after type-aware comparison.
- The 2,000-permutation module-preservation analysis reproduced brown `Zsummary = 2.46743` and turquoise `Zsummary = 12.85331`; all six leave-one-donor-out tables reproduced with zero numeric difference.
- The PRJNA813172 recurrence analysis reproduced all five final CSV files byte-for-byte, including `AUC = 0.7306122`, `AP = 0.3048002`, permutation P values, bootstrap intervals, and leave-one-recurrence-out results.
- Independent Python calculations reproduced the principal recurrence and external-composition point estimates.

## Code defects corrected in the release package

1. Legacy absolute/package-relative paths were replaced with repository-relative paths and documented environment variables.
2. `trajectory_root_and_donor_analysis.R` and `fibroblast_doublet_sensitivity.R` previously selected genes by the broad `dynamic_module` label. The release scripts now explicitly select `clean_module == clean_D1` or `clean_module == clean_D2` and assert 42/975 members.
3. `trajectory_ordering_gene_sensitivity.R` previously requested a nonexistent `ordering_gene` column from the final table. It now extracts `gene` from the `program_DE_all1586` rows and asserts 1,586 ordering genes.
4. `modulePreservation` previously attempted to write an undocumented WGCNA cache in the current directory. The release script disables this nonessential cache and supports checkpoint/resume without changing the 2,000/200 permutation design.
5. The module-preservation output field is standardized as `program_alignment`, matching the final tables and figure-building script.
6. Bulk raw-data acquisition, GENCODE v50 reference construction, Salmon quantification, and tximport aggregation were added as explicit public scripts.

## Result-file correction completed

`Supplementary_Data_S4_Trajectory_Dynamics.xlsx`, sheet `Donor_Trajectory_Summary`, originally labelled a donor median column as `median_clean_D2` although its six values had been calculated from all 1,073 D2 candidates. The released workbook now reports the corrected 975-gene clean D2 medians:

| sample | corrected median_clean_D2 |
|---|---:|
| K183156 | 0.395836368629022 |
| K183261 | 0.288458706493948 |
| K712771 | 0.343506740775223 |
| N176712 | 0.248634960880203 |
| N177536 | 0.215118331737006 |
| N387690 | 0.135701782579501 |

The clean D1 values, pseudotime values, state fractions, donor-level group direction, main figures, WGCNA, GSE191067, and bulk recurrence results are unaffected. The corrected workbook was visually inspected and its six donor values were verified after export.

## Remaining publication steps

- populate the created public GitHub repository from the final release archive;
- upload the four checksum-verified files and publish the Zenodo v1.0.0 record;
- insert the public URL and registered DOI into the manuscript availability statement.

The public repository is `https://github.com/y504071842/keloid-fibroblast-trajectory`, and `10.5281/zenodo.21503331` has been reserved for the v1.0.0 archive. The author list, affiliation, corresponding-author contact, MIT license, discovery-data access wording, and permission to archive the public-data-derived GSE191067 and PRJNA813172 files were confirmed on 2026-07-23.
