# Reproducibility audit report

Audit date: 2026-07-22 to 2026-07-24

## Scope

The audit used the final manuscript, supplementary files, frozen inputs, and release package as the authoritative sources. It covered code portability, frozen input identity, dynamic-program membership, donor-aware trajectory analyses, WGCNA, GSE191067 state localization and module preservation, PRJNA813172 recurrence projection, and final machine-readable outputs.

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

## Release integrity measures

1. Analysis scripts use repository-relative paths and documented environment variables.
2. Clean-program analyses explicitly select `clean_module == clean_D1` or `clean_module == clean_D2` and assert 42/975 members.
3. The trajectory sensitivity script extracts the 1,586 ordering genes from the frozen gene-definition table.
4. `modulePreservation` uses a documented checkpoint/resume design without an undocumented WGCNA cache.
5. The module-preservation output field is standardized as `program_alignment` across tables and figure-building code.
6. Public scripts cover raw-data acquisition, GENCODE v50 reference construction, Salmon quantification, and tximport aggregation for PRJNA813172.

## v1.0.1 clean D2 consistency update

All discovery donor-level clean D2 outputs were regenerated from the frozen 975-gene definition. Supplementary Data S4 and the corresponding machine-readable result files now report:

| sample | corrected median_clean_D2 |
|---|---:|
| K183156 | 0.395836368629022 |
| K183261 | 0.288458706493948 |
| K712771 | 0.343506740775223 |
| N176712 | 0.248634960880203 |
| N177536 | 0.215118331737006 |
| N387690 | 0.135701782579501 |

The corrected donor-level group summary has normal median `0.215118331737006`, keloid median `0.343506740775223`, mean difference `0.142782246900494`, exact Wilcoxon `P = 0.1`, exact permutation `P = 0.1`, and Cliff's delta `1.0`. The clean D1 values, pseudotime values, state fractions, group direction, main figures, WGCNA, GSE191067, and bulk recurrence results are unchanged. Supplementary Data S4 was visually inspected after export, and the release audit now verifies agreement among S4, donor-level result tables, doublet-sensitivity outputs, and the supplementary checksum index.

## Publication status

The public repository is `https://github.com/y504071842/keloid-fibroblast-trajectory`. Versioned archives are registered under the Zenodo concept DOI `10.5281/zenodo.21503330`.
