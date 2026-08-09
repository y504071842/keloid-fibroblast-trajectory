# Data access and frozen inputs

## Public source datasets

- **GSE191067**: external single-cell RNA-seq cohort, available from NCBI GEO.
- **PRJNA813172**: postoperative recurrence bulk RNA-seq cohort, available from NCBI SRA.

The bulk FASTQ files can be reconstructed from PRJNA813172 using the scripts in `analysis/`. GSE191067 should be acquired from GEO under the accession above and processed according to the manuscript and supplementary methods.

## Discovery-cohort sequencing data

Raw single-cell RNA-sequencing data from the discovery cohort have been deposited in the Genome Sequence Archive for Human (GSA-Human) under BioProject accession **PRJCA070046** and GSA-Human accession **HRA020030** (<https://ngdc.cncb.ac.cn/gsa-human/browse/HRA020030>).

## Frozen discovery inputs

The final discovery analyses start from three processed objects:

- `discovery_fibroblast_subclustered.rds`
- `discovery_trajectory_programDE_cds.rds`
- `discovery_trajectory_graph_pseudotime_cds.rds`

Expected sizes and SHA-256 checksums are recorded in `manifests/input_objects.csv`. These patient-derived objects are intentionally excluded from Git and are subject to the applicable ethics, consent, and institutional data-access requirements.

Place authorized copies in `data/frozen_inputs/` or set `KFT_INPUT_ROOT` to their directory.

## External processed object

External validation scripts use `GSE191067_standard_reclustered_scored.rds`. This derived object is not committed because it is large and can be reconstructed from the public accession. Set `KFT_GSE191067_OBJECT` to its local path.

## Bulk expression inputs

After Salmon and tximport processing, set:

- `PRJNA813172_LOG2TPM`: gene-symbol `log2(TPM + 1)` matrix;
- `PRJNA813172_METADATA`: recurrence metadata keyed by SRA Run ID.

No expression filtering is applied before score projection. Multiple transcripts mapping to the same gene symbol are aggregated by `tximport`; symbols are not made unique after matrix construction.
