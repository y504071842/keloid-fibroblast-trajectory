source(file.path(if (length(grep("^--file=", commandArgs(FALSE), value = TRUE))) dirname(chartr("\\", "/", sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))) else getwd(), "package_config.R"))
suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(WGCNA)
  library(dplyr)
})

allowWGCNAThreads(nThreads = 4)
set.seed(20260717)

wgcna_name <- Sys.getenv("WGCNA_DIR", unset = "WGCNA")
out_name <- Sys.getenv("PRES_OUT_DIR", unset = "GSE191067_standard_validation")
wgcna_dir <- file.path(audit_root, wgcna_name)
out_dir <- file.path(audit_root, out_name)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

val <- readRDS(external_object(
  "KFT_GSE191067_OBJECT",
  file.path(repo_root, "data", "frozen_inputs", "GSE191067_standard_reclustered_scored.rds"),
  "Processed GSE191067 fibroblast Seurat object"
))
DefaultAssay(val) <- "RNA"
val <- NormalizeData(val, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
val_data <- GetAssayData(val, assay = "RNA", layer = "data")

gene_modules <- read.csv(file.path(wgcna_dir, "WGCNA_gene_module_assignment.csv"), stringsAsFactors = FALSE)
colors <- setNames(gene_modules$module, gene_modules$gene)
selected <- read.csv(file.path(wgcna_dir, "WGCNA_selected_aligned_modules.csv"), stringsAsFactors = FALSE)

disc_raw <- read.csv(file.path(wgcna_dir, "matched_discovery_donor_state_pseudobulk.csv"), check.names = FALSE, stringsAsFactors = FALSE)
rownames(disc_raw) <- disc_raw$unit
disc_pb_all <- as.matrix(disc_raw[, setdiff(colnames(disc_raw), c("unit", "cohort")), drop = FALSE])
storage.mode(disc_pb_all) <- "double"

# Modules and their membership are frozen before validation is queried. Only
# here are reference genes intersected with genes measurable in GSE191067.
common_genes <- intersect(names(colors), rownames(val_data))
disc_pb <- disc_pb_all[, common_genes, drop = FALSE]
pres_colors <- colors[common_genes]

meta <- val@meta.data
meta$cell <- rownames(meta)
meta <- meta[!is.na(meta$sample) & !is.na(meta$standard_program), , drop = FALSE]
meta$state <- as.character(meta$standard_program)
val_group <- paste(meta$sample, meta$state, sep = "__")
val_pb <- t(sapply(split(meta$cell, val_group), function(cells) {
  Matrix::rowMeans(val_data[common_genes, cells, drop = FALSE])
}))

disc_eligibility <- read.csv(file.path(wgcna_dir, "WGCNA_metacell_donor_state_eligibility.csv"), stringsAsFactors = FALSE)
disc_unit_manifest <- transform(
  disc_eligibility,
  cohort = "Discovery",
  unit = paste(sample, state, sep = "__"),
  donor = sample,
  cell_n = input_cell_n
)[, c("cohort", "unit", "donor", "state", "cell_n")]
val_unit_manifest <- meta %>%
  count(sample, state, name = "cell_n") %>%
  transmute(
    cohort = "GSE191067",
    unit = paste(sample, state, sep = "__"),
    donor = sample,
    state = state,
    cell_n = cell_n
  )
unit_manifest <- bind_rows(disc_unit_manifest, val_unit_manifest) %>%
  arrange(cohort, donor, state)

module_retention <- merge(
  as.data.frame(table(module = colors), stringsAsFactors = FALSE),
  as.data.frame(table(module = pres_colors), stringsAsFactors = FALSE),
  by = "module", all.x = TRUE, suffixes = c("_reference", "_validation")
)
names(module_retention)[2:3] <- c("reference_gene_n", "validation_detected_gene_n")
module_retention$validation_detected_gene_n[is.na(module_retention$validation_detected_gene_n)] <- 0
module_retention$retention_fraction <- module_retention$validation_detected_gene_n / module_retention$reference_gene_n
module_retention$program_alignment <- ifelse(
  module_retention$module == selected$aligned_module[selected$program == "clean_D1"], "clean_D1-aligned",
  ifelse(module_retention$module == selected$aligned_module[selected$program == "clean_D2"], "clean_D2-aligned", NA)
)

run_preservation <- function(test_expr, n_perm, seed) {
  ref_var <- apply(disc_pb, 2, var, na.rm = TRUE)
  test_var <- apply(test_expr, 2, var, na.rm = TRUE)
  keep_genes <- intersect(
    names(ref_var)[is.finite(ref_var) & ref_var > 0],
    names(test_var)[is.finite(test_var) & test_var > 0]
  )
  result <- modulePreservation(
    multiData = list(
      Discovery = list(data = as.data.frame(disc_pb[, keep_genes, drop = FALSE])),
      Validation = list(data = as.data.frame(test_expr[, keep_genes, drop = FALSE]))
    ),
    multiColor = list(Discovery = pres_colors[keep_genes]),
    referenceNetworks = 1,
    nPermutations = n_perm,
    randomSeed = seed,
    savePermutedStatistics = FALSE,
    quickCor = 0,
    verbose = 1
  )
  list(result = result, genes = keep_genes)
}

resume <- tolower(Sys.getenv("PRES_RESUME", unset = "false")) %in% c("1", "true", "yes")
main_path <- file.path(out_dir, "standard_module_preservation_2000_permutations.csv")
if (resume && file.exists(main_path)) {
  preservation <- read.csv(main_path, stringsAsFactors = FALSE, check.names = FALSE)
  main_analysis_gene_n <- unique(preservation$main_analysis_gene_n)
  stopifnot(length(main_analysis_gene_n) == 1L)
} else {
  mp_run <- run_preservation(val_pb, 2000, 20260717)
  mp <- mp_run$result
  z <- as.data.frame(mp$preservation$Z[[1]][[2]], check.names = FALSE)
  z$module <- rownames(z)
  obs <- as.data.frame(mp$preservation$observed[[1]][[2]], check.names = FALSE)
  obs$module <- rownames(obs)
  preservation <- merge(z, obs, by = "module", suffixes = c("_Z", "_observed"))
  preservation <- merge(preservation, module_retention, by = "module", all.x = TRUE)
  main_analysis_gene_n <- length(mp_run$genes)
  preservation$main_analysis_gene_n <- main_analysis_gene_n
}

# Write the main 2,000-permutation result before entering the LODO loop so the
# computationally intensive primary result is checkpointed independently.
write.csv(preservation, main_path, row.names = FALSE)

val_donor <- sub("__.*$", "", rownames(val_pb))
donor_order <- sort(unique(val_donor))
loo_rows <- vector("list", length(donor_order))
loo_path <- file.path(out_dir, "standard_module_preservation_leave_one_donor_out.csv")
loo_existing <- if (resume && file.exists(loo_path)) {
  read.csv(loo_path, stringsAsFactors = FALSE, check.names = FALSE)
} else {
  data.frame()
}
for (i in seq_along(donor_order)) {
  drop <- donor_order[[i]]
  if (nrow(loo_existing) && drop %in% loo_existing$omitted_donor) {
    zz <- loo_existing[loo_existing$omitted_donor == drop, , drop = FALSE]
  } else {
    keep <- val_donor != drop
    loo_run <- run_preservation(val_pb[keep, , drop = FALSE], 200, 20260717 + i)
    result <- loo_run$result
    zz <- as.data.frame(result$preservation$Z[[1]][[2]], check.names = FALSE)
    zz$module <- rownames(zz)
    zz$omitted_donor <- drop
    zz$analysis_gene_n <- length(loo_run$genes)
  }
  loo_rows[[i]] <- zz
  write.csv(do.call(rbind, loo_rows[seq_len(i)]), loo_path, row.names = FALSE)
}
loo <- do.call(rbind, loo_rows)

write.csv(data.frame(unit = rownames(val_pb), cohort = "GSE191067", val_pb, check.names = FALSE), file.path(out_dir, "standard_matched_external_donor_state_pseudobulk.csv"), row.names = FALSE)
write.csv(unit_manifest, file.path(out_dir, "standard_module_preservation_unit_manifest.csv"), row.names = FALSE)
write.csv(module_retention, file.path(out_dir, "standard_module_gene_retention.csv"), row.names = FALSE)
write.csv(loo, loo_path, row.names = FALSE)
write.csv(data.frame(
  parameter = c(
    "seed", "reference_module_definition", "validation_used_before_module_freezing",
    "validation_gene_intersection_stage", "validation_unit", "discovery_units", "validation_units",
    "reference_network_genes", "validation_detected_network_genes", "main_permutations", "leave_one_donor_out_permutations",
    "unit_manifest"
  ),
  value = c(
    20260717, "discovery-only frozen WGCNA modules", FALSE,
    "after module detection and module membership freezing", "sample x standard marker-annotated state",
    nrow(disc_pb), nrow(val_pb), length(colors), main_analysis_gene_n, 2000, 200,
    "standard_module_preservation_unit_manifest.csv lists cohort, donor, state, and cell count for every actual unit"
  )
), file.path(out_dir, "standard_module_preservation_parameters.csv"), row.names = FALSE)

selected_modules <- selected$aligned_module
report <- capture.output({
  cat("Strict discovery-to-validation module preservation\n")
  cat("Validation units:", nrow(val_pb), "\n")
  cat("Actual donor-state units listed in standard_module_preservation_unit_manifest.csv:", nrow(unit_manifest), "\n")
  cat("Reference network genes:", length(colors), "validation-detected:", length(common_genes), "\n")
  print(module_retention[module_retention$module %in% selected_modules, ])
  print(preservation[preservation$module %in% selected_modules, ])
  cat("\nLeave-one-validation-donor-out aligned modules:\n")
  print(loo[loo$module %in% selected_modules, c("module", "omitted_donor", "Zsummary.pres")])
  print(sessionInfo())
})
writeLines(report, file.path(out_dir, "standard_module_preservation_report.txt"))
cat(paste(report, collapse = "\n"), "\n")
