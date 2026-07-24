source(file.path(if (length(grep("^--file=", commandArgs(FALSE), value=TRUE))) dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value=TRUE)[1]), winslash="/")) else getwd(), "package_config.R"))
suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(dplyr)
})

out_dir <- file.path(audit_root, "GSE191067_standard_validation")
obj <- readRDS(external_object(
  "KFT_GSE191067_OBJECT",
  file.path(repo_root, "data", "frozen_inputs", "GSE191067_standard_reclustered_scored.rds"),
  "Processed GSE191067 fibroblast Seurat object"
))
DefaultAssay(obj) <- "RNA"
x <- GetAssayData(obj, assay = "RNA", layer = "data")

marker_sets <- list(
  `Homeostatic/regulatory fibroblasts` = c("PI16", "DPT", "APOD", "CFD", "IGFBP3", "WIF1", "BMP7", "ABCA8", "ADH1B", "C7"),
  `Matrix-remodeling fibroblasts` = c("NGFR", "THBS4", "HAS2", "TIMP1", "COL6A5", "F13A1", "LOXL3", "CTSC", "FOXS1", "CHST1"),
  `ECM-producing fibroblasts` = c("COL1A1", "COL1A2", "COL3A1", "COL5A1", "COMP", "ASPN", "TAGLN", "CTHRC1", "LRRC15", "POSTN", "SPARC", "SERPINH1")
)
markers <- unique(unlist(marker_sets))
programs <- levels(obj$standard_program)

stats <- do.call(rbind, lapply(markers, function(gene) {
  if (!gene %in% rownames(x)) return(NULL)
  values <- as.numeric(x[gene, ])
  z <- as.numeric(scale(values))
  z[!is.finite(z)] <- 0
  do.call(rbind, lapply(programs, function(program) {
    keep <- obj$standard_program == program
    data.frame(
      gene = gene,
      annotation_marker_program = names(marker_sets)[vapply(marker_sets, function(gs) gene %in% gs, logical(1))][1],
      standard_program = program,
      n_cells = sum(keep),
      pct_expressing = mean(values[keep] > 0) * 100,
      mean_log_normalized_expression = mean(values[keep]),
      mean_gene_z = mean(z[keep]),
      stringsAsFactors = FALSE
    )
  }))
}))

write.csv(stats, file.path(out_dir, "standard_marker_dotplot_statistics.csv"), row.names = FALSE)
cat("standard validation summaries complete\n")
