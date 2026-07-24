source(file.path(if (length(grep("^--file=", commandArgs(FALSE), value=TRUE))) dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value=TRUE)[1]), winslash="/")) else getwd(), "package_config.R"))
suppressPackageStartupMessages(library(Seurat))

obj <- readRDS(external_object(
  "KFT_GSE191067_OBJECT",
  file.path(repo_root, "data", "frozen_inputs", "GSE191067_standard_reclustered_scored.rds"),
  "Processed GSE191067 fibroblast Seurat object"
))

cat("class:", class(obj), "\n")
cat("dimensions:", nrow(obj), "genes x", ncol(obj), "cells\n")
cat("assays:", paste(Assays(obj), collapse = ", "), "\n")
cat("reductions:", paste(Reductions(obj), collapse = ", "), "\n")
cat("graphs:", paste(names(obj@graphs), collapse = ", "), "\n")
cat("metadata columns:\n")
print(colnames(obj@meta.data))

for (field in c("sample", "sample_id", "condition", "original_validation_cluster",
                "seurat_clusters", "validation_program_3")) {
  if (field %in% colnames(obj@meta.data)) {
    cat("\n", field, ":\n", sep = "")
    print(table(obj@meta.data[[field]], useNA = "ifany"))
  }
}

if (all(c("original_validation_cluster", "validation_program_3") %in% colnames(obj@meta.data))) {
  cat("\ncluster-to-program mapping:\n")
  print(table(obj$original_validation_cluster, obj$validation_program_3, useNA = "ifany"))
}

cat("\ncommands:\n")
for (nm in names(obj@commands)) {
  cmd <- obj@commands[[nm]]
  cat("\n[", nm, "]\n", sep = "")
  print(cmd@params)
}

cat("\nsessionInfo:\n")
print(sessionInfo())
