script_dir <- if (length(grep("^--file=", commandArgs(FALSE), value = TRUE))) dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]), winslash = "/")) else getwd()
source(file.path(script_dir, "..", "analysis", "package_config.R"))
out <- file.path(package_root, "environment", "R_sessionInfo.txt")
pkgs <- c("Seurat", "SeuratObject", "harmony", "monocle", "DDRTree", "slingshot",
          "tradeSeq", "mgcv", "WGCNA", "clusterProfiler", "org.Hs.eg.db",
          "tximport", "pROC", "matrixStats", "dplyr", "ggplot2", "msigdbr",
          "scDblFinder")
versions <- vapply(pkgs, function(pkg) {
  if (requireNamespace(pkg, quietly = TRUE)) as.character(packageVersion(pkg)) else NA_character_
}, character(1))
lines <- c(
  paste("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "",
  "Package versions:",
  paste(names(versions), versions, sep = ": "),
  "",
  "Full sessionInfo():",
  capture.output(sessionInfo())
)
writeLines(lines, out, useBytes = TRUE)
