source(file.path(if (length(grep("^--file=", commandArgs(FALSE), value=TRUE))) dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value=TRUE)[1]), winslash="/")) else getwd(), "package_config.R"))
suppressPackageStartupMessages({
  library(monocle)
  library(Matrix)
  library(Biobase)
  library(matrixStats)
  library(slingshot)
  library(SingleCellExperiment)
  library(dplyr)
})

set.seed(20260711)
out_dir <- file.path(audit_root, "branch_trajectory")
cds <- readRDS(file.path(input_root, "discovery_trajectory_graph_pseudotime_cds.rds"))
pd <- pData(cds)
counts <- round(as.matrix(exprs(cds)))
sf <- pd$Size_Factor
names(sf) <- rownames(pd)
log_norm <- log1p(t(t(counts) / sf[colnames(counts)]))
rv <- matrixStats::rowVars(as.matrix(log_norm))
names(rv) <- rownames(log_norm)
pca_genes <- names(sort(rv, decreasing = TRUE))[seq_len(min(2000, length(rv)))]
pca <- prcomp(t(as.matrix(log_norm[pca_genes, , drop = FALSE])), center = TRUE, scale. = TRUE, rank. = 20)

sce <- SingleCellExperiment(assays = list(counts = counts[, rownames(pd), drop = FALSE]))
reducedDims(sce)$PCA <- pca$x[, seq_len(min(20, ncol(pca$x))), drop = FALSE]
colData(sce)$program <- factor(pd$major_fibroblast_subtype)
start_cluster <- "Homeostatic/regulatory fibroblasts"
end_clusters <- c("Matrix-remodeling fibroblasts", "ECM-producing fibroblasts")

forced <- slingshot(sce, clusterLabels = "program", reducedDim = "PCA", start.clus = start_cluster, end.clus = end_clusters, stretch = 0)
unsup <- slingshot(sce, clusterLabels = "program", reducedDim = "PCA", start.clus = start_cluster, stretch = 0)

extract_long <- function(obj, mode) {
  pt <- slingPseudotime(obj)
  w <- slingCurveWeights(obj)
  bind_rows(lapply(seq_len(ncol(pt)), function(j) data.frame(
    cell = rownames(pd), mode = mode, lineage = paste0("lineage", j),
    pseudotime = pt[, j], weight = w[, j], stringsAsFactors = FALSE
  )))
}
sling_long <- bind_rows(extract_long(forced, "forced_two_terminal"), extract_long(unsup, "start_only_unsupervised"))
old_pt <- setNames(pd$Pseudotime_graph, rownames(pd))
concordance <- sling_long %>%
  filter(is.finite(pseudotime), weight > 0.5) %>%
  group_by(mode, lineage) %>%
  summarise(n_cells = n(), spearman_with_DDRTree_graph_pseudotime = cor(pseudotime, old_pt[cell], method = "spearman"), .groups = "drop")

write.csv(sling_long, file.path(out_dir, "slingshot_PCA_pseudotime_long.csv"), row.names = FALSE)
write.csv(concordance, file.path(out_dir, "slingshot_DDRTree_pseudotime_concordance.csv"), row.names = FALSE)
writeLines(capture.output({
  cat("Forced Slingshot lineages:\n"); print(slingLineages(forced))
  cat("Start-only Slingshot lineages:\n"); print(slingLineages(unsup))
  cat("Concordance:\n"); print(concordance)
  print(sessionInfo())
}), file.path(out_dir, "slingshot_cross_method_sensitivity_report.txt"))
cat("slingshot reconstruction complete\n")
print(concordance)
