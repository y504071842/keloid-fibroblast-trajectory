source(file.path(if (length(grep("^--file=", commandArgs(FALSE), value = TRUE))) dirname(chartr("\\", "/", sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))) else getwd(), "package_config.R"))
suppressPackageStartupMessages({
  library(monocle)
  library(igraph)
  library(Matrix)
  library(Biobase)
  library(matrixStats)
  library(tradeSeq)
  library(slingshot)
  library(SingleCellExperiment)
  library(dplyr)
})

set.seed(20260711)
out_dir <- file.path(audit_root, "branch_trajectory")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
cds <- readRDS(file.path(input_root, "discovery_trajectory_graph_pseudotime_cds.rds"))
pd <- pData(cds)
counts_all <- round(as.matrix(exprs(cds)))
coords <- t(reducedDimS(cds))
vertex_coords <- t(reducedDimK(cds))

closest <- cds@auxOrderingData[["DDRTree"]][["pr_graph_cell_proj_closest_vertex"]]
closest_vertex <- paste0("Y_", as.character(closest[, 1]))
names(closest_vertex) <- rownames(closest)
closest_vertex <- closest_vertex[rownames(pd)]

mst <- minSpanningTree(cds)
if (is.null(E(mst)$weight)) {
  E(mst)$weight <- apply(ends(mst, E(mst)), 1, function(v) {
    sqrt(sum((vertex_coords[v[1], ] - vertex_coords[v[2], ])^2))
  })
}

root <- "Y_114"
terminals <- c(mixed_remodeling_branch = "Y_56", ecm_dominant_branch = "Y_1")
paths <- lapply(terminals, function(x) V(mst)$name[unlist(shortest_paths(mst, from = root, to = x, weights = E(mst)$weight)$vpath)])
vertex_distance <- distances(mst, v = root, to = V(mst), weights = E(mst)$weight)
vertex_distance <- as.numeric(vertex_distance[1, ])
names(vertex_distance) <- V(mst)$name
local_distance <- sqrt(rowSums((coords[rownames(pd), , drop = FALSE] - vertex_coords[closest_vertex, , drop = FALSE])^2))

n <- nrow(pd)
L <- length(paths)
pseudotime <- matrix(0, nrow = n, ncol = L, dimnames = list(rownames(pd), names(paths)))
cell_weights <- matrix(0, nrow = n, ncol = L, dimnames = list(rownames(pd), names(paths)))
for (j in seq_along(paths)) {
  member <- closest_vertex %in% paths[[j]]
  pseudotime[member, j] <- vertex_distance[closest_vertex[member]] + local_distance[member]
  cell_weights[member, j] <- 1
}
shared <- rowSums(cell_weights) > 1
cell_weights[shared, ] <- cell_weights[shared, , drop = FALSE] / rowSums(cell_weights[shared, , drop = FALSE])
if (any(rowSums(cell_weights) == 0)) stop("Some cells were not assigned to either graph lineage")
for (j in seq_len(L)) {
  member <- cell_weights[, j] > 0
  r <- range(pseudotime[member, j], na.rm = TRUE)
  pseudotime[member, j] <- (pseudotime[member, j] - r[1]) / diff(r)
}

dynamic <- read.csv(required_input(
  file.path(repo_root, "data", "derived", "dynamic_program_candidates.csv"),
  "Frozen dynamic-program candidate table"
), stringsAsFactors = FALSE)
genes <- intersect(dynamic$gene, rownames(counts_all))
counts <- counts_all[genes, , drop = FALSE]
counts <- counts[rowSums(counts) > 0, , drop = FALSE]
genes <- rownames(counts)
size_factor <- as.numeric(pd$Size_Factor)
if (any(!is.finite(size_factor)) || any(size_factor <= 0)) {
  stop("Invalid Monocle2 Size_Factor values in trajectory cells")
}
log_norm_counts <- log1p(sweep(counts, 2, size_factor, FUN = "/"))

lineage_rho <- do.call(rbind, lapply(seq_len(L), function(j) {
  member <- cell_weights[, j] > 0
  rho <- apply(log_norm_counts[, member, drop = FALSE], 1, function(x) cor(x, pseudotime[member, j], method = "spearman"))
  pval <- apply(log_norm_counts[, member, drop = FALSE], 1, function(x) suppressWarnings(cor.test(x, pseudotime[member, j], method = "spearman", exact = FALSE)$p.value))
  data.frame(
    gene = genes,
    lineage = colnames(pseudotime)[j],
    n_cells = sum(member),
    spearman_rho = rho,
    p_value = pval,
    fdr = p.adjust(pval, "BH"),
    stringsAsFactors = FALSE
  )
}))
write.csv(lineage_rho, file.path(out_dir, "lineage_specific_spearman_dynamic_tests.csv"), row.names = FALSE)

rho_wide <- reshape(lineage_rho[, c("gene", "lineage", "spearman_rho", "fdr")], idvar = "gene", timevar = "lineage", direction = "wide")
m_rho <- rho_wide[["spearman_rho.mixed_remodeling_branch"]]
e_rho <- rho_wide[["spearman_rho.ecm_dominant_branch"]]
m_fdr <- rho_wide[["fdr.mixed_remodeling_branch"]]
e_fdr <- rho_wide[["fdr.ecm_dominant_branch"]]
rho_wide$branch_pattern <- ifelse(m_fdr < 0.05 & e_fdr < 0.05 & m_rho < -0.20 & e_rho < -0.20, "common_decreasing",
  ifelse(m_fdr < 0.05 & e_fdr < 0.05 & m_rho > 0.20 & e_rho > 0.20, "common_increasing",
    ifelse(m_fdr < 0.05 & abs(m_rho) > 0.20 & !(e_fdr < 0.05 & abs(e_rho) > 0.20), "mixed_remodeling_branch_biased",
      ifelse(e_fdr < 0.05 & abs(e_rho) > 0.20 & !(m_fdr < 0.05 & abs(m_rho) > 0.20), "ecm_dominant_branch_biased", "other"))))
rho_wide$original_module <- dynamic$dynamic_module[match(rho_wide$gene, dynamic$gene)]
write.csv(rho_wide, file.path(out_dir, "lineage_dynamic_gene_classification.csv"), row.names = FALSE)

lineage_meta <- data.frame(
  cell = rownames(pd), sample = pd$sample, condition = pd$condition,
  program = pd$major_fibroblast_subtype, closest_vertex = closest_vertex,
  pseudotime_mixed = pseudotime[, 1], pseudotime_ecm_dominant = pseudotime[, 2],
  weight_mixed = cell_weights[, 1], weight_ecm_dominant = cell_weights[, 2],
  stringsAsFactors = FALSE
)
write.csv(lineage_meta, file.path(out_dir, "graph_two_lineage_cell_metadata.csv"), row.names = FALSE)
write.csv(do.call(rbind, lapply(names(paths), function(nm) data.frame(lineage = nm, vertex_order = seq_along(paths[[nm]]), vertex = paths[[nm]]))), file.path(out_dir, "graph_two_lineage_vertex_paths.csv"), row.names = FALSE)

cat("Fitting two-lineage tradeSeq model for", nrow(counts), "genes and", ncol(counts), "cells\n")
donor <- factor(pd$sample)
# tradeSeq uses y ~ -1 + U + smooth + offset, so U must include an intercept.
donor_design <- model.matrix(~ donor)
fit <- fitGAM(
  counts = counts,
  pseudotime = pseudotime,
  cellWeights = cell_weights,
  U = donor_design,
  nknots = 6,
  verbose = TRUE,
  parallel = FALSE
)

assoc_global <- as.data.frame(associationTest(fit, global = TRUE, lineages = FALSE))
assoc_global$gene <- rownames(assoc_global)
assoc_global$fdr <- p.adjust(assoc_global$pvalue, "BH")
assoc_lineage <- as.data.frame(associationTest(fit, global = FALSE, lineages = TRUE))
assoc_lineage$gene <- rownames(assoc_lineage)
pattern <- as.data.frame(patternTest(fit, nPoints = 50))
pattern$gene <- rownames(pattern)
pattern$fdr <- p.adjust(pattern$pvalue, "BH")
diff_end <- as.data.frame(diffEndTest(fit))
diff_end$gene <- rownames(diff_end)
diff_end$fdr <- p.adjust(diff_end$pvalue, "BH")
early <- as.data.frame(earlyDETest(fit, knots = c(1, 2)))
early$gene <- rownames(early)
early$fdr <- p.adjust(early$pvalue, "BH")

write.csv(assoc_global, file.path(out_dir, "tradeSeq_two_lineage_association_global.csv"), row.names = FALSE)
write.csv(assoc_lineage, file.path(out_dir, "tradeSeq_two_lineage_association_by_lineage.csv"), row.names = FALSE)
write.csv(pattern, file.path(out_dir, "tradeSeq_two_lineage_patternTest.csv"), row.names = FALSE)
write.csv(diff_end, file.path(out_dir, "tradeSeq_two_lineage_diffEndTest.csv"), row.names = FALSE)
write.csv(early, file.path(out_dir, "tradeSeq_two_lineage_earlyDETest.csv"), row.names = FALSE)
test_status <- dplyr::bind_rows(
  data.frame(test = "association", input_genes = nrow(assoc_global), estimable_genes = sum(is.finite(assoc_global$fdr)), fdr_lt_0_05 = sum(is.finite(assoc_global$fdr) & assoc_global$fdr < 0.05)),
  data.frame(test = "pattern", input_genes = nrow(pattern), estimable_genes = sum(is.finite(pattern$fdr)), fdr_lt_0_05 = sum(is.finite(pattern$fdr) & pattern$fdr < 0.05)),
  data.frame(test = "terminal", input_genes = nrow(diff_end), estimable_genes = sum(is.finite(diff_end$fdr)), fdr_lt_0_05 = sum(is.finite(diff_end$fdr) & diff_end$fdr < 0.05)),
  data.frame(test = "early", input_genes = nrow(early), estimable_genes = sum(is.finite(early$fdr)), fdr_lt_0_05 = sum(is.finite(early$fdr) & early$fdr < 0.05))
) %>%
  mutate(reported_fraction = ifelse(estimable_genes > 0, fdr_lt_0_05 / estimable_genes, NA_real_))
write.csv(test_status, file.path(out_dir, "tradeSeq_test_estimability_summary.csv"), row.names = FALSE)
write.csv(data.frame(
  parameter = c(
    "lineage_spearman_expression", "lineage_spearman_multiple_testing",
    "tradeSeq_input", "tradeSeq_family", "tradeSeq_nknots", "tradeSeq_nknots_basis",
    "tradeSeq_donor_design", "tradeSeq_library_size_adjustment", "random_seed"
  ),
  value = c(
    "log1p(raw_count / Monocle2 Size_Factor)",
    "Benjamini-Hochberg within each lineage across 1,118 dynamic genes",
    "integer raw counts", "negative binomial", "6",
    "prespecified tradeSeq default; evaluateK sensitivity compared 3-8 knots in 500 fixed-seed sampled genes",
    "model.matrix(~ donor), including intercept and five donor indicators",
    "default TMM-derived effective-library-size offset", "20260711"
  ),
  stringsAsFactors = FALSE
), file.path(out_dir, "branch_analysis_parameters.csv"), row.names = FALSE)

# Slingshot cross-method sensitivity reconstruction on expression-derived principal components.
sf <- pd$Size_Factor
names(sf) <- rownames(pd)
log_norm <- log1p(t(t(counts_all) / sf[colnames(counts_all)]))
rv <- matrixStats::rowVars(as.matrix(log_norm))
names(rv) <- rownames(log_norm)
pca_genes <- names(sort(rv, decreasing = TRUE))[seq_len(min(2000, length(rv)))]
pca <- prcomp(t(as.matrix(log_norm[pca_genes, , drop = FALSE])), center = TRUE, scale. = TRUE, rank. = 20)
clusters <- factor(pd$major_fibroblast_subtype)
start_cluster <- "Homeostatic/regulatory fibroblasts"
end_clusters <- c("Matrix-remodeling fibroblasts", "ECM-producing fibroblasts")

sce <- SingleCellExperiment(assays = list(counts = counts_all[, rownames(pd), drop = FALSE]))
reducedDims(sce)$PCA <- pca$x[, seq_len(min(20, ncol(pca$x))), drop = FALSE]
colData(sce)$program <- clusters
sce_forced <- slingshot(sce, clusterLabels = "program", reducedDim = "PCA", start.clus = start_cluster, end.clus = end_clusters, stretch = 0)
sce_unsup <- slingshot(sce, clusterLabels = "program", reducedDim = "PCA", start.clus = start_cluster, stretch = 0)

extract_slingshot <- function(obj, mode) {
  pt <- slingPseudotime(obj)
  w <- slingCurveWeights(obj)
  out <- data.frame(cell = rownames(pd), mode = mode, stringsAsFactors = FALSE)
  for (j in seq_len(ncol(pt))) {
    out[[paste0("pseudotime_lineage", j)]] <- pt[, j]
    out[[paste0("weight_lineage", j)]] <- w[, j]
  }
  out
}
sling_meta <- dplyr::bind_rows(
  extract_slingshot(sce_forced, "forced_two_terminal"),
  extract_slingshot(sce_unsup, "start_only_unsupervised")
)
write.csv(sling_meta, file.path(out_dir, "slingshot_independent_pca_pseudotime.csv"), row.names = FALSE)
writeLines(capture.output({
  cat("Forced Slingshot lineages:\n"); print(slingLineages(sce_forced))
  cat("Unsupervised Slingshot lineages:\n"); print(slingLineages(sce_unsup))
  cat("Graph lineage pattern counts:\n"); print(table(rho_wide$branch_pattern, rho_wide$original_module))
  cat("Lineage Spearman expression: log1p(raw_count / Monocle2 Size_Factor)\n")
  cat("tradeSeq donor design columns:", paste(colnames(donor_design), collapse = ", "), "\n")
  cat("tradeSeq global FDR < 0.05:", sum(assoc_global$fdr < 0.05, na.rm = TRUE), "\n")
  cat("patternTest FDR < 0.05:", sum(pattern$fdr < 0.05, na.rm = TRUE), "\n")
  cat("diffEndTest FDR < 0.05:", sum(diff_end$fdr < 0.05, na.rm = TRUE), "\n")
  cat("earlyDETest FDR < 0.05:", sum(early$fdr < 0.05, na.rm = TRUE), "\n")
  print(sessionInfo())
}), file.path(out_dir, "branch_trajectory_report.txt"))

cat("branch trajectory reanalysis complete\n")
