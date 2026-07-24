source(file.path(if (length(grep("^--file=", commandArgs(FALSE), value = TRUE))) dirname(chartr("\\", "/", sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))) else getwd(), "package_config.R"))
suppressPackageStartupMessages({
  library(monocle)
  library(igraph)
  library(Matrix)
  library(Biobase)
  library(tradeSeq)
})

set.seed(20260717)
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
vertex_distance <- as.numeric(distances(mst, v = root, to = V(mst), weights = E(mst)$weight)[1, ])
names(vertex_distance) <- V(mst)$name
local_distance <- sqrt(rowSums((coords[rownames(pd), , drop = FALSE] - vertex_coords[closest_vertex, , drop = FALSE])^2))

pseudotime <- matrix(0, nrow = nrow(pd), ncol = length(paths), dimnames = list(rownames(pd), names(paths)))
cell_weights <- matrix(0, nrow = nrow(pd), ncol = length(paths), dimnames = list(rownames(pd), names(paths)))
for (j in seq_along(paths)) {
  member <- closest_vertex %in% paths[[j]]
  pseudotime[member, j] <- vertex_distance[closest_vertex[member]] + local_distance[member]
  cell_weights[member, j] <- 1
}
shared <- rowSums(cell_weights) > 1
cell_weights[shared, ] <- cell_weights[shared, , drop = FALSE] / rowSums(cell_weights[shared, , drop = FALSE])
if (any(rowSums(cell_weights) == 0)) stop("Some cells were not assigned to either graph lineage")
for (j in seq_len(ncol(pseudotime))) {
  member <- cell_weights[, j] > 0
  rr <- range(pseudotime[member, j], na.rm = TRUE)
  pseudotime[member, j] <- (pseudotime[member, j] - rr[1]) / diff(rr)
}

dynamic <- read.csv(file.path(audit_root, "clean_programs", "clean_D1D2_gene_level_audit.csv"), stringsAsFactors = FALSE)
genes <- intersect(dynamic$gene, rownames(counts_all))
counts <- counts_all[genes, , drop = FALSE]
counts <- counts[rowSums(counts) > 0, , drop = FALSE]
donor_design <- model.matrix(~ factor(pd$sample))

k_values <- 3:8
n_genes <- min(500L, nrow(counts))
aic <- evaluateK(
  counts = counts,
  pseudotime = pseudotime,
  cellWeights = cell_weights,
  U = donor_design,
  nGenes = n_genes,
  k = k_values,
  plot = FALSE,
  verbose = TRUE,
  parallel = FALSE
)

colnames(aic) <- as.character(k_values)
gene_min_aic <- apply(aic, 1, min, na.rm = TRUE)
aic_summary <- data.frame(
  nknots = k_values,
  evaluable_gene_n = colSums(is.finite(aic)),
  median_AIC = apply(aic, 2, median, na.rm = TRUE),
  mean_AIC = colMeans(aic, na.rm = TRUE),
  genes_within_2_AIC_of_gene_minimum = vapply(seq_along(k_values), function(i) {
    sum(aic[, i] <= gene_min_aic + 2, na.rm = TRUE)
  }, integer(1)),
  fraction_within_2_AIC_of_gene_minimum = vapply(seq_along(k_values), function(i) {
    mean(aic[, i] <= gene_min_aic + 2, na.rm = TRUE)
  }, numeric(1))
)

write.csv(aic_summary, file.path(out_dir, "tradeSeq_knot_sensitivity_summary.csv"), row.names = FALSE)
write.csv(data.frame(gene = rownames(aic), aic, check.names = FALSE), file.path(out_dir, "tradeSeq_knot_sensitivity_gene_AIC.csv"), row.names = FALSE)
write.csv(data.frame(
  parameter = c("seed", "candidate_knots", "evaluated_gene_n", "gene_sampling", "donor_covariates", "selection_role"),
  value = c(20260717, "3,4,5,6,7,8", n_genes, "tradeSeq evaluateK random subset with fixed seed", TRUE,
            "sensitivity assessment of the prespecified tradeSeq default nknots=6")
), file.path(out_dir, "tradeSeq_knot_sensitivity_parameters.csv"), row.names = FALSE)
writeLines(capture.output({
  cat("tradeSeq knot sensitivity\n")
  cat("tradeSeq default nknots=6 was prespecified; evaluateK assessed k=3:8 on", n_genes, "genes.\n")
  print(aic_summary)
  print(sessionInfo())
}), file.path(out_dir, "tradeSeq_knot_sensitivity_report.txt"))

print(aic_summary)
