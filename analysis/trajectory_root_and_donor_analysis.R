source(file.path(if (length(grep("^--file=", commandArgs(FALSE), value=TRUE))) dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value=TRUE)[1]), winslash="/")) else getwd(), "package_config.R"))
suppressPackageStartupMessages({
  library(monocle)
  library(igraph)
  library(Matrix)
  library(Biobase)
  library(dplyr)
  library(readr)
  library(tidyr)
  library(ggplot2)
})

set.seed(20260711)

audit_dir <- audit_root
dir.create(audit_dir, recursive = TRUE, showWarnings = FALSE)

cds_path <- required_input(
  file.path(input_root, "discovery_trajectory_graph_pseudotime_cds.rds"),
  "Discovery trajectory CellDataSet"
)
gene_audit_path <- required_input(
  file.path(audit_root, "clean_programs", "clean_D1D2_gene_level_audit.csv"),
  "Frozen clean-program gene audit"
)
dynamic_path <- required_input(
  file.path(repo_root, "data", "derived", "dynamic_program_candidates.csv"),
  "Frozen dynamic-program candidate table"
)

cds <- readRDS(cds_path)
pd <- pData(cds)
counts <- exprs(cds)
coords <- t(reducedDimS(cds))
vertex_coords <- t(reducedDimK(cds))
colnames(coords) <- c("Component_1", "Component_2")
colnames(vertex_coords) <- c("Component_1", "Component_2")

closest <- cds@auxOrderingData[["DDRTree"]][["pr_graph_cell_proj_closest_vertex"]]
closest_vertex <- paste0("Y_", as.character(closest[, 1]))
names(closest_vertex) <- rownames(closest)
closest_vertex <- closest_vertex[rownames(pd)]

sf <- pd$Size_Factor
names(sf) <- rownames(pd)
norm_log <- log1p(t(t(counts) / sf[colnames(counts)]))

score_mean <- function(genes) {
  present <- intersect(genes, rownames(norm_log))
  if (!length(present)) return(rep(NA_real_, ncol(norm_log)))
  Matrix::colMeans(norm_log[present, , drop = FALSE])
}

homeostatic_genes <- c("PI16", "APOD", "ABCA8", "IGFBP3", "WIF1")
homeostatic_score <- score_mean(homeostatic_genes)
names(homeostatic_score) <- colnames(norm_log)

clean <- read.csv(gene_audit_path, stringsAsFactors = FALSE, check.names = FALSE)
d1_clean <- unique(clean$gene[clean$clean_module == "clean_D1" & !is.na(clean$clean_module)])
d2_clean <- unique(clean$gene[clean$clean_module == "clean_D2" & !is.na(clean$clean_module)])
stopifnot(length(d1_clean) == 42L, length(d2_clean) == 975L)
d1_score <- score_mean(d1_clean)
d2_score <- score_mean(d2_clean)
names(d1_score) <- names(d2_score) <- colnames(norm_log)

mst <- minSpanningTree(cds)
if (is.null(E(mst)$weight)) {
  E(mst)$weight <- apply(ends(mst, E(mst)), 1, function(v) {
    sqrt(sum((vertex_coords[v[1], ] - vertex_coords[v[2], ])^2))
  })
}
deg <- degree(mst)

cell_meta <- data.frame(
  cell = rownames(pd),
  sample = as.character(pd$sample),
  condition = as.character(pd$condition),
  program = as.character(pd$major_fibroblast_subtype),
  closest_vertex = unname(closest_vertex),
  homeostatic_score = unname(homeostatic_score[rownames(pd)]),
  clean_D1_score = unname(d1_score[rownames(pd)]),
  clean_D2_score = unname(d2_score[rownames(pd)]),
  old_pseudotime = as.numeric(pd$Pseudotime_graph),
  stringsAsFactors = FALSE
)

vertex_audit <- cell_meta %>%
  group_by(closest_vertex) %>%
  summarise(
    n_cells = n(),
    homeostatic_fraction = mean(program == "Homeostatic/regulatory fibroblasts"),
    matrix_fraction = mean(program == "Matrix-remodeling fibroblasts"),
    ecm_fraction = mean(program == "ECM-producing fibroblasts"),
    mean_homeostatic_score = mean(homeostatic_score, na.rm = TRUE),
    median_homeostatic_score = median(homeostatic_score, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    degree = unname(deg[closest_vertex]),
    is_topological_endpoint = degree == 1,
    old_selected_root = closest_vertex == unique(pd$Root_vertex_graph)[1]
  ) %>%
  arrange(desc(is_topological_endpoint), desc(mean_homeostatic_score))

endpoint_candidates <- vertex_audit %>%
  filter(is_topological_endpoint, n_cells >= 8) %>%
  arrange(desc(mean_homeostatic_score), desc(homeostatic_fraction), desc(n_cells))
if (nrow(endpoint_candidates) < 3) {
  endpoint_candidates <- vertex_audit %>%
    filter(is_topological_endpoint) %>%
    arrange(desc(mean_homeostatic_score), desc(homeostatic_fraction), desc(n_cells))
}
candidate_roots <- head(endpoint_candidates$closest_vertex, 3)
old_root <- unique(pd$Root_vertex_graph)[1]
all_roots <- unique(c(candidate_roots, old_root))

calc_pt <- function(root) {
  vd <- distances(mst, v = root, to = V(mst), weights = E(mst)$weight)
  vd <- as.numeric(vd[1, ])
  names(vd) <- V(mst)$name
  local <- sqrt(rowSums((coords[rownames(pd), , drop = FALSE] -
                           vertex_coords[closest_vertex, , drop = FALSE])^2))
  x <- vd[closest_vertex] + local
  unname(x - min(x, na.rm = TRUE))
}

pt_mat <- sapply(all_roots, calc_pt)
colnames(pt_mat) <- all_roots
rownames(pt_mat) <- rownames(pd)

root_concordance <- expand.grid(root_a = all_roots, root_b = all_roots, stringsAsFactors = FALSE) %>%
  rowwise() %>%
  mutate(spearman_rho = suppressWarnings(cor(pt_mat[, root_a], pt_mat[, root_b], method = "spearman"))) %>%
  ungroup()

dynamic <- read.csv(dynamic_path, stringsAsFactors = FALSE, check.names = FALSE)
dyn_genes <- intersect(dynamic$gene, rownames(norm_log))
expr_fraction <- Matrix::rowMeans(counts[dyn_genes, , drop = FALSE] > 0)

direction_audit <- bind_rows(lapply(all_roots, function(root) {
  pt <- pt_mat[, root]
  rho <- apply(as.matrix(norm_log[dyn_genes, , drop = FALSE]), 1, function(x) {
    suppressWarnings(cor(x, pt, method = "spearman", use = "complete.obs"))
  })
  pval <- apply(as.matrix(norm_log[dyn_genes, , drop = FALSE]), 1, function(x) {
    suppressWarnings(cor.test(x, pt, method = "spearman", exact = FALSE)$p.value)
  })
  data.frame(
    root = root,
    gene = dyn_genes,
    original_module = dynamic$dynamic_module[match(dyn_genes, dynamic$gene)],
    original_rho = dynamic$spearman_rho[match(dyn_genes, dynamic$gene)],
    rho = rho,
    q_value = p.adjust(pval, "BH"),
    expressed_fraction = unname(expr_fraction[dyn_genes]),
    same_direction = sign(rho) == sign(dynamic$spearman_rho[match(dyn_genes, dynamic$gene)]),
    passes_original_threshold = p.adjust(pval, "BH") < 0.05 & abs(rho) > 0.20 & unname(expr_fraction[dyn_genes]) > 0.10,
    stringsAsFactors = FALSE
  )
}))

root_summary <- direction_audit %>%
  group_by(root, original_module) %>%
  summarise(
    genes = n(),
    direction_concordance = mean(same_direction, na.rm = TRUE),
    threshold_retained_fraction = mean(same_direction & passes_original_threshold, na.rm = TRUE),
    median_rho = median(rho, na.rm = TRUE),
    .groups = "drop"
  )

score_trends <- bind_rows(lapply(all_roots, function(root) {
  pt <- pt_mat[, root]
  data.frame(
    root = root,
    score = c("clean_D1", "clean_D2", "homeostatic_marker"),
    spearman_rho = c(
      cor(d1_score[rownames(pd)], pt, method = "spearman", use = "complete.obs"),
      cor(d2_score[rownames(pd)], pt, method = "spearman", use = "complete.obs"),
      cor(homeostatic_score[rownames(pd)], pt, method = "spearman", use = "complete.obs")
    ),
    stringsAsFactors = FALSE
  )
}))

sample_stats <- cell_meta %>%
  group_by(sample, condition) %>%
  summarise(
    n_cells = n(),
    median_old_pseudotime = median(old_pseudotime, na.rm = TRUE),
    median_clean_D1 = median(clean_D1_score, na.rm = TRUE),
    median_clean_D2 = median(clean_D2_score, na.rm = TRUE),
    homeostatic_fraction = mean(program == "Homeostatic/regulatory fibroblasts"),
    matrix_fraction = mean(program == "Matrix-remodeling fibroblasts"),
    ecm_fraction = mean(program == "ECM-producing fibroblasts"),
    .groups = "drop"
  )
for (root in all_roots) {
  tmp <- data.frame(cell = rownames(pd), pt = pt_mat[, root]) %>%
    left_join(cell_meta[, c("cell", "sample")], by = "cell") %>%
    group_by(sample) %>% summarise(value = median(pt), .groups = "drop")
  sample_stats[[paste0("median_pt_", root)]] <- tmp$value[match(sample_stats$sample, tmp$sample)]
}

all_labelings <- combn(seq_len(nrow(sample_stats)), 3, simplify = FALSE)
permutation_test <- function(values, labels) {
  obs <- mean(values[labels == "Keloid"]) - mean(values[labels == "Normal"])
  null <- vapply(all_labelings, function(k_idx) {
    lab <- rep("Normal", length(values)); lab[k_idx] <- "Keloid"
    mean(values[lab == "Keloid"]) - mean(values[lab == "Normal"])
  }, numeric(1))
  c(effect = obs, exact_permutation_p = mean(abs(null) >= abs(obs) - 1e-12))
}

cliffs_delta <- function(values, labels) {
  keloid <- values[labels == "Keloid"]
  normal <- values[labels == "Normal"]
  pair_sign <- as.vector(outer(keloid, normal, FUN = function(x, y) sign(x - y)))
  c(
    cliffs_delta = mean(pair_sign),
    keloid_higher_pair_fraction = mean(pair_sign > 0),
    tied_pair_fraction = mean(pair_sign == 0),
    cross_group_pair_count = length(pair_sign)
  )
}

sample_tests <- bind_rows(lapply(setdiff(names(sample_stats), c("sample", "condition", "n_cells")), function(metric) {
  x <- sample_stats[[metric]]
  out <- permutation_test(x, sample_stats$condition)
  ordinal_effect <- cliffs_delta(x, sample_stats$condition)
  data.frame(
    metric = metric,
    normal_median = median(x[sample_stats$condition == "Normal"]),
    keloid_median = median(x[sample_stats$condition == "Keloid"]),
    keloid_minus_normal_mean = unname(out["effect"]),
    wilcoxon_exact_p = suppressWarnings(wilcox.test(x ~ sample_stats$condition, exact = TRUE)$p.value),
    exact_permutation_p = unname(out["exact_permutation_p"]),
    cliffs_delta = unname(ordinal_effect["cliffs_delta"]),
    keloid_higher_pair_fraction = unname(ordinal_effect["keloid_higher_pair_fraction"]),
    tied_pair_fraction = unname(ordinal_effect["tied_pair_fraction"]),
    cross_group_pair_count = unname(ordinal_effect["cross_group_pair_count"]),
    stringsAsFactors = FALSE
  )
}))

write.csv(vertex_audit, file.path(audit_dir, "root_endpoint_condition_blind_audit.csv"), row.names = FALSE)
write.csv(data.frame(rank = seq_along(candidate_roots), root = candidate_roots), file.path(audit_dir, "condition_blind_candidate_roots.csv"), row.names = FALSE)
write.csv(root_concordance, file.path(audit_dir, "root_pseudotime_rank_concordance.csv"), row.names = FALSE)
write.csv(direction_audit, file.path(audit_dir, "root_dynamic_gene_direction_sensitivity.csv"), row.names = FALSE)
write.csv(root_summary, file.path(audit_dir, "root_dynamic_module_retention_summary.csv"), row.names = FALSE)
write.csv(score_trends, file.path(audit_dir, "root_clean_module_score_trends.csv"), row.names = FALSE)
write.csv(sample_stats, file.path(audit_dir, "discovery_donor_level_summary.csv"), row.names = FALSE)
write.csv(sample_tests, file.path(audit_dir, "discovery_donor_level_tests.csv"), row.names = FALSE)
write.csv(cell_meta, file.path(audit_dir, "discovery_cell_metadata_sensitivity.csv"), row.names = FALSE)

metadata_columns <- data.frame(column = colnames(pd), class = vapply(pd, function(x) paste(class(x), collapse = ";"), character(1)))
write.csv(metadata_columns, file.path(audit_dir, "trajectory_object_metadata_columns.csv"), row.names = FALSE)
writeLines(capture.output({
  cat("CDS:", cds_path, "\n")
  cat("cells:", ncol(cds), "genes:", nrow(cds), "\n")
  cat("old root:", old_root, "\n")
  cat("condition-blind endpoint candidates:", paste(candidate_roots, collapse = ", "), "\n")
  print(endpoint_candidates)
  cat("\nRoot module sensitivity:\n")
  print(root_summary)
  cat("\nDonor-level tests:\n")
  print(sample_tests)
  cat("\nSession info:\n")
  print(sessionInfo())
}), file.path(audit_dir, "phase1_reanalysis_report.txt"))

cat("phase1 complete\n")
print(endpoint_candidates)
print(root_summary)
print(sample_tests)
