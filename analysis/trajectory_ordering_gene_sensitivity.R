source(file.path(if (length(grep("^--file=", commandArgs(FALSE), value=TRUE))) dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value=TRUE)[1]), winslash="/")) else getwd(), "package_config.R"))
suppressPackageStartupMessages({
  library(monocle)
  library(Biobase)
  library(Matrix)
  library(igraph)
  library(dplyr)
})

set.seed(20260711)
out_dir <- file.path(audit_root, "ordering_gene_sensitivity")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

patch_dplyr_underscored <- function() {
  select_compat <- function(.data, ...) { dots <- list(...); exprs <- lapply(dots, function(x) if (is.character(x)) rlang::parse_expr(x) else x); names(exprs) <- names(dots); dplyr::select(.data, !!!exprs) }
  mutate_compat <- function(.data, ...) { dots <- list(...); exprs <- lapply(dots, function(x) if (is.character(x)) rlang::parse_expr(x) else x); names(exprs) <- names(dots); dplyr::mutate(.data, !!!exprs) }
  for (env in list(asNamespace("dplyr"), parent.env(asNamespace("monocle")), asNamespace("monocle"))) {
    for (nm in c("select_", "mutate_")) {
      fn <- if (nm == "select_") select_compat else mutate_compat
      if (exists(nm, envir = env, inherits = FALSE)) { try(unlockBinding(nm, env), silent = TRUE); assign(nm, fn, envir = env); try(lockBinding(nm, env), silent = TRUE) }
    }
  }
}
patch_dplyr_underscored()

base_cds <- readRDS(file.path(input_root, "discovery_trajectory_programDE_cds.rds"))
current_cds <- readRDS(file.path(input_root, "discovery_trajectory_graph_pseudotime_cds.rds"))
pd <- pData(base_cds)
counts <- exprs(base_cds)
sf <- pd$Size_Factor; names(sf) <- rownames(pd)
norm_log <- log1p(t(t(counts) / sf[colnames(counts)]))
current_pt <- setNames(pData(current_cds)$Pseudotime_graph, rownames(pData(current_cds)))

ordering_gene_table <- read.csv(
  file.path(audit_root, "ordering_gene_sensitivity", "ordering_gene_sets.csv"),
  stringsAsFactors = FALSE
)
if (all(c("ordering_set", "gene") %in% colnames(ordering_gene_table))) {
  program_de <- ordering_gene_table$gene[ordering_gene_table$ordering_set == "program_DE_all1586"]
} else if ("ordering_gene" %in% colnames(ordering_gene_table)) {
  program_de <- ordering_gene_table$ordering_gene
} else {
  stop("ordering_gene_sets.csv must contain ordering_set/gene or ordering_gene columns")
}
program_de <- unique(program_de[!is.na(program_de) & nzchar(program_de)])
stopifnot(length(program_de) == 1586L)
rv <- matrixStats::rowVars(as.matrix(norm_log)); names(rv) <- rownames(norm_log)
hvg <- names(sort(rv, decreasing = TRUE))
ordering_sets <- list(
  program_DE_top500 = head(program_de, 500),
  program_DE_top1000 = head(program_de, 1000),
  program_DE_all1586 = program_de,
  variance_top1000 = head(hvg, 1000),
  variance_top2000 = head(hvg, 2000)
)

homeostatic_genes <- intersect(c("PI16", "APOD", "ABCA8", "IGFBP3", "WIF1"), rownames(norm_log))
homeo_score <- Matrix::colMeans(norm_log[homeostatic_genes, , drop = FALSE])
names(homeo_score) <- colnames(norm_log)
dynamic <- read.csv(file.path(audit_root, "clean_programs", "clean_D1D2_gene_level_audit.csv"), stringsAsFactors = FALSE)
dyn_genes <- intersect(dynamic$gene, rownames(norm_log))

run_one <- function(name, genes) {
  cat("Running", name, "with", length(genes), "ordering genes\n")
  cds <- setOrderingFilter(base_cds, intersect(genes, rownames(base_cds)))
  set.seed(20260711)
  cds <- reduceDimension(cds, max_components = 2, method = "DDRTree")
  coords <- t(reducedDimS(cds)); vcoords <- t(reducedDimK(cds))
  closest <- cds@auxOrderingData[["DDRTree"]][["pr_graph_cell_proj_closest_vertex"]]
  cv <- paste0("Y_", as.character(closest[, 1])); names(cv) <- rownames(closest); cv <- cv[rownames(pd)]
  mst <- minSpanningTree(cds)
  if (is.null(E(mst)$weight)) E(mst)$weight <- apply(ends(mst, E(mst)), 1, function(v) sqrt(sum((vcoords[v[1], ] - vcoords[v[2], ])^2)))
  deg <- degree(mst)
  va <- data.frame(cell = rownames(pd), vertex = cv, homeo = homeo_score[rownames(pd)], program = pd$major_fibroblast_subtype) %>%
    group_by(vertex) %>% summarise(n_cells = n(), mean_homeostatic_score = mean(homeo), homeostatic_fraction = mean(program == "Homeostatic/regulatory fibroblasts"), .groups = "drop") %>%
    mutate(degree = unname(deg[vertex]), endpoint = degree == 1) %>%
    arrange(desc(endpoint), desc(mean_homeostatic_score), desc(homeostatic_fraction))
  cand <- va %>% filter(endpoint, n_cells >= 8)
  if (!nrow(cand)) cand <- va %>% filter(endpoint)
  root <- cand$vertex[1]
  vd <- distances(mst, v = root, to = V(mst), weights = E(mst)$weight); vd <- as.numeric(vd[1, ]); names(vd) <- V(mst)$name
  local <- sqrt(rowSums((coords[rownames(pd), , drop = FALSE] - vcoords[cv, , drop = FALSE])^2))
  pt <- unname(vd[cv] + local); pt <- pt - min(pt)
  names(pt) <- rownames(pd)
  rho <- apply(as.matrix(norm_log[dyn_genes, , drop = FALSE]), 1, function(x) cor(x, pt, method = "spearman"))
  direction <- sign(rho) == sign(dynamic$spearman_rho[match(dyn_genes, dynamic$gene)])
  summary <- data.frame(
    ordering_set = name, ordering_gene_n = length(intersect(genes, rownames(base_cds))),
    root_vertex_internal = root, graph_endpoints = sum(deg == 1), graph_branch_vertices = sum(deg > 2),
    spearman_with_current_graph_pseudotime = cor(pt[rownames(pd)], current_pt[rownames(pd)], method = "spearman"),
    dynamic_gene_direction_concordance = mean(direction),
    D1_direction_concordance = mean(direction[dynamic$dynamic_module[match(dyn_genes, dynamic$gene)] == "D1_early_decreasing"]),
    D2_direction_concordance = mean(direction[dynamic$dynamic_module[match(dyn_genes, dynamic$gene)] == "D2_late_increasing"]),
    stringsAsFactors = FALSE
  )
  list(summary = summary, cell = data.frame(cell = rownames(pd), ordering_set = name, pseudotime = pt[rownames(pd)], program = pd$major_fibroblast_subtype, sample = pd$sample, condition = pd$condition), root = cbind(ordering_set = name, va))
}

res <- lapply(names(ordering_sets), function(nm) run_one(nm, ordering_sets[[nm]]))
summary <- bind_rows(lapply(res, `[[`, "summary"))
cell <- bind_rows(lapply(res, `[[`, "cell"))
root <- bind_rows(lapply(res, `[[`, "root"))
write.csv(summary, file.path(out_dir, "ordering_gene_trajectory_sensitivity_summary.csv"), row.names = FALSE)
write.csv(cell, file.path(out_dir, "ordering_gene_trajectory_sensitivity_cell_pseudotime.csv"), row.names = FALSE)
write.csv(root, file.path(out_dir, "ordering_gene_trajectory_condition_blind_root_audit.csv"), row.names = FALSE)
write.csv(do.call(rbind, lapply(names(ordering_sets), function(nm) data.frame(ordering_set = nm, rank = seq_along(ordering_sets[[nm]]), gene = ordering_sets[[nm]]))), file.path(out_dir, "ordering_gene_sets.csv"), row.names = FALSE)
writeLines(capture.output({ print(summary); print(sessionInfo()) }), file.path(out_dir, "ordering_gene_sensitivity_report.txt"))
cat("ordering-gene sensitivity complete\n"); print(summary)
