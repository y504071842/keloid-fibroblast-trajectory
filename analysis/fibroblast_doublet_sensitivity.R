source(file.path(if (length(grep("^--file=", commandArgs(FALSE), value=TRUE))) dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value=TRUE)[1]), winslash="/")) else getwd(), "package_config.R"))
suppressPackageStartupMessages({ library(monocle); library(Biobase); library(Matrix); library(dplyr) })
set.seed(20260711)
out_dir <- file.path(audit_root, "QC_Harmony")
cds <- readRDS(file.path(input_root, "discovery_trajectory_graph_pseudotime_cds.rds"))
pd <- pData(cds); counts <- exprs(cds); sf <- pd$Size_Factor; names(sf) <- rownames(pd)
norm <- log1p(t(t(counts) / sf[colnames(counts)]))
dbl <- read.csv(file.path(out_dir, "candidate_fibroblast_scDblFinder_audit.csv"), stringsAsFactors = FALSE)
predicted <- dbl$cell[dbl$scDblFinder_class == "doublet"]
keep <- !rownames(pd) %in% predicted
clean <- read.csv(file.path(audit_root, "clean_programs", "clean_D1D2_gene_level_audit.csv"), stringsAsFactors = FALSE)
d1 <- unique(clean$gene[clean$clean_module == "clean_D1" & !is.na(clean$clean_module)])
d2 <- unique(clean$gene[clean$clean_module == "clean_D2" & !is.na(clean$clean_module)])
stopifnot(length(d1) == 42L, length(d2) == 975L)
score <- function(g) Matrix::colMeans(norm[intersect(g, rownames(norm)), , drop = FALSE])
meta <- data.frame(cell = rownames(pd), sample = pd$sample, condition = pd$condition, program = pd$major_fibroblast_subtype, pseudotime = pd$Pseudotime_graph, D1 = score(d1)[rownames(pd)], D2 = score(d2)[rownames(pd)], predicted_doublet = !keep)

summarize_version <- function(df, version) df %>% group_by(sample, condition) %>% summarise(version = version, n_cells = n(), median_pseudotime = median(pseudotime), median_D1 = median(D1), median_D2 = median(D2), homeostatic_fraction = mean(program == "Homeostatic/regulatory fibroblasts"), matrix_fraction = mean(program == "Matrix-remodeling fibroblasts"), ecm_fraction = mean(program == "ECM-producing fibroblasts"), .groups = "drop")
donor <- bind_rows(summarize_version(meta, "all_current_cells"), summarize_version(meta[keep, ], "exclude_scDblFinder_predicted"))

exact_perm <- function(x, label) {
  obs <- mean(x[label == "Keloid"]) - mean(x[label == "Normal"])
  idx <- combn(seq_along(x), sum(label == "Keloid"), simplify = FALSE)
  null <- vapply(idx, function(k) { z <- rep("Normal", length(x)); z[k] <- "Keloid"; mean(x[z == "Keloid"]) - mean(x[z == "Normal"]) }, numeric(1))
  c(effect = obs, p = mean(abs(null) >= abs(obs) - 1e-12))
}
tests <- bind_rows(lapply(split(donor, donor$version), function(df) bind_rows(lapply(c("median_pseudotime", "median_D1", "median_D2", "homeostatic_fraction", "matrix_fraction", "ecm_fraction"), function(v) {
  z <- exact_perm(df[[v]], df$condition); data.frame(version = df$version[1], metric = v, keloid_minus_normal = z[1], exact_permutation_p = z[2])
}))))

dynamic <- read.csv(file.path(audit_root, "clean_programs", "clean_D1D2_gene_level_audit.csv"), stringsAsFactors = FALSE)
genes <- intersect(dynamic$gene, rownames(norm))
direction <- bind_rows(lapply(list(all_current_cells = rep(TRUE, nrow(pd)), exclude_scDblFinder_predicted = keep), function(k) {
  rho <- apply(as.matrix(norm[genes, k, drop = FALSE]), 1, function(x) cor(x, pd$Pseudotime_graph[k], method = "spearman"))
  data.frame(gene = genes, rho = rho, direction_concordant = sign(rho) == sign(dynamic$spearman_rho[match(genes, dynamic$gene)]))
}), .id = "version")
direction_summary <- direction %>% left_join(dynamic[, c("gene", "dynamic_module")], by = "gene") %>% group_by(version, dynamic_module) %>% summarise(gene_n = n(), direction_concordance = mean(direction_concordant), median_rho = median(rho), .groups = "drop")

write.csv(donor, file.path(out_dir, "doublet_sensitivity_donor_level_summary.csv"), row.names = FALSE)
write.csv(tests, file.path(out_dir, "doublet_sensitivity_donor_level_tests.csv"), row.names = FALSE)
write.csv(direction_summary, file.path(out_dir, "doublet_sensitivity_dynamic_direction_summary.csv"), row.names = FALSE)
writeLines(capture.output({ print(tests); print(direction_summary); print(sessionInfo()) }), file.path(out_dir, "doublet_sensitivity_report.txt"))
cat("doublet sensitivity complete\n"); print(tests); print(direction_summary)
