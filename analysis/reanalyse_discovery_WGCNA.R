source(file.path(if (length(grep("^--file=", commandArgs(FALSE), value = TRUE))) dirname(chartr("\\", "/", sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))) else getwd(), "package_config.R"))
suppressPackageStartupMessages({
  library(monocle)
  library(Biobase)
  library(Matrix)
  library(WGCNA)
  library(dplyr)
})

options(stringsAsFactors = FALSE)
allowWGCNAThreads(nThreads = 4)
set.seed(20260717)

out_name <- Sys.getenv("WGCNA_OUT_DIR", unset = "WGCNA")
out_dir <- file.path(audit_root, out_name)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

cds <- readRDS(file.path(input_root, "discovery_trajectory_graph_pseudotime_cds.rds"))
pd <- pData(cds)
pd$sample <- as.character(pd$sample)
pd$condition <- as.character(pd$condition)
pd$major_fibroblast_subtype <- as.character(pd$major_fibroblast_subtype)
disc_counts <- exprs(cds)
disc_sf <- as.numeric(pd$Size_Factor)
names(disc_sf) <- rownames(pd)
disc_data <- log1p(t(t(disc_counts) / disc_sf[colnames(disc_counts)]))

clean <- read.csv(file.path(audit_root, "clean_programs", "clean_D1D2_gene_level_audit.csv"), stringsAsFactors = FALSE)
d1 <- unique(clean$gene[clean$clean_module == "clean_D1" & !is.na(clean$clean_module)])
d2 <- unique(clean$gene[clean$clean_module == "clean_D2" & !is.na(clean$clean_module)])
score_mean <- function(mat, genes) {
  g <- intersect(genes, rownames(mat))
  Matrix::colMeans(mat[g, , drop = FALSE])
}

meta <- data.frame(
  cell = rownames(pd), sample = pd$sample, condition = pd$condition,
  state = pd$major_fibroblast_subtype, pseudotime = as.numeric(pd$Pseudotime_graph),
  clean_D1 = score_mean(disc_data, d1)[rownames(pd)],
  clean_D2 = score_mean(disc_data, d2)[rownames(pd)],
  stringsAsFactors = FALSE
) %>%
  group_by(sample, state) %>%
  arrange(pseudotime, .by_group = TRUE) %>%
  mutate(bin = ceiling(row_number() / 25)) %>%
  ungroup() %>%
  mutate(metacell = paste(sample, state, bin, sep = "__"))

mc_candidates <- meta %>%
  group_by(metacell, sample, condition, state, bin) %>%
  summarise(
    n_cells = n(), mean_pseudotime = mean(pseudotime),
    clean_D1 = mean(clean_D1), clean_D2 = mean(clean_D2), .groups = "drop"
  )
mc_meta <- mc_candidates %>% filter(n_cells >= 8)

all_discovery_donors <- sort(unique(meta$sample))
all_discovery_states <- sort(unique(meta$state))
donor_state_grid <- expand.grid(
  sample = all_discovery_donors,
  state = all_discovery_states,
  stringsAsFactors = FALSE
)
donor_state_eligibility <- donor_state_grid %>%
  left_join(meta %>% count(sample, state, name = "input_cell_n"), by = c("sample", "state")) %>%
  left_join(
    mc_candidates %>% group_by(sample, state) %>% summarise(candidate_bin_n = n(), .groups = "drop"),
    by = c("sample", "state")
  ) %>%
  left_join(
    mc_meta %>% group_by(sample, state) %>% summarise(
      retained_metacell_n = n(), retained_cell_n = sum(n_cells), .groups = "drop"
    ),
    by = c("sample", "state")
  ) %>%
  mutate(
    across(c(input_cell_n, candidate_bin_n, retained_metacell_n, retained_cell_n), ~ tidyr::replace_na(.x, 0L)),
    eligibility_reason = ifelse(
      retained_metacell_n > 0,
      "at least one within-donor state-specific pseudotime bin contained >=8 cells",
      "all within-donor state-specific pseudotime bins contained <8 cells"
    )
  )

# The discovery network gene space is defined without loading or querying any
# external cohort. Validation detectability is considered only after modules
# have been frozen.
disc_detect <- Matrix::rowMeans(disc_counts > 0)
eligible_genes <- names(disc_detect)[is.finite(disc_detect) & disc_detect >= 0.05]
mc_expr_all <- sapply(mc_meta$metacell, function(mc) {
  cells <- meta$cell[meta$metacell == mc]
  Matrix::rowMeans(disc_data[eligible_genes, cells, drop = FALSE])
})
gene_var <- matrixStats::rowVars(as.matrix(mc_expr_all))
names(gene_var) <- rownames(mc_expr_all)
ranked_genes <- names(sort(gene_var, decreasing = TRUE, na.last = NA))
gene_pool <- head(ranked_genes, min(3000L, length(ranked_genes)))

datExpr <- t(as.matrix(mc_expr_all[gene_pool, , drop = FALSE]))
rownames(datExpr) <- mc_meta$metacell
good <- goodSamplesGenes(datExpr, verbose = 0)
datExpr <- datExpr[good$goodSamples, good$goodGenes, drop = FALSE]
mc_meta <- mc_meta[match(rownames(datExpr), mc_meta$metacell), , drop = FALSE]

powers <- c(1:10, 12, 14, 16, 18, 20)
sft <- pickSoftThreshold(datExpr, powerVector = powers, networkType = "signed", verbose = 0)
fit <- sft$fitIndices
ok <- which(fit$SFT.R.sq >= 0.80)
soft_power <- if (length(ok)) fit$Power[min(ok)] else fit$Power[which.max(fit$SFT.R.sq)]

net <- blockwiseModules(
  datExpr, power = soft_power, networkType = "signed", TOMType = "signed",
  minModuleSize = 25, mergeCutHeight = 0.25, numericLabels = TRUE,
  pamRespectsDendro = FALSE, maxBlockSize = ncol(datExpr), verbose = 2
)
colors <- labels2colors(net$colors)
names(colors) <- colnames(datExpr)
MEs <- orderMEs(moduleEigengenes(datExpr, colors = colors)$eigengenes)
traits <- mc_meta[, c("mean_pseudotime", "clean_D1", "clean_D2")]

weighted_rank_cor <- function(x, y, w) {
  keep <- is.finite(x) & is.finite(y) & is.finite(w) & w > 0
  x <- x[keep]; y <- y[keep]; w <- w[keep]
  if (length(x) < 4L || sd(x) == 0 || sd(y) == 0) return(NA_real_)
  xr <- rank(x, ties.method = "average")
  yr <- rank(y, ties.method = "average")
  mx <- sum(w * xr) / sum(w)
  my <- sum(w * yr) / sum(w)
  denom <- sqrt(sum(w * (xr - mx)^2) * sum(w * (yr - my)^2))
  if (!is.finite(denom) || denom == 0) return(NA_real_)
  sum(w * (xr - mx) * (yr - my)) / denom
}

donor_equal_cor <- function(x, y, donor) {
  donor <- as.character(donor)
  tab <- table(donor)
  w <- 1 / as.numeric(tab[donor])
  weighted_rank_cor(x, y, w)
}

discovery_donors <- sort(unique(mc_meta$sample))
module_names <- setdiff(sub("^ME", "", colnames(MEs)), "grey")
module_trait <- do.call(rbind, lapply(module_names, function(module) {
  x <- MEs[[paste0("ME", module)]]
  do.call(rbind, lapply(colnames(traits), function(trait) {
    y <- traits[[trait]]
    lodo <- vapply(discovery_donors, function(drop) {
      keep <- mc_meta$sample != drop
      donor_equal_cor(x[keep], y[keep], mc_meta$sample[keep])
    }, numeric(1))
    data.frame(
      module = module, trait = trait,
      full_metacell_rho_descriptive = suppressWarnings(cor(x, y, method = "spearman")),
      donor_equal_rho = donor_equal_cor(x, y, mc_meta$sample),
      donor_equal_lodo_median_rho = median(lodo, na.rm = TRUE),
      donor_equal_lodo_min_abs_rho = min(abs(lodo), na.rm = TRUE),
      donor_equal_lodo_sign_consistent = length(unique(sign(lodo[is.finite(lodo)]))) == 1L,
      stringsAsFactors = FALSE
    )
  }))
}))

select_module <- function(trait) {
  z <- module_trait[module_trait$trait == trait & is.finite(module_trait$donor_equal_rho), , drop = FALSE]
  z$module[which.max(z$donor_equal_rho)]
}
d1_module <- select_module("clean_D1")
d2_module <- select_module("clean_D2")
selected_pairs <- data.frame(
  program = c("clean_D1", "clean_D2"),
  aligned_module = c(d1_module, d2_module),
  selection_metric = "maximum positive donor-equal Spearman correlation",
  stringsAsFactors = FALSE
)

cluster_bootstrap <- function(x, y, donor, n_boot = 2000L, seed = 20260717L) {
  donors <- sort(unique(as.character(donor)))
  set.seed(seed)
  replicate(n_boot, {
    draws <- sample(donors, length(donors), replace = TRUE)
    pieces <- lapply(seq_along(draws), function(i) {
      idx <- which(donor == draws[[i]])
      data.frame(x = x[idx], y = y[idx], draw = paste0("draw", i))
    })
    boot <- do.call(rbind, pieces)
    donor_equal_cor(boot$x, boot$y, boot$draw)
  })
}

selected_sensitivity <- do.call(rbind, lapply(seq_len(nrow(selected_pairs)), function(i) {
  module <- selected_pairs$aligned_module[[i]]
  trait <- selected_pairs$program[[i]]
  x <- MEs[[paste0("ME", module)]]
  y <- traits[[trait]]
  if (donor_equal_cor(x, y, mc_meta$sample) < 0) x <- -x
  lodo <- vapply(discovery_donors, function(drop) {
    keep <- mc_meta$sample != drop
    donor_equal_cor(x[keep], y[keep], mc_meta$sample[keep])
  }, numeric(1))
  boot <- cluster_bootstrap(x, y, mc_meta$sample, n_boot = 2000L, seed = 20260717L + i)
  data.frame(
    program = trait, aligned_module = module,
    metacells_n = length(x), represented_donors_n = length(discovery_donors),
    full_metacell_rho_descriptive = cor(x, y, method = "spearman"),
    donor_equal_rho = donor_equal_cor(x, y, mc_meta$sample),
    donor_equal_lodo_min_rho = min(lodo, na.rm = TRUE),
    donor_equal_lodo_median_rho = median(lodo, na.rm = TRUE),
    donor_equal_lodo_positive_fraction = mean(lodo > 0, na.rm = TRUE),
    donor_equal_cluster_bootstrap_median = median(boot, na.rm = TRUE),
    donor_equal_cluster_bootstrap_q025 = unname(quantile(boot, 0.025, na.rm = TRUE)),
    donor_equal_cluster_bootstrap_q975 = unname(quantile(boot, 0.975, na.rm = TRUE)),
    donor_equal_cluster_bootstrap_positive_fraction = mean(boot > 0, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}))

fisher_overlap <- function(genes, module) {
  universe <- names(colors)
  in_set <- universe %in% genes
  in_module <- colors == module
  ft <- fisher.test(table(factor(in_set, levels = c(FALSE, TRUE)), factor(in_module, levels = c(FALSE, TRUE))), alternative = "greater")
  data.frame(
    module = module,
    overlap_n = sum(in_set & in_module),
    set_n_in_network = sum(in_set),
    set_total_n = length(unique(genes)),
    module_n = sum(in_module),
    network_gene_n = length(universe),
    odds_ratio_exploratory = unname(ft$estimate),
    p_value_exploratory = ft$p.value,
    inference_role = "descriptive post-selection overlap",
    stringsAsFactors = FALSE
  )
}
overlap <- bind_rows(
  cbind(program = "clean_D1", fisher_overlap(d1, d1_module)),
  cbind(program = "clean_D2", fisher_overlap(d2, d2_module))
)
overlap$fdr_exploratory <- p.adjust(overlap$p_value_exploratory, "BH")

disc_group <- paste(pd$sample, pd$major_fibroblast_subtype, sep = "__")
disc_pb <- t(sapply(split(rownames(pd), disc_group), function(cells) {
  Matrix::rowMeans(disc_data[names(colors), cells, drop = FALSE])
}))

gene_selection <- data.frame(
  gene = names(disc_detect),
  discovery_detection_fraction = as.numeric(disc_detect),
  discovery_detection_eligible = names(disc_detect) %in% eligible_genes,
  metacell_variance = as.numeric(gene_var[names(disc_detect)]),
  selected_top_variable = names(disc_detect) %in% colnames(datExpr),
  clean_program = ifelse(names(disc_detect) %in% d1, "clean_D1", ifelse(names(disc_detect) %in% d2, "clean_D2", NA)),
  stringsAsFactors = FALSE
)

write.csv(rename(mc_meta, program = state), file.path(out_dir, "discovery_metacell_composition.csv"), row.names = FALSE)
write.csv(donor_state_eligibility, file.path(out_dir, "WGCNA_metacell_donor_state_eligibility.csv"), row.names = FALSE)
write.csv(gene_selection, file.path(out_dir, "WGCNA_discovery_only_gene_selection.csv"), row.names = FALSE)
write.csv(data.frame(Power = fit$Power, SFT_R2 = fit$SFT.R.sq, mean_connectivity = fit$mean.k.), file.path(out_dir, "WGCNA_soft_threshold.csv"), row.names = FALSE)
write.csv(data.frame(gene = names(colors), module = unname(colors)), file.path(out_dir, "WGCNA_gene_module_assignment.csv"), row.names = FALSE)
write.csv(module_trait, file.path(out_dir, "WGCNA_module_trait_correlations.csv"), row.names = FALSE)
write.csv(selected_pairs, file.path(out_dir, "WGCNA_selected_aligned_modules.csv"), row.names = FALSE)
write.csv(data.frame(
  donor = all_discovery_donors,
  trajectory_cell_n = as.integer(table(factor(meta$sample, levels = all_discovery_donors))),
  metacells_n = as.integer(table(factor(mc_meta$sample, levels = all_discovery_donors))),
  represented_in_WGCNA = all_discovery_donors %in% discovery_donors
), file.path(out_dir, "WGCNA_metacell_donor_counts.csv"), row.names = FALSE)
write.csv(selected_sensitivity, file.path(out_dir, "WGCNA_selected_module_donor_sensitivity.csv"), row.names = FALSE)
write.csv(overlap, file.path(out_dir, "WGCNA_clean_D1D2_overlap.csv"), row.names = FALSE)
write.csv(data.frame(unit = rownames(disc_pb), cohort = "Discovery", disc_pb, check.names = FALSE), file.path(out_dir, "matched_discovery_donor_state_pseudobulk.csv"), row.names = FALSE)
write.csv(data.frame(
  parameter = c(
    "seed", "gene_space", "validation_used_in_gene_selection", "clean_genes_forced_into_network",
    "discovery_detection_threshold", "variance_gene_count", "metacell_definition", "metacell_min_cells",
    "soft_power", "network_gene_n", "discovery_metacells", "trajectory_donors", "represented_donors",
    "unrepresented_donor", "unrepresented_donor_reason",
    "module_trait_primary_metric", "module_trait_LODO", "donor_cluster_bootstrap",
    "D1_aligned_module", "D2_aligned_module"
  ),
  value = c(
    20260717, "discovery cohort only", FALSE, FALSE, 0.05, 3000,
    "sample x state x consecutive 25-cell pseudotime bin", 8, soft_power, ncol(datExpr), nrow(datExpr),
    length(all_discovery_donors), length(discovery_donors),
    paste(setdiff(all_discovery_donors, discovery_donors), collapse = ";"),
    "11 trajectory cells distributed across three states (1/6/4); no within-state bin reached the >=8-cell metacell threshold",
    "donor-equal weighted rank correlation", "donor-equal leave-one-represented-donor-out",
    "2000 donor-cluster replicates with equal total weight per sampled donor draw", d1_module, d2_module
  )
), file.path(out_dir, "WGCNA_parameters.csv"), row.names = FALSE)

saveRDS(list(datExpr = datExpr, colors = colors, MEs = MEs, mc_meta = mc_meta, net = net, selected_pairs = selected_pairs), file.path(out_dir, "WGCNA_discovery_only_result.rds"), compress = FALSE)
writeLines(capture.output({
  cat("Discovery-only WGCNA\n")
  cat("D1-aligned module:", d1_module, "D2-aligned module:", d2_module, "\n")
  cat("Trajectory donors:", length(all_discovery_donors), "represented in WGCNA:", length(discovery_donors), "\n")
  cat("Unrepresented donor(s):", paste(setdiff(all_discovery_donors, discovery_donors), collapse = ", "), "\n")
  print(donor_state_eligibility[donor_state_eligibility$sample %in% setdiff(all_discovery_donors, discovery_donors), ])
  cat("Metacells per represented donor:\n"); print(table(mc_meta$sample))
  cat("Selected-module donor sensitivity:\n"); print(selected_sensitivity)
  cat("Clean-program overlap (post-selection descriptive):\n"); print(overlap)
  print(sessionInfo())
}), file.path(out_dir, "WGCNA_report.txt"))

cat("discovery-only WGCNA complete\n")
print(selected_sensitivity)
print(overlap)
