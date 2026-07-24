source(file.path(if (length(grep("^--file=", commandArgs(FALSE), value = TRUE))) dirname(chartr("\\", "/", sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))) else getwd(), "package_config.R"))
suppressPackageStartupMessages({
  library(Biobase)
  library(Matrix)
  library(monocle)
  library(mgcv)
  library(tradeSeq)
})

options(stringsAsFactors = FALSE)
set.seed(20260716)

out_dir <- file.path(audit_root, "donor_aware_trajectory")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

cds <- readRDS(file.path(input_root, "discovery_trajectory_graph_pseudotime_cds.rds"))
pd <- pData(cds)
audit <- read.csv(file.path(audit_root, "clean_programs", "clean_D1D2_gene_level_audit.csv"), stringsAsFactors = FALSE)
raw_all <- exprs(cds)

required_pd <- c("sample", "condition", "Size_Factor", "Pseudotime_graph")
if (!all(required_pd %in% colnames(pd))) stop("Missing phenotype columns: ", paste(setdiff(required_pd, colnames(pd)), collapse = ", "))
if (!all(audit$gene %in% rownames(raw_all))) stop("Some frozen dynamic genes are absent from the trajectory count matrix")

genes <- audit$gene
raw <- as.matrix(raw_all[genes, rownames(pd), drop = FALSE])
sf <- as.numeric(pd$Size_Factor)
if (any(!is.finite(sf)) || any(sf <= 0)) stop("Invalid Size_Factor values")
lognorm <- log1p(sweep(raw, 2, sf, FUN = "/"))
pt <- as.numeric(pd$Pseudotime_graph)
pt_scaled <- (pt - min(pt)) / diff(range(pt))
donor <- factor(pd$sample)

gene_info <- data.frame(
  gene = genes,
  dynamic_module = audit$dynamic_module,
  clean_module = ifelse(is.na(audit$clean_module), "none", audit$clean_module),
  expected_sign = ifelse(audit$dynamic_module == "D1_early_decreasing", -1, 1),
  stringsAsFactors = FALSE
)

safe_rho <- function(y, x) {
  if (length(y) < 3 || stats::sd(y) == 0 || stats::sd(x) == 0) return(NA_real_)
  suppressWarnings(stats::cor(y, x, method = "spearman", use = "complete.obs"))
}

sets <- list(
  dynamic_D1 = gene_info$gene[gene_info$dynamic_module == "D1_early_decreasing"],
  dynamic_D2 = gene_info$gene[gene_info$dynamic_module == "D2_late_increasing"],
  clean_D1 = gene_info$gene[gene_info$clean_module == "clean_D1"],
  clean_D2 = gene_info$gene[gene_info$clean_module == "clean_D2"]
)

summarize_support <- function(df, support_col, analysis) {
  do.call(rbind, lapply(names(sets), function(s) {
    z <- df[df$gene %in% sets[[s]], , drop = FALSE]
    evaluable <- !is.na(z[[support_col]])
    data.frame(
      analysis = analysis, gene_set = s, total_n = length(sets[[s]]),
      evaluable_n = sum(evaluable), supported_n = sum(z[[support_col]] %in% TRUE, na.rm = TRUE),
      supported_fraction_total = mean(z[[support_col]] %in% TRUE),
      supported_fraction_evaluable = if (sum(evaluable)) mean(z[[support_col]][evaluable] %in% TRUE) else NA_real_,
      stringsAsFactors = FALSE
    )
  }))
}

donor_counts <- as.data.frame(table(donor), stringsAsFactors = FALSE)
colnames(donor_counts) <- c("donor", "n_cells")
donor_counts$condition <- vapply(donor_counts$donor, function(d) unique(as.character(pd$condition[donor == d]))[1], character(1))
write.csv(donor_counts, file.path(out_dir, "donor_cell_counts.csv"), row.names = FALSE)

# Gene-level within-donor direction audit. The 11-cell donor remains in module-level
# summaries but is not used for gene-level direction calls.
per_donor <- do.call(rbind, lapply(levels(donor), function(d) {
  idx <- which(donor == d)
  min_detect <- max(5L, ceiling(0.10 * length(idx)))
  detected_cells <- rowSums(raw[, idx, drop = FALSE] > 0)
  rho <- if (length(unique(pt_scaled[idx])) >= 10) apply(lognorm[, idx, drop = FALSE], 1, safe_rho, x = pt_scaled[idx]) else rep(NA_real_, length(genes))
  evaluable <- length(idx) >= 20 & detected_cells >= min_detect & apply(lognorm[, idx, drop = FALSE], 1, sd) > 0 & is.finite(rho)
  rho[!evaluable] <- NA_real_
  data.frame(
    gene = genes, donor = d, n_cells = length(idx), detected_cells = detected_cells,
    rho = rho, evaluable = evaluable, expected_sign = gene_info$expected_sign,
    concordant = ifelse(evaluable, sign(rho) == gene_info$expected_sign, NA), stringsAsFactors = FALSE
  )
}))
write.csv(per_donor, file.path(out_dir, "gene_by_donor_spearman.csv"), row.names = FALSE)

per_donor_summary <- do.call(rbind, lapply(genes, function(g) {
  z <- per_donor[per_donor$gene == g & per_donor$evaluable, , drop = FALSE]
  data.frame(
    gene = g, n_evaluable_donors = nrow(z), n_concordant_donors = sum(z$concordant, na.rm = TRUE),
    concordant_fraction = if (nrow(z)) mean(z$concordant) else NA_real_,
    median_within_donor_rho = if (nrow(z)) median(z$rho) else NA_real_,
    majority_concordant = if (nrow(z) >= 3) mean(z$concordant) > 0.5 else NA,
    two_thirds_concordant = if (nrow(z) >= 3) mean(z$concordant) >= 2/3 else NA,
    stringsAsFactors = FALSE
  )
}))
per_donor_summary <- merge(gene_info, per_donor_summary, by = "gene", all.x = TRUE, sort = FALSE)
write.csv(per_donor_summary, file.path(out_dir, "per_donor_direction_summary.csv"), row.names = FALSE)

# Fixed-trajectory leave-one-donor-out association audit.
lodo <- do.call(rbind, lapply(levels(donor), function(d) {
  idx <- which(donor != d)
  rho <- apply(lognorm[, idx, drop = FALSE], 1, safe_rho, x = pt_scaled[idx])
  n <- length(idx)
  t_value <- rho * sqrt((n - 2) / pmax(1e-12, 1 - rho^2))
  p_value <- 2 * pt(-abs(t_value), df = n - 2)
  fdr <- p.adjust(p_value, "BH")
  data.frame(
    gene = genes, left_out_donor = d, n_cells = n, rho = rho, p_value = p_value, fdr = fdr,
    expected_sign = gene_info$expected_sign,
    retained = sign(rho) == gene_info$expected_sign & abs(rho) >= 0.20 & fdr < 0.05,
    stringsAsFactors = FALSE
  )
}))
write.csv(lodo, file.path(out_dir, "fixed_trajectory_LODO_gene_results.csv"), row.names = FALSE)
lodo_summary <- do.call(rbind, lapply(genes, function(g) {
  z <- lodo[lodo$gene == g, , drop = FALSE]
  data.frame(
    gene = g, retained_n = sum(z$retained), retained_fraction = mean(z$retained),
    retained_all_6 = all(z$retained), retained_at_least_5 = sum(z$retained) >= 5,
    retained_majority = sum(z$retained) >= 4, stringsAsFactors = FALSE
  )
}))
lodo_summary <- merge(gene_info, lodo_summary, by = "gene", all.x = TRUE, sort = FALSE)
write.csv(lodo_summary, file.path(out_dir, "fixed_trajectory_LODO_summary.csv"), row.names = FALSE)

# Equal-donor downsampling evaluates direction without using cell-level p-values.
run_downsampling <- function(use_donors, n_each, n_rep, label) {
  rho_mat <- matrix(NA_real_, nrow = length(genes), ncol = n_rep)
  for (b in seq_len(n_rep)) {
    idx <- unlist(lapply(use_donors, function(d) sample(which(donor == d), n_each, replace = FALSE)), use.names = FALSE)
    rho_mat[, b] <- apply(lognorm[, idx, drop = FALSE], 1, safe_rho, x = pt_scaled[idx])
  }
  data.frame(
    gene = genes, scheme = label, n_repeats = n_rep,
    median_rho = apply(rho_mat, 1, median, na.rm = TRUE),
    expected_direction_fraction = rowMeans(sign(rho_mat) == gene_info$expected_sign, na.rm = TRUE),
    strong_expected_fraction = rowMeans(sign(rho_mat) == gene_info$expected_sign & abs(rho_mat) >= 0.20, na.rm = TRUE),
    direction_retained_in_at_least_80pct_repeats =
      rowMeans(sign(rho_mat) == gene_info$expected_sign, na.rm = TRUE) >= 0.80,
    stringsAsFactors = FALSE
  )
}
downsampling <- rbind(
  run_downsampling(levels(donor), min(table(donor)), 100, "all_6_donors_n11"),
  run_downsampling(names(table(donor))[table(donor) >= 43], 43, 100, "five_donors_n43")
)
downsampling <- merge(gene_info, downsampling, by = "gene", all.y = TRUE, sort = FALSE)
write.csv(downsampling, file.path(out_dir, "equal_donor_downsampling_results.csv"), row.names = FALSE)

fit_gam_set <- function(model_type) {
  message("Fitting ", model_type, " GAMs")
  donor_weight <- 1 / as.numeric(table(donor)[donor])
  donor_weight <- donor_weight / mean(donor_weight)
  ans <- lapply(seq_along(genes), function(i) {
    dat <- data.frame(y = lognorm[i, ], pt_scaled = pt_scaled, donor = donor, donor_weight = donor_weight)
    fit <- tryCatch({
      if (model_type == "donor_fixed") mgcv::gam(y ~ donor + s(pt_scaled, k = 6), data = dat, method = "REML")
      else if (model_type == "donor_random") mgcv::gam(y ~ s(pt_scaled, k = 6) + s(donor, bs = "re"), data = dat, method = "REML")
      else mgcv::gam(y ~ s(pt_scaled, k = 6), data = dat, weights = donor_weight, method = "REML")
    }, error = function(e) NULL)
    if (is.null(fit)) return(data.frame(gene = genes[i], smooth_p = NA_real_, edf = NA_real_, endpoint_delta = NA_real_, fit_ok = FALSE))
    st <- summary(fit)$s.table
    rr <- grep("^s\\(pt_scaled\\)", rownames(st))[1]
    nd <- data.frame(pt_scaled = c(0, 1), donor = factor(rep(levels(donor)[1], 2), levels = levels(donor)), donor_weight = 1)
    pred <- tryCatch(as.numeric(predict(fit, newdata = nd)), error = function(e) c(NA_real_, NA_real_))
    data.frame(gene = genes[i], smooth_p = st[rr, "p-value"], edf = st[rr, "edf"], endpoint_delta = pred[2] - pred[1], fit_ok = TRUE)
  })
  z <- do.call(rbind, ans)
  z$fdr <- p.adjust(z$smooth_p, "BH")
  z <- merge(gene_info, z, by = "gene", all.x = TRUE, sort = FALSE)
  z$direction_concordant <- sign(z$endpoint_delta) == z$expected_sign
  z$supported <- z$fit_ok & z$fdr < 0.05 & z$direction_concordant
  z$model <- model_type
  z
}

gam_all <- rbind(fit_gam_set("donor_fixed"), fit_gam_set("donor_random"), fit_gam_set("donor_equal_weight"))
write.csv(gam_all, file.path(out_dir, "donor_aware_GAM_gene_results.csv"), row.names = FALSE)

message("Fitting single-lineage tradeSeq with donor covariates")
# tradeSeq uses y ~ -1 + U + smooth + offset, so the design must include an intercept.
donor_design <- model.matrix(~ donor)
trade_fit <- tradeSeq::fitGAM(
  counts = round(raw), pseudotime = matrix(pt_scaled, ncol = 1),
  cellWeights = matrix(1, nrow = length(pt_scaled), ncol = 1),
  U = donor_design, nknots = 6, verbose = TRUE, parallel = FALSE
)
trade_res <- as.data.frame(tradeSeq::associationTest(trade_fit))
trade_res$gene <- rownames(trade_res)
trade_res$tradeSeq_fdr <- p.adjust(trade_res$pvalue, "BH")
trade_res <- merge(gene_info, trade_res, by = "gene", all.x = TRUE, sort = FALSE)
trade_res$supported <- is.finite(trade_res$tradeSeq_fdr) & trade_res$tradeSeq_fdr < 0.05
write.csv(trade_res, file.path(out_dir, "tradeSeq_with_donor_covariates.csv"), row.names = FALSE)

# Module-level within-donor trends use all six donors, including the 11-cell donor.
score_mean <- function(g) colMeans(lognorm[intersect(g, rownames(lognorm)), , drop = FALSE])
module_dat <- data.frame(
  donor = as.character(donor), condition = as.character(pd$condition), pseudotime = pt,
  clean_D1 = score_mean(sets$clean_D1), clean_D2 = score_mean(sets$clean_D2), stringsAsFactors = FALSE
)
module_by_donor <- do.call(rbind, lapply(split(module_dat, module_dat$donor), function(z) data.frame(
  donor = z$donor[1], condition = z$condition[1], n_cells = nrow(z),
  pt_min = min(z$pseudotime), pt_median = median(z$pseudotime), pt_max = max(z$pseudotime),
  clean_D1_rho = safe_rho(z$clean_D1, z$pseudotime), clean_D2_rho = safe_rho(z$clean_D2, z$pseudotime),
  clean_D1_expected = safe_rho(z$clean_D1, z$pseudotime) < 0,
  clean_D2_expected = safe_rho(z$clean_D2, z$pseudotime) > 0, stringsAsFactors = FALSE
)))
write.csv(module_by_donor, file.path(out_dir, "module_score_within_donor_trends.csv"), row.names = FALSE)
module_lodo <- do.call(rbind, lapply(levels(donor), function(d) {
  z <- module_dat[module_dat$donor != d, ]
  data.frame(left_out_donor = d, n_cells = nrow(z), clean_D1_rho = safe_rho(z$clean_D1, z$pseudotime), clean_D2_rho = safe_rho(z$clean_D2, z$pseudotime))
}))
write.csv(module_lodo, file.path(out_dir, "module_score_fixed_trajectory_LODO.csv"), row.names = FALSE)

summary_rows <- list(
  summarize_support(per_donor_summary, "majority_concordant", "Per-donor majority direction among genes evaluable in >=3 donors"),
  summarize_support(per_donor_summary, "two_thirds_concordant", "Per-donor >=2/3 direction among genes evaluable in >=3 donors"),
  summarize_support(lodo_summary, "retained_at_least_5", "Fixed-trajectory LODO retained in >=5/6"),
  summarize_support(
    downsampling[downsampling$scheme == "all_6_donors_n11", ],
    "direction_retained_in_at_least_80pct_repeats",
    "Equal-donor downsampling: 6 donors x 11 cells; expected direction in >=80/100 repeats"
  ),
  summarize_support(
    downsampling[downsampling$scheme == "five_donors_n43", ],
    "direction_retained_in_at_least_80pct_repeats",
    "Equal-donor downsampling: 5 donors x 43 cells; expected direction in >=80/100 repeats"
  )
)
for (m in unique(gam_all$model)) summary_rows[[length(summary_rows) + 1]] <- summarize_support(gam_all[gam_all$model == m, ], "supported", paste0("GAM: ", m, ", FDR<0.05 and expected endpoint direction"))
summary_rows[[length(summary_rows) + 1]] <- summarize_support(trade_res, "supported", "tradeSeq with donor covariates, association FDR<0.05")
summary_table <- do.call(rbind, summary_rows)
write.csv(summary_table, file.path(out_dir, "gene_set_robustness_summary.csv"), row.names = FALSE)

parameters <- data.frame(
  parameter = c(
    "seed", "response", "pseudotime", "gene_level_evaluable_donor",
    "per_donor_direction_support", "LODO_definition", "downsampling_schemes", "downsampling_summary_rule",
    "donor_fixed_GAM", "donor_random_GAM", "donor_equal_weight_GAM",
    "GAM_support", "tradeSeq_input", "tradeSeq_donor_design",
    "tradeSeq_library_size_adjustment", "tradeSeq_nknots", "tradeSeq_nknots_basis"
  ),
  value = c(
    "20260716", "log1p(raw_count / Monocle2 Size_Factor)", "Pseudotime_graph scaled to [0,1]",
    "donor n>=20; detected cells>=max(5, ceiling(0.10*n)); expression variance>0; >=10 unique pseudotime values",
    "expected sign in >50% or >=2/3 of at least three evaluable donors",
    "fixed graph/pseudotime; expected sign; |rho|>=0.20; BH FDR<0.05; retained in >=5/6 deletions",
    "100 repeats: all six donors x 11 cells; five donors with n>=43 x 43 cells",
    "gene retains expected direction in at least 80 of 100 repeats",
    "response ~ donor + s(scaled_pseudotime, k=6)",
    "response ~ s(scaled_pseudotime, k=6) + s(donor, bs='re')",
    "response ~ s(scaled_pseudotime, k=6), cell weight proportional to 1/n_cells_in_donor",
    "smooth-term BH FDR<0.05 and expected endpoint direction",
    "integer raw counts", "model.matrix(~ donor), including intercept and five donor indicators",
    "default TMM-derived effective-library-size offset", "6",
    "prespecified tradeSeq default; two-lineage evaluateK sensitivity compared 3-8 knots in 500 fixed-seed sampled genes"
  ), stringsAsFactors = FALSE
)
write.csv(parameters, file.path(out_dir, "donor_aware_trajectory_parameters.csv"), row.names = FALSE)

writeLines(capture.output({
  cat("DONOR-AWARE TRAJECTORY ROBUSTNESS\n")
  print(donor_counts)
  cat("\nMODULE-LEVEL WITHIN-DONOR TRENDS\n"); print(module_by_donor)
  cat("\nGENE-SET SUMMARY\n"); print(summary_table)
  cat("\nSESSION INFO\n"); print(sessionInfo())
}), file.path(out_dir, "donor_aware_trajectory_report.txt"))

cat("donor-aware trajectory robustness complete\n")
