source(file.path(if (length(grep("^--file=", commandArgs(FALSE), value = TRUE))) dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]), winslash = "/")) else getwd(), "package_config.R"))
suppressPackageStartupMessages({
  library(pROC)
  library(dplyr)
  library(ggplot2)
})

set.seed(20260711)
out_dir <- file.path(audit_root, "bulk_recurrence")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

expr_path <- Sys.getenv("PRJNA813172_LOG2TPM")
meta_path <- Sys.getenv("PRJNA813172_METADATA")
expr <- read.csv(expr_path, row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)
meta <- read.csv(meta_path, check.names = FALSE, stringsAsFactors = FALSE)
meta <- meta[match(colnames(expr), meta$Run), ]
stopifnot(all(meta$Run == colnames(expr)))
y <- as.integer(meta$recurrence_binary_clean)
names(y) <- meta$Run

module_audit <- read.csv(
  file.path(audit_root, "clean_programs", "clean_D1D2_gene_level_audit.csv"),
  stringsAsFactors = FALSE
)

sets <- list(
  clean_D1 = unique(module_audit$gene[!is.na(module_audit$clean_module) & module_audit$clean_module == "clean_D1"]),
  clean_D2 = unique(module_audit$gene[!is.na(module_audit$clean_module) & module_audit$clean_module == "clean_D2"]),
  five_hub = c("CHI3L1", "IL1RN", "MMP7", "TNFAIP3", "TNFAIP6")
)

score_set <- function(genes) {
  detected <- intersect(genes, rownames(expr))
  x <- as.matrix(expr[detected, , drop = FALSE])
  z <- t(scale(t(x)))
  z[!is.finite(z)] <- 0
  colMeans(z)
}
scores <- as.data.frame(lapply(sets, score_set), check.names = FALSE)
rownames(scores) <- colnames(expr)
scores$clean_D1_loss <- -scores$clean_D1

auc_rank <- function(labels, score) {
  n1 <- sum(labels == 1); n0 <- sum(labels == 0)
  (sum(rank(score, ties.method = "average")[labels == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}
average_precision <- function(labels, score) {
  o <- order(score, decreasing = TRUE)
  yy <- labels[o]
  precision <- cumsum(yy) / seq_along(yy)
  sum(precision[yy == 1]) / sum(yy)
}
strat_boot_idx <- function(labels) {
  c(sample(which(labels == 1), sum(labels == 1), replace = TRUE),
    sample(which(labels == 0), sum(labels == 0), replace = TRUE))
}

score_names <- c("clean_D1_loss", "clean_D2", "five_hub")
n_boot <- 5000
n_perm <- 10000
stats <- lapply(score_names, function(nm) {
  s <- scores[[nm]]
  base_nm <- sub("_loss$", "", nm)
  gene_set <- if (base_nm %in% names(sets)) sets[[base_nm]] else character()
  auc <- auc_rank(y, s)
  ap <- average_precision(y, s)
  boot_auc <- boot_ap <- numeric(n_boot)
  for (b in seq_len(n_boot)) {
    idx <- strat_boot_idx(y)
    boot_auc[b] <- auc_rank(y[idx], s[idx])
    boot_ap[b] <- average_precision(y[idx], s[idx])
  }
  null_auc <- null_ap <- numeric(n_perm)
  for (b in seq_len(n_perm)) {
    yp <- sample(y)
    null_auc[b] <- auc_rank(yp, s)
    null_ap[b] <- average_precision(yp, s)
  }
  wt <- suppressWarnings(wilcox.test(s[y == 1], s[y == 0], exact = FALSE))
  data.frame(
    score = nm,
    gene_n_defined = length(gene_set),
    gene_n_detected = length(intersect(gene_set, rownames(expr))),
    recurrence_median = median(s[y == 1]),
    nonrecurrence_median = median(s[y == 0]),
    wilcoxon_p = wt$p.value,
    auc = auc,
    auc_ci_low = quantile(boot_auc, 0.025, na.rm = TRUE),
    auc_ci_high = quantile(boot_auc, 0.975, na.rm = TRUE),
    auc_permutation_p = (1 + sum(null_auc >= auc)) / (n_perm + 1),
    auc_permutation_two_sided_p = (1 + sum(abs(null_auc - 0.5) >= abs(auc - 0.5))) / (n_perm + 1),
    average_precision = ap,
    ap_ci_low = quantile(boot_ap, 0.025, na.rm = TRUE),
    ap_ci_high = quantile(boot_ap, 0.975, na.rm = TRUE),
    ap_permutation_p = (1 + sum(null_ap >= ap)) / (n_perm + 1),
    prevalence_baseline_ap = mean(y),
    stringsAsFactors = FALSE
  )
})
stats <- do.call(rbind, stats)

primary <- c("clean_D1_loss", "clean_D2", "five_hub")
stats$primary_translation_fdr <- NA_real_
stats$primary_translation_fdr[match(primary, stats$score)] <- p.adjust(stats$wilcoxon_p[match(primary, stats$score)], "BH")

comparisons <- list(
  clean_D1loss_vs_fivehub = c("clean_D1_loss", "five_hub"),
  clean_D1loss_vs_clean_D2 = c("clean_D1_loss", "clean_D2")
)
pairwise <- do.call(rbind, lapply(names(comparisons), function(nm) {
  pair <- comparisons[[nm]]
  delta <- numeric(n_boot)
  for (b in seq_len(n_boot)) {
    idx <- strat_boot_idx(y)
    delta[b] <- auc_rank(y[idx], scores[[pair[1]]][idx]) - auc_rank(y[idx], scores[[pair[2]]][idx])
  }
  obs <- auc_rank(y, scores[[pair[1]]]) - auc_rank(y, scores[[pair[2]]])
  data.frame(
    comparison = nm, score_a = pair[1], score_b = pair[2], observed_auc_difference = obs,
    bootstrap_ci_low = quantile(delta, 0.025, na.rm = TRUE),
    bootstrap_ci_high = quantile(delta, 0.975, na.rm = TRUE),
    paired_bootstrap_two_sided_p = 2 * min(mean(delta <= 0), mean(delta >= 0)),
    stringsAsFactors = FALSE
  )
}))

positive_idx <- which(y == 1)
leave_one_positive <- do.call(rbind, lapply(positive_idx, function(drop) {
  keep <- setdiff(seq_along(y), drop)
  do.call(rbind, lapply(primary, function(nm) data.frame(
    omitted_run = names(y)[drop], score = nm,
    auc = auc_rank(y[keep], scores[[nm]][keep]),
    average_precision = average_precision(y[keep], scores[[nm]][keep]),
    stringsAsFactors = FALSE
  )))
}))

coverage <- do.call(rbind, lapply(names(sets), function(nm) data.frame(
  gene_set = nm, defined_n = length(sets[[nm]]), detected_n = length(intersect(sets[[nm]], rownames(expr))),
  missing_genes = paste(setdiff(sets[[nm]], rownames(expr)), collapse = ";"), stringsAsFactors = FALSE
)))

score_out <- cbind(
  data.frame(Run = rownames(scores), recurrence = y, stringsAsFactors = FALSE),
  scores[, primary, drop = FALSE]
)
write.csv(score_out, file.path(out_dir, "bulk_primary_scores.csv"), row.names = FALSE)
write.csv(stats, file.path(out_dir, "bulk_ROC_AP_permutation_statistics.csv"), row.names = FALSE)
write.csv(pairwise, file.path(out_dir, "bulk_paired_bootstrap_AUC_comparisons.csv"), row.names = FALSE)
write.csv(leave_one_positive, file.path(out_dir, "bulk_leave_one_recurrence_case_out.csv"), row.names = FALSE)
write.csv(coverage, file.path(out_dir, "bulk_gene_set_coverage.csv"), row.names = FALSE)

writeLines(capture.output({
  cat("Samples:", length(y), "recurrence:", sum(y), "nonrecurrence:", sum(y == 0), "\n")
  cat("Prevalence AP baseline:", mean(y), "\n")
  print(stats)
  print(pairwise)
  cat("Leave-one-positive AUC ranges:\n")
  print(leave_one_positive %>% group_by(score) %>% summarise(min_auc = min(auc), max_auc = max(auc), min_ap = min(average_precision), max_ap = max(average_precision)))
  print(sessionInfo())
}), file.path(out_dir, "bulk_primary_analysis_report.txt"))

cat("bulk primary translation analysis complete\n")
print(stats)
print(pairwise)
