source(file.path(if (length(grep("^--file=", commandArgs(FALSE), value = TRUE))) dirname(chartr("\\", "/", sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))) else getwd(), "package_config.R"))
suppressPackageStartupMessages({
  library(monocle)
  library(Biobase)
  library(Matrix)
  library(Seurat)
  library(WGCNA)
  library(dplyr)
})

options(stringsAsFactors = FALSE)
allowWGCNAThreads(nThreads = 4)
set.seed(20260717)

primary_dir <- file.path(audit_root, Sys.getenv("WGCNA_PRIMARY_DIR", unset = "WGCNA"))
out_name <- Sys.getenv("WGCNA_SENS_OUT_DIR", unset = "WGCNA_PCA_metacell_sensitivity")
out_dir <- file.path(audit_root, out_name)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

cds <- readRDS(file.path(input_root, "discovery_trajectory_graph_pseudotime_cds.rds"))
obj <- readRDS(file.path(input_root, "discovery_fibroblast_subclustered.rds"))
pd <- pData(cds)
pd$sample <- as.character(pd$sample)
pd$condition <- as.character(pd$condition)
pd$major_fibroblast_subtype <- as.character(pd$major_fibroblast_subtype)
cells <- intersect(rownames(pd), colnames(obj))
if (length(cells) != nrow(pd)) stop("Not all trajectory cells were found in the discovery Seurat object.")

counts <- exprs(cds)
sf <- as.numeric(pd$Size_Factor); names(sf) <- rownames(pd)
logexpr <- log1p(t(t(counts) / sf[colnames(counts)]))
pca <- Embeddings(obj, reduction = "pca")[rownames(pd), seq_len(min(20L, ncol(Embeddings(obj, reduction = "pca")))), drop = FALSE]

clean <- read.csv(file.path(audit_root, "clean_programs", "clean_D1D2_gene_level_audit.csv"), stringsAsFactors = FALSE)
d1 <- unique(clean$gene[clean$clean_module == "clean_D1" & !is.na(clean$clean_module)])
d2 <- unique(clean$gene[clean$clean_module == "clean_D2" & !is.na(clean$clean_module)])
score_mean <- function(mat, genes) Matrix::colMeans(mat[intersect(genes, rownames(mat)), , drop = FALSE])

cell_meta <- data.frame(
  cell = rownames(pd), sample = pd$sample, condition = pd$condition,
  program = pd$major_fibroblast_subtype, pseudotime = as.numeric(pd$Pseudotime_graph),
  clean_D1 = score_mean(logexpr, d1)[rownames(pd)], clean_D2 = score_mean(logexpr, d2)[rownames(pd)],
  stringsAsFactors = FALSE
)

groups <- split(cell_meta$cell, paste(cell_meta$sample, cell_meta$program, sep = "__"))
assignments <- lapply(names(groups), function(group_name) {
  group_cells <- groups[[group_name]]
  n <- length(group_cells)
  if (n < 8L) return(NULL)
  centers <- max(1L, round(n / 25))
  centers <- min(centers, n)
  if (centers == 1L) {
    cluster <- rep(1L, n)
  } else {
    set.seed(20260717 + match(group_name, sort(names(groups))))
    cluster <- kmeans(pca[group_cells, , drop = FALSE], centers = centers, nstart = 100, iter.max = 200)$cluster
    repeat {
      sizes <- table(cluster)
      small <- as.integer(names(sizes)[sizes < 8L])
      if (!length(small) || length(unique(cluster)) == 1L) break
      centroids <- rowsum(pca[group_cells, , drop = FALSE], cluster) / as.numeric(table(cluster))
      for (cl in small) {
        other <- setdiff(as.integer(rownames(centroids)), cl)
        if (!length(other)) next
        distance <- rowSums((centroids[as.character(other), , drop = FALSE] - matrix(centroids[as.character(cl), ], nrow = length(other), ncol = ncol(centroids), byrow = TRUE))^2)
        cluster[cluster == cl] <- other[which.min(distance)]
      }
      cluster <- match(cluster, sort(unique(cluster)))
    }
  }
  data.frame(cell = group_cells, pca_cluster = cluster, stringsAsFactors = FALSE)
})
assignments <- do.call(rbind, assignments)
assignments <- merge(assignments, cell_meta, by = "cell", all.x = TRUE, sort = FALSE)
assignments$metacell <- paste(assignments$sample, assignments$program, paste0("PCA", assignments$pca_cluster), sep = "__")

mc_meta <- assignments %>%
  group_by(metacell, sample, condition, program, pca_cluster) %>%
  summarise(n_cells = n(), mean_pseudotime = mean(pseudotime), clean_D1 = mean(clean_D1), clean_D2 = mean(clean_D2), .groups = "drop") %>%
  filter(n_cells >= 8)

selection <- read.csv(file.path(primary_dir, "WGCNA_discovery_only_gene_selection.csv"), stringsAsFactors = FALSE)
gene_pool <- selection$gene[selection$selected_top_variable %in% TRUE]
gene_pool <- intersect(gene_pool, rownames(logexpr))
mc_expr <- sapply(mc_meta$metacell, function(mc) {
  mc_cells <- assignments$cell[assignments$metacell == mc]
  Matrix::rowMeans(logexpr[gene_pool, mc_cells, drop = FALSE])
})
datExpr <- t(as.matrix(mc_expr))
good <- goodSamplesGenes(datExpr, verbose = 0)
datExpr <- datExpr[good$goodSamples, good$goodGenes, drop = FALSE]
mc_meta <- mc_meta[match(rownames(datExpr), mc_meta$metacell), , drop = FALSE]

primary_parameters <- read.csv(file.path(primary_dir, "WGCNA_parameters.csv"), stringsAsFactors = FALSE)
soft_power <- as.numeric(primary_parameters$value[primary_parameters$parameter == "soft_power"])
net <- blockwiseModules(
  datExpr, power = soft_power, networkType = "signed", TOMType = "signed",
  minModuleSize = 25, mergeCutHeight = 0.25, numericLabels = TRUE,
  pamRespectsDendro = FALSE, maxBlockSize = ncol(datExpr), verbose = 2
)
colors <- labels2colors(net$colors); names(colors) <- colnames(datExpr)
MEs <- orderMEs(moduleEigengenes(datExpr, colors = colors)$eigengenes)

weighted_rank_cor <- function(x, y, w) {
  keep <- is.finite(x) & is.finite(y) & is.finite(w) & w > 0
  x <- x[keep]; y <- y[keep]; w <- w[keep]
  if (length(x) < 4L || sd(x) == 0 || sd(y) == 0) return(NA_real_)
  xr <- rank(x); yr <- rank(y)
  mx <- sum(w * xr) / sum(w); my <- sum(w * yr) / sum(w)
  sum(w * (xr - mx) * (yr - my)) / sqrt(sum(w * (xr - mx)^2) * sum(w * (yr - my)^2))
}
donor_equal_cor <- function(x, y, donor) {
  tab <- table(donor); weighted_rank_cor(x, y, 1 / as.numeric(tab[donor]))
}

traits <- mc_meta[, c("clean_D1", "clean_D2")]
modules <- setdiff(unique(colors), "grey")
trait_cor <- do.call(rbind, lapply(modules, function(module) do.call(rbind, lapply(colnames(traits), function(trait) {
  data.frame(module = module, trait = trait, donor_equal_rho = donor_equal_cor(MEs[[paste0("ME", module)]], traits[[trait]], mc_meta$sample))
}))))
selected <- trait_cor %>% group_by(trait) %>% slice_max(donor_equal_rho, n = 1, with_ties = FALSE) %>% ungroup()

primary_modules <- read.csv(file.path(primary_dir, "WGCNA_gene_module_assignment.csv"), stringsAsFactors = FALSE)
primary_selected <- read.csv(file.path(primary_dir, "WGCNA_selected_aligned_modules.csv"), stringsAsFactors = FALSE)
module_concordance <- do.call(rbind, lapply(seq_len(nrow(primary_selected)), function(i) {
  program <- primary_selected$program[[i]]
  primary_module <- primary_selected$aligned_module[[i]]
  primary_genes <- primary_modules$gene[primary_modules$module == primary_module]
  overlaps <- vapply(modules, function(module) sum(primary_genes %in% names(colors)[colors == module]), numeric(1))
  best <- names(which.max(overlaps))
  aligned <- selected$module[selected$trait == program]
  aligned_genes <- names(colors)[colors == aligned]
  best_genes <- names(colors)[colors == best]
  aligned_overlap <- sum(primary_genes %in% aligned_genes)
  data.frame(
    program = program, primary_module = primary_module,
    pca_metacell_aligned_module = aligned,
    pca_metacell_aligned_rho = selected$donor_equal_rho[selected$trait == program],
    aligned_module_overlap_n = aligned_overlap,
    aligned_module_jaccard = aligned_overlap / length(union(primary_genes, aligned_genes)),
    pca_metacell_best_overlap_module = best,
    best_overlap_module_rho = trait_cor$donor_equal_rho[trait_cor$trait == program & trait_cor$module == best],
    best_overlap_n = max(overlaps),
    best_overlap_jaccard = max(overlaps) / length(union(primary_genes, best_genes)),
    same_module_by_trait_and_overlap = aligned == best,
    stringsAsFactors = FALSE
  )
}))

write.csv(assignments, file.path(out_dir, "PCA_metacell_cell_assignments.csv"), row.names = FALSE)
write.csv(mc_meta, file.path(out_dir, "PCA_metacell_composition.csv"), row.names = FALSE)
write.csv(data.frame(gene = names(colors), module = unname(colors)), file.path(out_dir, "PCA_metacell_gene_module_assignment.csv"), row.names = FALSE)
write.csv(trait_cor, file.path(out_dir, "PCA_metacell_module_trait_correlations.csv"), row.names = FALSE)
write.csv(selected, file.path(out_dir, "PCA_metacell_selected_aligned_modules.csv"), row.names = FALSE)
write.csv(module_concordance, file.path(out_dir, "PCA_metacell_primary_module_concordance.csv"), row.names = FALSE)
write.csv(data.frame(
  parameter = c("seed", "embedding", "pseudotime_used_for_metacell_construction", "grouping", "target_cells_per_metacell", "minimum_cells", "network_gene_space", "soft_power"),
  value = c(20260717, "first 20 discovery PCA dimensions", FALSE, "sample x fibroblast program", 25, 8, "same discovery-only top-variable genes as primary network", soft_power)
), file.path(out_dir, "PCA_metacell_parameters.csv"), row.names = FALSE)
saveRDS(list(datExpr = datExpr, colors = colors, MEs = MEs, mc_meta = mc_meta, net = net, selected = selected), file.path(out_dir, "PCA_metacell_WGCNA_result.rds"), compress = FALSE)
writeLines(capture.output({
  cat("Pseudotime-independent PCA-metacell sensitivity analysis\n")
  cat("Metacells:", nrow(datExpr), "represented donors:", length(unique(mc_meta$sample)), "\n")
  print(selected)
  print(module_concordance)
  print(sessionInfo())
}), file.path(out_dir, "PCA_metacell_report.txt"))

cat("PCA-metacell WGCNA sensitivity complete\n")
print(selected)
print(module_concordance)
