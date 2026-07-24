source(file.path(if (length(grep("^--file=", commandArgs(FALSE), value=TRUE))) dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value=TRUE)[1]), winslash="/")) else getwd(), "package_config.R"))
suppressPackageStartupMessages({
  library(monocle)
  library(Biobase)
  library(Seurat)
  library(Matrix)
  library(harmony)
  library(scDblFinder)
  library(SingleCellExperiment)
  library(RANN)
  library(cluster)
  library(dplyr)
})

set.seed(20260711)
out_dir <- file.path(audit_root, "QC_Harmony")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

adjusted_rand_index <- function(x, y) {
  tab <- table(x, y); ch2 <- function(z) z * (z - 1) / 2; n <- sum(tab)
  a <- sum(ch2(tab)); b <- sum(ch2(rowSums(tab))); c <- sum(ch2(colSums(tab)))
  e <- b * c / ch2(n); m <- (b + c) / 2
  if (m == e) 1 else (a - e) / (m - e)
}
neighbor_metrics <- function(embedding, sample, program, k = 30) {
  nn <- RANN::nn2(embedding, k = min(k + 1, nrow(embedding)))$nn.idx[, -1, drop = FALSE]
  entropy <- apply(nn, 1, function(idx) {
    p <- as.numeric(prop.table(table(sample[idx])))
    p <- p[p > 0]
    -sum(p * log(p)) / log(length(unique(sample)))
  })
  purity <- vapply(seq_len(nrow(nn)), function(i) mean(program[nn[i, ]] == program[i]), numeric(1))
  c(mean_sample_entropy = mean(entropy), mean_program_neighbor_purity = mean(purity))
}

cds <- readRDS(file.path(input_root, "discovery_trajectory_graph_pseudotime_cds.rds"))
pd <- pData(cds)
counts <- exprs(cds)
obj <- CreateSeuratObject(counts = counts, meta.data = pd)
obj$program <- pd$major_fibroblast_subtype
if (!"percent.mt" %in% colnames(obj@meta.data)) obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")

obj <- NormalizeData(obj, verbose = FALSE)
obj <- FindVariableFeatures(obj, nfeatures = 3000, verbose = FALSE)
regress_vars <- intersect(c("nCount_RNA", "percent.mt"), colnames(obj@meta.data))
obj <- ScaleData(obj, features = VariableFeatures(obj), vars.to.regress = regress_vars, verbose = FALSE)
obj <- RunPCA(obj, features = VariableFeatures(obj), npcs = 50, verbose = FALSE)

resolutions <- c(0.1, 0.2, 0.3, 0.4)
dims_values <- c(20, 30, 40)
sens <- list()
coords_out <- list()
for (nd in dims_values) {
  no <- FindNeighbors(obj, reduction = "pca", dims = seq_len(nd), graph.name = c(paste0("pca_nn_", nd), paste0("pca_snn_", nd)), verbose = FALSE)
  no <- FindClusters(no, graph.name = paste0("pca_snn_", nd), resolution = 0.3, algorithm = 1, random.seed = 20260711, verbose = FALSE)
  no_cluster <- as.character(Idents(no))
  no_met <- neighbor_metrics(Embeddings(no, "pca")[, seq_len(nd), drop = FALSE], no$sample, no$program)
  sens[[length(sens) + 1]] <- data.frame(method = "PCA_no_Harmony", dims = nd, resolution = 0.3, clusters = length(unique(no_cluster)), ARI_program = adjusted_rand_index(no_cluster, no$program), mean_sample_entropy = no_met[1], mean_program_neighbor_purity = no_met[2])

  ha <- RunHarmony(obj, group.by.vars = "sample", reduction.use = "pca", dims.use = seq_len(nd), reduction.save = paste0("harmony_", nd), plot_convergence = FALSE, verbose = FALSE)
  ha <- FindNeighbors(ha, reduction = paste0("harmony_", nd), dims = seq_len(nd), graph.name = c(paste0("h_nn_", nd), paste0("h_snn_", nd)), verbose = FALSE)
  ha <- FindClusters(ha, graph.name = paste0("h_snn_", nd), resolution = 0.3, algorithm = 1, random.seed = 20260711, verbose = FALSE)
  h_cluster <- as.character(Idents(ha))
  h_met <- neighbor_metrics(Embeddings(ha, paste0("harmony_", nd))[, seq_len(nd), drop = FALSE], ha$sample, ha$program)
  sens[[length(sens) + 1]] <- data.frame(method = "Harmony_sample", dims = nd, resolution = 0.3, clusters = length(unique(h_cluster)), ARI_program = adjusted_rand_index(h_cluster, ha$program), mean_sample_entropy = h_met[1], mean_program_neighbor_purity = h_met[2])

  if (nd == 30) {
    no <- RunUMAP(no, reduction = "pca", dims = 1:30, reduction.name = "umap_no_harmony", seed.use = 20260711, verbose = FALSE)
    ha <- RunUMAP(ha, reduction = "harmony_30", dims = 1:30, reduction.name = "umap_harmony", seed.use = 20260711, verbose = FALSE)
    u1 <- Embeddings(no, "umap_no_harmony"); u2 <- Embeddings(ha, "umap_harmony")
    coords_out[[1]] <- data.frame(cell = rownames(u1), method = "PCA_no_Harmony", UMAP_1 = u1[, 1], UMAP_2 = u1[, 2], sample = no$sample, condition = no$condition, program = no$program, cluster = no_cluster)
    coords_out[[2]] <- data.frame(cell = rownames(u2), method = "Harmony_sample", UMAP_1 = u2[, 1], UMAP_2 = u2[, 2], sample = ha$sample, condition = ha$condition, program = ha$program, cluster = h_cluster)
  }
}
sens <- bind_rows(sens)
write.csv(sens, file.path(out_dir, "Harmony_dimension_sensitivity_metrics.csv"), row.names = FALSE)
write.csv(bind_rows(coords_out), file.path(out_dir, "Harmony_before_after_UMAP_metadata.csv"), row.names = FALSE)

# Resolution sensitivity at 30 dimensions for both representations.
obj_h <- RunHarmony(obj, group.by.vars = "sample", reduction.use = "pca", dims.use = 1:30, reduction.save = "harmony_final", plot_convergence = FALSE, verbose = FALSE)
resolution_audit <- bind_rows(lapply(resolutions, function(res) {
  a <- FindNeighbors(obj, reduction = "pca", dims = 1:30, graph.name = c("tmp_p_nn", "tmp_p_snn"), verbose = FALSE)
  a <- FindClusters(a, graph.name = "tmp_p_snn", resolution = res, algorithm = 1, random.seed = 20260711, verbose = FALSE)
  b <- FindNeighbors(obj_h, reduction = "harmony_final", dims = 1:30, graph.name = c("tmp_h_nn", "tmp_h_snn"), verbose = FALSE)
  b <- FindClusters(b, graph.name = "tmp_h_snn", resolution = res, algorithm = 1, random.seed = 20260711, verbose = FALSE)
  bind_rows(
    data.frame(method = "PCA_no_Harmony", resolution = res, n_clusters = length(unique(Idents(a))), ARI_program = adjusted_rand_index(Idents(a), a$program)),
    data.frame(method = "Harmony_sample", resolution = res, n_clusters = length(unique(Idents(b))), ARI_program = adjusted_rand_index(Idents(b), b$program))
  )
}))
write.csv(resolution_audit, file.path(out_dir, "Harmony_resolution_sensitivity_metrics.csv"), row.names = FALSE)

# Leave-one-sample-out structure audit.
loo <- bind_rows(lapply(unique(obj$sample), function(drop) {
  sub <- subset(obj, cells = colnames(obj)[obj$sample != drop])
  sub <- RunPCA(sub, features = intersect(VariableFeatures(obj), rownames(sub)), npcs = 30, verbose = FALSE)
  sub_h <- RunHarmony(sub, group.by.vars = "sample", reduction.use = "pca", dims.use = 1:30, reduction.save = "harmony_loo", plot_convergence = FALSE, verbose = FALSE)
  get_one <- function(x, reduction, method) {
    x <- FindNeighbors(x, reduction = reduction, dims = 1:30, graph.name = c("loo_nn", "loo_snn"), verbose = FALSE)
    x <- FindClusters(x, graph.name = "loo_snn", resolution = 0.3, algorithm = 1, random.seed = 20260711, verbose = FALSE)
    met <- neighbor_metrics(Embeddings(x, reduction)[, 1:30, drop = FALSE], x$sample, x$program)
    data.frame(omitted_sample = drop, method = method, n_cells = ncol(x), n_clusters = length(unique(Idents(x))), ARI_program = adjusted_rand_index(Idents(x), x$program), mean_sample_entropy = met[1], mean_program_neighbor_purity = met[2])
  }
  bind_rows(get_one(sub, "pca", "PCA_no_Harmony"), get_one(sub_h, "harmony_loo", "Harmony_sample"))
}))
write.csv(loo, file.path(out_dir, "Harmony_leave_one_sample_out_structure_audit.csv"), row.names = FALSE)

# Mitochondrial-threshold sensitivity within the final trajectory population.
mt <- bind_rows(lapply(c(10, 15, 20, 25), function(th) {
  keep <- obj$percent.mt <= th
  meta <- obj@meta.data[keep, , drop = FALSE]
  bind_rows(lapply(split(seq_len(nrow(meta)), meta$sample), function(idx) data.frame(
    threshold_percent_mt = th,
    sample = meta$sample[idx[1]], condition = meta$condition[idx[1]],
    retained_cells = length(idx),
    retained_fraction = length(idx) / sum(obj$sample == meta$sample[idx[1]]),
    homeostatic_fraction = mean(meta$program[idx] == "Homeostatic/regulatory fibroblasts"),
    matrix_fraction = mean(meta$program[idx] == "Matrix-remodeling fibroblasts"),
    ecm_fraction = mean(meta$program[idx] == "ECM-producing fibroblasts")
  )))
}))
write.csv(mt, file.path(out_dir, "fibroblast_percent_mt_threshold_sensitivity.csv"), row.names = FALSE)

# Retrospective doublet audit in the 1,350 candidate-fibroblast object, stratified by sample.
cand <- readRDS(file.path(input_root, "discovery_fibroblast_subclustered.rds"))
DefaultAssay(cand) <- "RNA"
cand_counts <- GetAssayData(cand, assay = "RNA", layer = "counts")
cand_meta <- cand@meta.data
sce <- SingleCellExperiment(assays = list(counts = cand_counts), colData = S4Vectors::DataFrame(cand_meta))
sce <- scDblFinder(sce, samples = "sample", verbose = FALSE)
dbl <- data.frame(
  cell = colnames(sce), sample = colData(sce)$sample,
  scDblFinder_score = colData(sce)$scDblFinder.score,
  scDblFinder_class = colData(sce)$scDblFinder.class,
  included_in_trajectory_analysis = colnames(sce) %in% colnames(obj),
  stringsAsFactors = FALSE
)
write.csv(dbl, file.path(out_dir, "candidate_fibroblast_scDblFinder_audit.csv"), row.names = FALSE)
dbl_summary <- dbl %>% group_by(sample, included_in_trajectory_analysis) %>% summarise(n_cells = n(), predicted_doublets = sum(scDblFinder_class == "doublet"), predicted_doublet_fraction = mean(scDblFinder_class == "doublet"), .groups = "drop")
write.csv(dbl_summary, file.path(out_dir, "candidate_fibroblast_scDblFinder_summary.csv"), row.names = FALSE)

writeLines(capture.output({
  print(sens); print(resolution_audit); print(loo); print(dbl_summary); print(sessionInfo())
}), file.path(out_dir, "QC_Harmony_report.txt"))
cat("QC and Harmony sensitivity complete\n")
print(sens); print(dbl_summary)
