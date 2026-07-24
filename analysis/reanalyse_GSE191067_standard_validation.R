source(file.path(if (length(grep("^--file=", commandArgs(FALSE), value=TRUE))) dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value=TRUE)[1]), winslash="/")) else getwd(), "package_config.R"))
suppressPackageStartupMessages({
  library(Seurat)
  library(harmony)
  library(Matrix)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(cluster)
})

set.seed(20260712)

out_dir <- file.path(audit_root, "GSE191067_standard_validation")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

obj_path <- external_object(
  "KFT_GSE191067_OBJECT",
  file.path(repo_root, "data", "frozen_inputs", "GSE191067_standard_reclustered_scored.rds"),
  "Processed GSE191067 fibroblast Seurat object"
)
audit_path <- file.path(audit_root, "clean_programs", "clean_D1D2_gene_level_audit.csv")

obj <- readRDS(obj_path)
DefaultAssay(obj) <- "RNA"
obj$sample <- as.character(obj$sample)
obj$condition <- ifelse(grepl("^HK", obj$sample), "Keloid", "Normal")

# Standard full-HVG reclustering of the purified validation fibroblasts.
obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 3000, verbose = FALSE)
obj <- ScaleData(obj, features = VariableFeatures(obj), vars.to.regress = c("percent.mt", "nCount_RNA"), verbose = FALSE)
obj <- RunPCA(obj, features = VariableFeatures(obj), npcs = 30, seed.use = 20260712, verbose = FALSE)
obj <- RunHarmony(obj, group.by.vars = "sample", reduction = "pca", dims.use = 1:30,
                  reduction.save = "standard_harmony", plot_convergence = FALSE, verbose = FALSE)
obj <- FindNeighbors(obj, reduction = "standard_harmony", dims = 1:30,
                     graph.name = c("standard_nn", "standard_snn"), k.param = 20, verbose = FALSE)
obj <- RunUMAP(obj, reduction = "standard_harmony", dims = 1:30, reduction.name = "standard_umap",
               reduction.key = "sUMAP_", n.neighbors = 30, min.dist = 0.3,
               metric = "cosine", seed.use = 20260712, verbose = FALSE)

resolutions <- c(0.02, 0.05, 0.08, 0.10, 0.15, 0.20, 0.30)
set.seed(20260712)
sil_cells <- sample(colnames(obj), min(4000, ncol(obj)))
sil_emb <- Embeddings(obj, "standard_harmony")[sil_cells, 1:15, drop = FALSE]
sil_dist <- dist(sil_emb)
resolution_audit <- list()
for (res in resolutions) {
  cname <- paste0("standard_res_", gsub("\\.", "_", sprintf("%.2f", res)))
  obj <- FindClusters(obj, graph.name = "standard_snn", resolution = res, algorithm = 1,
                      random.seed = 20260712, cluster.name = cname, verbose = FALSE)
  labs <- as.character(obj@meta.data[sil_cells, cname])
  sil <- if (length(unique(labs)) > 1) mean(cluster::silhouette(as.integer(factor(labs)), sil_dist)[, "sil_width"]) else NA_real_
  sizes <- table(obj@meta.data[[cname]])
  resolution_audit[[length(resolution_audit) + 1]] <- data.frame(
    resolution = res, cluster_column = cname, n_clusters = length(sizes),
    min_cluster_n = min(sizes), median_cluster_n = median(as.numeric(sizes)),
    mean_silhouette = sil, stringsAsFactors = FALSE
  )
}
resolution_audit <- bind_rows(resolution_audit)
eligible <- resolution_audit %>% filter(n_clusters >= 3, n_clusters <= 10, min_cluster_n >= 50, is.finite(mean_silhouette))
if (nrow(eligible) == 0) stop("No eligible standard reclustering resolution.")
chosen <- eligible %>% arrange(desc(mean_silhouette), resolution) %>% slice(1)
cluster_col <- chosen$cluster_column[[1]]
obj$standard_cluster <- as.character(obj@meta.data[[cluster_col]])

marker_sets <- list(
  `Homeostatic/regulatory fibroblasts` = c("PI16", "DPT", "APOD", "CFD", "IGFBP3", "WIF1", "BMP7", "ABCA8", "ADH1B", "C7"),
  `Matrix-remodeling fibroblasts` = c("NGFR", "THBS4", "HAS2", "TIMP1", "COL6A5", "F13A1", "LOXL3", "CTSC", "FOXS1", "CHST1"),
  `ECM-producing fibroblasts` = c("COL1A1", "COL1A2", "COL3A1", "COL5A1", "COMP", "ASPN", "TAGLN", "CTHRC1", "LRRC15", "POSTN", "SPARC", "SERPINH1")
)

x <- GetAssayData(obj, assay = "RNA", layer = "data")
score_gene_set <- function(genes) {
  present <- intersect(genes, rownames(x))
  mat <- as.matrix(x[present, , drop = FALSE])
  z <- t(scale(t(mat)))
  z[!is.finite(z)] <- 0
  as.numeric(colMeans(z))
}

marker_score_names <- character()
for (nm in names(marker_sets)) {
  col <- paste0("marker_", make.names(nm), "_score")
  obj[[col]] <- score_gene_set(marker_sets[[nm]])
  marker_score_names[nm] <- col
}

marker_meta <- data.frame(
  cell = rownames(obj@meta.data),
  standard_cluster = obj$standard_cluster,
  obj@meta.data[, unname(marker_score_names), drop = FALSE],
  check.names = FALSE,
  stringsAsFactors = FALSE
)
marker_long <- marker_meta %>%
  pivot_longer(cols = all_of(unname(marker_score_names)), names_to = "score_column", values_to = "marker_score") %>%
  mutate(program = names(marker_score_names)[match(score_column, marker_score_names)])
cluster_marker_scores <- marker_long %>%
  group_by(standard_cluster, program) %>%
  summarise(mean_marker_score = mean(marker_score), median_marker_score = median(marker_score), .groups = "drop")
cluster_mapping <- cluster_marker_scores %>%
  group_by(standard_cluster) %>%
  arrange(desc(mean_marker_score), program, .by_group = TRUE) %>%
  mutate(rank = row_number(), margin_to_second = mean_marker_score - lead(mean_marker_score)) %>%
  slice(1) %>%
  ungroup() %>%
  select(standard_cluster, standard_program = program, top_mean_marker_score = mean_marker_score, margin_to_second)
obj$standard_program <- cluster_mapping$standard_program[match(obj$standard_cluster, cluster_mapping$standard_cluster)]

# Freeze and score the current 42-gene clean D1 and 975-gene clean D2 definitions.
audit <- read.csv(audit_path, stringsAsFactors = FALSE)
clean_d1 <- unique(audit$gene[!is.na(audit$clean_module) & audit$clean_module == "clean_D1"])
clean_d2 <- unique(audit$gene[!is.na(audit$clean_module) & audit$clean_module == "clean_D2"])
obj$clean_D1_score <- score_gene_set(clean_d1)
obj$clean_D2_score <- score_gene_set(clean_d2)
obj$clean_D1_loss <- -obj$clean_D1_score

program_order <- c("Homeostatic/regulatory fibroblasts", "Matrix-remodeling fibroblasts", "ECM-producing fibroblasts")
obj$standard_program <- factor(obj$standard_program, levels = program_order)

umap <- as.data.frame(Embeddings(obj, "standard_umap"))
colnames(umap) <- c("UMAP_1", "UMAP_2")
meta <- bind_cols(
  data.frame(cell = rownames(obj@meta.data), stringsAsFactors = FALSE),
  obj@meta.data %>% select(sample, condition, standard_cluster, standard_program,
                           clean_D1_score, clean_D2_score, clean_D1_loss),
  umap
)

clean_coverage <- data.frame(
  program = c("clean_D1", "clean_D2"),
  defined_n = c(length(clean_d1), length(clean_d2)),
  detected_n = c(length(intersect(clean_d1, rownames(x))), length(intersect(clean_d2, rownames(x)))),
  coverage = c(length(intersect(clean_d1, rownames(x))) / length(clean_d1),
               length(intersect(clean_d2, rownames(x))) / length(clean_d2))
)

state_cell_summary <- meta %>%
  group_by(standard_program) %>%
  summarise(n_cells = n(), median_D1 = median(clean_D1_score), median_D2 = median(clean_D2_score), .groups = "drop")
donor_state_scores <- meta %>%
  group_by(sample, condition, standard_program) %>%
  summarise(n_cells = n(), median_D1 = median(clean_D1_score), median_D2 = median(clean_D2_score), .groups = "drop")
sample_fractions <- meta %>% count(sample, condition, standard_program, name = "n_cells") %>%
  group_by(sample) %>% mutate(fraction = n_cells / sum(n_cells)) %>% ungroup()

# Directional summaries are descriptive; no cell-level inferential P values are generated.
state_direction <- state_cell_summary %>%
  mutate(state_rank = match(standard_program, program_order) - 1) %>%
  summarise(
    D1_spearman_rho = cor(state_rank, median_D1, method = "spearman"),
    D2_spearman_rho = cor(state_rank, median_D2, method = "spearman"),
    D1_homeostatic_to_ECM_change = median_D1[state_rank == 2] - median_D1[state_rank == 0],
    D2_homeostatic_to_ECM_change = median_D2[state_rank == 2] - median_D2[state_rank == 0]
  )

write.csv(resolution_audit, file.path(out_dir, "standard_reclustering_resolution_audit.csv"), row.names = FALSE)
write.csv(cluster_marker_scores, file.path(out_dir, "standard_cluster_marker_scores.csv"), row.names = FALSE)
write.csv(cluster_mapping, file.path(out_dir, "standard_cluster_to_program_mapping.csv"), row.names = FALSE)
write.csv(data.frame(program = rep(names(marker_sets), lengths(marker_sets)), gene = unlist(marker_sets),
                     detected = unlist(marker_sets) %in% rownames(x)),
          file.path(out_dir, "standard_annotation_marker_genes.csv"), row.names = FALSE)
write.csv(clean_coverage, file.path(out_dir, "standard_clean_D1D2_gene_coverage.csv"), row.names = FALSE)
write.csv(meta, file.path(out_dir, "standard_validation_cell_metadata.csv"), row.names = FALSE)
write.csv(state_cell_summary, file.path(out_dir, "standard_state_cell_score_summary.csv"), row.names = FALSE)
write.csv(donor_state_scores, file.path(out_dir, "standard_donor_state_scores.csv"), row.names = FALSE)
write.csv(sample_fractions, file.path(out_dir, "standard_program_fractions_by_sample.csv"), row.names = FALSE)
write.csv(state_direction, file.path(out_dir, "standard_state_direction_summary.csv"), row.names = FALSE)
unlink(file.path(out_dir, c("standard_disease_sample_scores.csv", "standard_disease_effect_summary.csv",
                            "D2_disease_statistical_method_audit.csv",
                            "standard_state_composition_disease_effects.csv")))
saveRDS(obj, file.path(out_dir, "GSE191067_standard_reclustered_scored.rds"), compress = FALSE)

report <- capture.output({
  cat("Chosen resolution\n"); print(chosen)
  cat("\nCluster mapping\n"); print(cluster_mapping)
  cat("\nProgram counts\n"); print(table(obj$standard_program, useNA = "ifany"))
  cat("\nClean coverage\n"); print(clean_coverage)
  cat("\nState cell summaries\n"); print(state_cell_summary)
  cat("\nState direction\n"); print(state_direction)
  cat("\nProgram fractions by donor\n"); print(sample_fractions)
  cat("\nSession info\n"); print(sessionInfo())
})
writeLines(report, file.path(out_dir, "standard_validation_report.txt"))
cat(paste(report, collapse = "\n"), "\n")
