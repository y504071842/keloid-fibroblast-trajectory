script_path <- if (length(grep("^--file=", commandArgs(FALSE), value = TRUE))) {
  sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
} else {
  file.path(getwd(), "define_clean_programs_and_enrichment.R")
}
source(file.path(dirname(chartr("\\", "/", script_path)), "package_config.R"))
suppressPackageStartupMessages({
  library(dplyr)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(msigdbr)
})

set.seed(20260711)
out_dir <- file.path(audit_root, "clean_programs")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dynamic_source <- required_input(
  file.path(repo_root, "data", "derived", "dynamic_program_candidates.csv"),
  "Frozen dynamic-program candidate table"
)
dynamic <- read.csv(dynamic_source, stringsAsFactors = FALSE)
branch <- read.csv(file.path(audit_root, "branch_trajectory", "lineage_dynamic_gene_classification.csv"), stringsAsFactors = FALSE)

is_minimal_technical <- function(g) grepl("^MT-", g) | grepl("^(RP[SL]|MRP[SL])[0-9A-Z]", g) | grepl("^HB[ABDEGMQZ][0-9]", g)
is_histone_or_keratin <- function(g) grepl("^HIST[0-9]", g) | grepl("^H[1234]-", g) | grepl("^KRT[0-9]", g)

audit <- dynamic %>%
  mutate(
    expected_direction = ifelse(dynamic_module == "D1_early_decreasing", "decreasing", "increasing"),
    observed_direction = ifelse(spearman_rho < 0, "decreasing", "increasing"),
    direction_concordant = expected_direction == observed_direction,
    minimal_technical_flag = is_minimal_technical(gene),
    histone_or_keratin_flag = is_histone_or_keratin(gene),
    exclusion_reason = case_when(
      !direction_concordant ~ "direction_discordant_with_kmeans_program",
      minimal_technical_flag ~ "mitochondrial_ribosomal_or_hemoglobin",
      histone_or_keratin_flag ~ "histone_or_keratin_nonspecific",
      TRUE ~ "retained"
    ),
    clean_module = case_when(
      exclusion_reason == "retained" & dynamic_module == "D1_early_decreasing" ~ "clean_D1",
      exclusion_reason == "retained" & dynamic_module == "D2_late_increasing" ~ "clean_D2",
      TRUE ~ NA_character_
    )
  ) %>%
  left_join(branch[, c("gene", "branch_pattern")], by = "gene")

write.csv(audit, file.path(out_dir, "clean_D1D2_gene_level_audit.csv"), row.names = FALSE)
summary <- audit %>% count(dynamic_module, observed_direction, exclusion_reason, clean_module, name = "gene_n")
write.csv(summary, file.path(out_dir, "clean_D1D2_summary.csv"), row.names = FALSE)

conf <- as.data.frame.matrix(table(kmeans_program = dynamic$dynamic_module, spearman_direction = ifelse(dynamic$spearman_rho < 0, "negative", "positive"))) %>% tibble::rownames_to_column("kmeans_program")
write.csv(conf, file.path(out_dir, "kmeans_vs_Spearman_direction_confusion.csv"), row.names = FALSE)

sets <- list(
  clean_D1 = audit$gene[!is.na(audit$clean_module) & audit$clean_module == "clean_D1"],
  clean_D2 = audit$gene[!is.na(audit$clean_module) & audit$clean_module == "clean_D2"]
)
for (bp in c("common_increasing", "mixed_remodeling_branch_biased", "ecm_dominant_branch_biased", "common_decreasing")) {
  g <- audit$gene[audit$branch_pattern == bp & !is.na(audit$clean_module)]
  if (length(g) >= 5) sets[[paste0("branch_", bp)]] <- unique(g)
}
write.csv(do.call(rbind, lapply(names(sets), function(nm) data.frame(gene_set = nm, gene = sets[[nm]], stringsAsFactors = FALSE))), file.path(out_dir, "clean_D1D2_and_branch_gene_sets.csv"), row.names = FALSE)

universe <- unique(dynamic$gene)
id_map <- bitr(universe, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
universe_entrez <- unique(id_map$ENTREZID)

run_cp <- function(nm, genes) {
  ids <- unique(id_map$ENTREZID[id_map$SYMBOL %in% genes])
  go <- tryCatch(enrichGO(ids, OrgDb = org.Hs.eg.db, keyType = "ENTREZID", ont = "BP", universe = universe_entrez, pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1, readable = TRUE), error = function(e) NULL)
  list(
    GO_BP = if (is.null(go)) data.frame() else as.data.frame(go)
  )
}

ora <- lapply(names(sets), function(nm) {
  x <- run_cp(nm, sets[[nm]])
  bind_rows(lapply(names(x), function(db) {
    if (!nrow(x[[db]])) return(data.frame())
    cbind(gene_set = nm, database = db, x[[db]], stringsAsFactors = FALSE)
  }))
}) %>% bind_rows()
write.csv(ora, file.path(out_dir, "clean_D1D2_GO_KEGG_enrichment.csv"), row.names = FALSE)

coll <- msigdbr_collections()
hallmark <- msigdbr(species = "Homo sapiens", collection = "H") %>% dplyr::select(gs_name, gene_symbol)
msig_kegg <- msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "CP:KEGG_LEGACY") %>% dplyr::select(gs_name, gene_symbol)
reactome_row <- coll %>% filter(gs_collection == "C2", grepl("REACTOME", gs_subcollection, ignore.case = TRUE)) %>% slice_head(n = 1)
reactome <- if (nrow(reactome_row)) msigdbr(species = "Homo sapiens", collection = "C2", subcollection = reactome_row$gs_subcollection) %>% dplyr::select(gs_name, gene_symbol) else data.frame(gs_name = character(), gene_symbol = character())

run_symbol_ora <- function(nm, genes, term2gene, db) {
  if (!nrow(term2gene)) return(data.frame())
  e <- tryCatch(enricher(genes, universe = universe, TERM2GENE = term2gene, pvalueCutoff = 1, pAdjustMethod = "BH", qvalueCutoff = 1), error = function(e) NULL)
  if (is.null(e)) return(data.frame())
  cbind(gene_set = nm, database = db, as.data.frame(e), stringsAsFactors = FALSE)
}
msig_ora <- bind_rows(lapply(names(sets), function(nm) bind_rows(
  run_symbol_ora(nm, sets[[nm]], hallmark, "MSigDB_Hallmark"),
  run_symbol_ora(nm, sets[[nm]], msig_kegg, "MSigDB_KEGG_Legacy"),
  run_symbol_ora(nm, sets[[nm]], reactome, "MSigDB_Reactome")
)))
write.csv(msig_ora, file.path(out_dir, "clean_D1D2_Hallmark_Reactome_enrichment.csv"), row.names = FALSE)

writeLines(capture.output({
  print(summary)
  cat("Set sizes:\n"); print(vapply(sets, length, integer(1)))
  cat("Top GO:\n"); print(ora %>% arrange(p.adjust) %>% dplyr::select(gene_set, database, Description, p.adjust, Count) %>% head(30))
  cat("Top Hallmark/KEGG/Reactome:\n"); print(msig_ora %>% arrange(p.adjust) %>% dplyr::select(gene_set, database, Description, p.adjust, Count) %>% head(30))
  print(sessionInfo())
}), file.path(out_dir, "clean_program_report.txt"))
cat("objective module audit and enrichment complete\n")
print(summary)
