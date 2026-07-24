source(file.path(if (length(grep("^--file=", commandArgs(FALSE), value=TRUE))) dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value=TRUE)[1]), winslash="/")) else getwd(), "package_config.R"))
suppressPackageStartupMessages({
  library(dplyr)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(msigdbr)
})
set.seed(20260711)
out_dir <- file.path(audit_root, "program_annotation")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
markers <- read.csv(file.path(audit_root, "program_annotation", "fibroblast_program_top50_positive_markers.csv"), stringsAsFactors = FALSE)
markers <- markers %>% filter(p_val_adj < 0.05, avg_log2FC > 0.10, pct.1 >= 0.10) %>% arrange(cluster, p_val_adj, desc(avg_log2FC))
program_genes <- split(markers$gene, markers$cluster)
program_genes <- lapply(program_genes, function(x) head(unique(x), 300))
top_markers <- markers %>% group_by(cluster) %>% distinct(gene, .keep_all = TRUE) %>% slice_head(n = 50) %>% ungroup()
write.csv(top_markers, file.path(out_dir, "fibroblast_program_top50_positive_markers.csv"), row.names = FALSE)
write.csv(do.call(rbind, lapply(names(program_genes), function(nm) data.frame(program = nm, gene = program_genes[[nm]], rank = seq_along(program_genes[[nm]])))), file.path(out_dir, "fibroblast_program_marker_gene_sets.csv"), row.names = FALSE)

universe <- unique(markers$gene)
id_map <- bitr(universe, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
universe_entrez <- unique(id_map$ENTREZID)
go <- bind_rows(lapply(names(program_genes), function(nm) {
  ids <- unique(id_map$ENTREZID[id_map$SYMBOL %in% program_genes[[nm]]])
  e <- enrichGO(ids, OrgDb = org.Hs.eg.db, keyType = "ENTREZID", ont = "BP", universe = universe_entrez, pvalueCutoff = 1, qvalueCutoff = 1, pAdjustMethod = "BH", readable = TRUE)
  cbind(program = nm, database = "GO_BP", as.data.frame(e), stringsAsFactors = FALSE)
}))
write.csv(go, file.path(out_dir, "fibroblast_program_GO_BP_enrichment.csv"), row.names = FALSE)

hallmark <- msigdbr(species = "Homo sapiens", collection = "H") %>% dplyr::select(gs_name, gene_symbol)
kegg <- msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "CP:KEGG_LEGACY") %>% dplyr::select(gs_name, gene_symbol)
reactome <- msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "CP:REACTOME") %>% dplyr::select(gs_name, gene_symbol)
run_msig <- function(nm, term2gene, db) {
  e <- enricher(program_genes[[nm]], universe = universe, TERM2GENE = term2gene, pvalueCutoff = 1, qvalueCutoff = 1, pAdjustMethod = "BH")
  cbind(program = nm, database = db, as.data.frame(e), stringsAsFactors = FALSE)
}
msig <- bind_rows(lapply(names(program_genes), function(nm) bind_rows(
  run_msig(nm, hallmark, "MSigDB_Hallmark"),
  run_msig(nm, kegg, "MSigDB_KEGG_Legacy"),
  run_msig(nm, reactome, "MSigDB_Reactome")
)))
write.csv(msig, file.path(out_dir, "fibroblast_program_Hallmark_KEGG_Reactome_enrichment.csv"), row.names = FALSE)

cross <- data.frame(gene = unlist(program_genes), stringsAsFactors = FALSE) %>% count(gene, name = "program_membership_count") %>% arrange(desc(program_membership_count), gene)
write.csv(cross, file.path(out_dir, "fibroblast_program_marker_cross_program_frequency.csv"), row.names = FALSE)
writeLines(capture.output({
  cat("Program set sizes:\n"); print(vapply(program_genes, length, integer(1)))
  cat("Top GO:\n"); print(go %>% arrange(program, p.adjust) %>% dplyr::select(program, Description, p.adjust, Count) %>% group_by(program) %>% slice_head(n = 10))
  cat("Top MSigDB:\n"); print(msig %>% arrange(program, database, p.adjust) %>% dplyr::select(program, database, Description, p.adjust, Count) %>% group_by(program, database) %>% slice_head(n = 5))
  print(sessionInfo())
}), file.path(out_dir, "program_annotation_report.txt"))
cat("program annotation enrichment complete\n")
