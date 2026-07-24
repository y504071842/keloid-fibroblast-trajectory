args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  stop("Usage: Rscript bulk_04_tximport_gene_symbol.R SALMON_QUANT_DIR GENCODE_V50_GTF OUTPUT_DIR")
}

suppressPackageStartupMessages({
  library(tximport)
  library(rtracklayer)
})

quant_dir <- normalizePath(args[[1]], winslash = "/", mustWork = TRUE)
gtf_path <- normalizePath(args[[2]], winslash = "/", mustWork = TRUE)
out_dir <- normalizePath(args[[3]], winslash = "/", mustWork = FALSE)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

quant_files <- list.files(quant_dir, pattern = "^quant\\.sf$", recursive = TRUE, full.names = TRUE)
names(quant_files) <- basename(dirname(quant_files))
quant_files <- quant_files[order(names(quant_files))]
stopifnot(length(quant_files) == 77L, all(file.exists(quant_files)))

gtf <- rtracklayer::import(gtf_path)
tx <- as.data.frame(gtf[gtf$type == "transcript"])
tx2gene <- unique(data.frame(TXNAME = tx$transcript_id, GENEID = tx$gene_name, stringsAsFactors = FALSE))
tx2gene <- tx2gene[complete.cases(tx2gene) & tx2gene$TXNAME != "" & tx2gene$GENEID != "", ]
write.csv(tx2gene, file.path(out_dir, "gencode_v50_tx2gene_symbol.csv"), row.names = FALSE, quote = FALSE)

txi <- tximport(
  files = quant_files,
  type = "salmon",
  tx2gene = tx2gene,
  countsFromAbundance = "no",
  ignoreTxVersion = FALSE
)

write_matrix <- function(x, filename) {
  write.csv(data.frame(Gene = rownames(x), x, check.names = FALSE),
            file.path(out_dir, filename), row.names = FALSE, quote = FALSE)
}

write_matrix(txi$abundance, "gene_symbol_TPM_matrix.csv")
write_matrix(txi$counts, "gene_symbol_est_counts_matrix.csv")
write_matrix(txi$length, "gene_symbol_effective_length_matrix.csv")
write_matrix(log2(txi$abundance + 1), "gene_symbol_log2TPM_plus1_matrix.csv")

cat("Generated gene-symbol matrices for", ncol(txi$abundance), "samples and", nrow(txi$abundance), "genes.\n")

