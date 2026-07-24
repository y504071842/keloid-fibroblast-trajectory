source(file.path(if (length(grep("^--file=", commandArgs(FALSE), value = TRUE))) dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]), winslash = "/")) else getwd(), "package_config.R"))
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(scales)
  library(pROC)
})
set.seed(20260711)
root <- package_root
main_root <- file.path(root, "01_Main_Figures")
supp_root <- file.path(root, "02_Supplementary_Figures")
dir.create(main_root, recursive = TRUE, showWarnings = FALSE)
dir.create(supp_root, recursive = TRUE, showWarnings = FALSE)

save_panel <- function(p, stem, w, h) {
  dir.create(dirname(stem), recursive = TRUE, showWarnings = FALSE)
  cairo_png <- function(filename, width, height, bg = "white", ...) {
    grDevices::png(filename, width = width, height = height, units = "in", res = 450,
                   type = "cairo", bg = bg)
  }
  ggsave(paste0(stem, ".png"), p, width = w, height = h, units = "in", device = cairo_png, bg = "white")
  ggsave(paste0(stem, ".pdf"), p, width = w, height = h, units = "in", device = function(filename, width, height, ...) grDevices::pdf(filename, width = width, height = height, useDingbats = FALSE, bg = "white"))
}
theme_pub <- function(base = 10.5) theme_classic(base_size = base) + theme(plot.title = element_blank(), plot.subtitle = element_blank(), axis.text = element_text(color = "#171717"), strip.background = element_rect(fill = "#F0F0F0", color = NA), strip.text = element_text(face = "bold"), legend.title = element_text(size = base - 1))
# Figure 1: discovery atlas and fibroblast states, unchanged analytically.
f1 <- file.path(main_root, "Figure1"); dir.create(f1, recursive = TRUE, showWarnings = FALSE)
workflow <- data.frame(
  x = c(1, 2, 3, 3, 2, 1), y = c(2, 2, 2, 1, 1, 1),
  heading = c("Discovery scRNA-seq", "Purified fibroblasts", "Branch-aware trajectory",
              "Clean-program freezing", "External single-cell evaluation", "Exploratory bulk translation"),
  detail = c("3 keloids; 3 normal skin\n70,759 post-QC cells\nmajor-lineage annotation",
             "1,350 candidates\nboundary-cluster audit\n1,117 purified cells",
             "state-marker-anchored root\ndonor-aware robustness\ndonor-adjusted tradeSeq",
             "direction concordance\ntechnical gene-symbol rules\nfrozen D1: 42; D2: 975 genes",
             "GSE191067; 17,020 fibroblasts\nstate localization\nnetwork preservation",
             "PRJNA813172; n = 77\n7 recurrence events\nD1-loss candidate feature"),
  stringsAsFactors = FALSE
)
p1a <- ggplot(workflow) +
  geom_rect(aes(xmin = x - .39, xmax = x + .39, ymin = y - .31, ymax = y + .31),
            fill = "#F6F8F7", color = "#78958E", linewidth = .45) +
  geom_text(aes(x, y + .16, label = heading), fontface = "bold", size = 3.25, color = "#233B37") +
  geom_text(aes(x, y - .07, label = detail), size = 2.55, lineheight = .98, color = "#4A5553") +
  annotate("segment", x = 1.40, xend = 1.60, y = 2, yend = 2, arrow = arrow(length = grid::unit(.12, "in")), color = "#56756E") +
  annotate("segment", x = 2.40, xend = 2.60, y = 2, yend = 2, arrow = arrow(length = grid::unit(.12, "in")), color = "#56756E") +
  annotate("segment", x = 3, xend = 3, y = 1.68, yend = 1.32, arrow = arrow(length = grid::unit(.12, "in")), color = "#56756E") +
  annotate("segment", x = 2.60, xend = 2.40, y = 1, yend = 1, arrow = arrow(length = grid::unit(.12, "in")), color = "#56756E") +
  annotate("segment", x = 1.60, xend = 1.40, y = 1, yend = 1, arrow = arrow(length = grid::unit(.12, "in")), color = "#56756E") +
  coord_cartesian(xlim = c(.52, 3.48), ylim = c(.55, 2.45), clip = "off") + theme_void()
save_panel(p1a, file.path(f1, "Figure1A_study_design_workflow"), 10.0, 4.8)

# Figure 2: condition-blind trajectory, donor-level inference, branching dynamics.
f2 <- file.path(main_root, "Figure2"); dir.create(f2, recursive = TRUE, showWarnings = FALSE)

donor <- read.csv(file.path(root, "results/discovery_donor_level_summary.csv"), stringsAsFactors = FALSE)
donor$condition <- factor(donor$condition, levels = c("Normal", "Keloid"))
cond_cols <- c(Normal = "#4C78A8", Keloid = "#C44E52")
normal_pt <- donor$median_old_pseudotime[donor$condition == "Normal"]
keloid_pt <- donor$median_old_pseudotime[donor$condition == "Keloid"]
pt_delta <- mean(as.vector(outer(keloid_pt, normal_pt, FUN = function(x, y) sign(x - y))))
p2c <- ggplot(donor, aes(condition, median_old_pseudotime, color = condition)) +
  stat_summary(fun = median, fun.min = median, fun.max = median, geom = "crossbar", width = .48, linewidth = .55, color = "#333333") +
  geom_point(size = 2.8, position = position_jitter(width = .07, seed = 11)) +
  annotate("text", x = 1.5, y = max(donor$median_old_pseudotime) * 1.08,
           label = sprintf("Cliff's delta = %.2f", pt_delta), size = 3.4, color = "#222222") +
  scale_color_manual(values = cond_cols, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(.03, .16))) +
  labs(x = NULL, y = "Donor median graph pseudotime") + theme_pub()
save_panel(p2c, file.path(f2, "Figure2C_donor_level_pseudotime"), 4.4, 4.2)


branch <- read.csv(file.path(root, "results/branch_trajectory/lineage_dynamic_gene_classification.csv"), stringsAsFactors = FALSE)
branch_counts <- branch %>% count(original_module, branch_pattern) %>% group_by(original_module) %>% mutate(frac = n / sum(n)) %>% ungroup()
branch_labels <- c(common_decreasing = "Common decreasing", common_increasing = "Common increasing", ecm_dominant_branch_biased = "ECM-dominant branch biased", mixed_remodeling_branch_biased = "Mixed-remodeling branch biased", other = "Other")
branch_cols <- c(common_decreasing = "#2F8F83", common_increasing = "#B84A62", mixed_remodeling_branch_biased = "#D39A32", ecm_dominant_branch_biased = "#6F5AA8", other = "#BDBDBD")
p2f <- ggplot(branch_counts, aes(original_module, frac, fill = branch_pattern)) + geom_col(width = .62, color = "white", linewidth = .2) +
  scale_fill_manual(values = branch_cols, labels = branch_labels, name = "Two-lineage pattern") + scale_y_continuous(labels = percent) +
  scale_x_discrete(labels = c(D1_early_decreasing = "D1", D2_late_increasing = "D2")) + labs(x = NULL, y = "Gene fraction") + theme_pub() + theme(legend.position = "right")
save_panel(p2f, file.path(f2, "Figure2G_two_lineage_dynamic_composition"), 6.0, 4.4)

trade_files <- c(association = "tradeSeq_two_lineage_association_global.csv", pattern = "tradeSeq_two_lineage_patternTest.csv", terminal = "tradeSeq_two_lineage_diffEndTest.csv", early = "tradeSeq_two_lineage_earlyDETest.csv")
trade_sum <- bind_rows(lapply(names(trade_files), function(nm) {
  x <- read.csv(file.path(root, "results/branch_trajectory", trade_files[[nm]]), stringsAsFactors = FALSE)
  data.frame(test = nm, tested = sum(is.finite(x$fdr)), significant = sum(is.finite(x$fdr) & x$fdr < .05), fraction = sum(is.finite(x$fdr) & x$fdr < .05) / sum(is.finite(x$fdr)))
}))
p2g <- ggplot(trade_sum, aes(test, fraction, fill = test)) + geom_col(width = .62, color = "#333333", linewidth = .25) + geom_text(aes(label = paste0(significant, "/", tested)), vjust = -.35, size = 3) +
  scale_y_continuous(labels = percent, limits = c(0, min(1, max(trade_sum$fraction) * 1.18 + .03))) + scale_fill_manual(values = c(association = "#356D8C", pattern = "#8D5B8C", terminal = "#B76A45", early = "#5E8C61"), guide = "none") +
  labs(x = NULL, y = "FDR < 0.05 fraction") + theme_pub() + theme(axis.text.x = element_text(angle = 25, hjust = 1))
save_panel(p2g, file.path(f2, "Figure2H_two_lineage_tradeSeq_tests"), 5.4, 4.3)

sling <- read.csv(file.path(root, "results/branch_trajectory/slingshot_DDRTree_pseudotime_concordance.csv"), stringsAsFactors = FALSE)
sling$label <- ifelse(sling$mode == "start_only_unsupervised", "Start-only", paste("Forced", sling$lineage))
p2h <- ggplot(sling, aes(label, spearman_with_DDRTree_graph_pseudotime, fill = mode)) + geom_col(width = .62, color = "#333333", linewidth = .25) + geom_text(aes(label = sprintf("rho = %.2f", spearman_with_DDRTree_graph_pseudotime)), vjust = -.35, size = 3) +
  scale_fill_manual(values = c(forced_two_terminal = "#D39A32", start_only_unsupervised = "#4C78A8"), guide = "none") + ylim(0, 1) + labs(x = NULL, y = "Spearman correlation with DDRTree") + theme_pub() + theme(axis.text.x = element_text(angle = 20, hjust = 1))
save_panel(p2h, file.path(f2, "Figure2D_cross_method_Slingshot_sensitivity"), 5.3, 4.3)

# Figure 3: clean-program audit, functional structure and discovery-only WGCNA.
f3 <- file.path(main_root, "Figure3"); dir.create(f3, recursive = TRUE, showWarnings = FALSE)
audit <- read.csv(file.path(root, "results/clean_programs/clean_D1D2_gene_level_audit.csv"), stringsAsFactors = FALSE)
funnel <- bind_rows(
  data.frame(program = "D1", step = factor(c("k-means", "Direction-concordant", "Clean"), levels = c("k-means", "Direction-concordant", "Clean")), n = c(sum(audit$dynamic_module == "D1_early_decreasing"), sum(audit$dynamic_module == "D1_early_decreasing" & audit$direction_concordant), sum(audit$clean_module == "clean_D1", na.rm = TRUE))),
  data.frame(program = "D2", step = factor(c("k-means", "Direction-concordant", "Clean"), levels = c("k-means", "Direction-concordant", "Clean")), n = c(sum(audit$dynamic_module == "D2_late_increasing"), sum(audit$dynamic_module == "D2_late_increasing" & audit$direction_concordant), sum(audit$clean_module == "clean_D2", na.rm = TRUE)))
)
p3a <- ggplot(funnel, aes(step, n, group = program, color = program)) + geom_line(linewidth = 1) + geom_point(size = 2.5) + geom_text(aes(label = n), vjust = -0.6, size = 3) + facet_wrap(~program, scales = "free_y") + scale_color_manual(values = c(D1 = "#2F8F83", D2 = "#B84A62"), guide = "none") + labs(x = NULL, y = "Genes") + theme_pub() + theme(axis.text.x = element_text(angle = 20, hjust = 1))
save_panel(p3a, file.path(f3, "Figure3A_clean_program_audit"), 6.2, 4.2)

excl <- audit %>% filter(exclusion_reason != "retained") %>% count(dynamic_module, exclusion_reason) %>%
  mutate(reason = recode(exclusion_reason,
    direction_discordant_with_kmeans_program = "Direction discordance",
    histone_or_keratin_nonspecific = "Histone or keratin",
    mitochondrial_ribosomal_or_hemoglobin = "Mito/ribo/HB"))
p3b <- ggplot(excl, aes(reason, n, fill = dynamic_module)) + geom_col(position = "dodge", width = .68) + scale_fill_manual(values = c(D1_early_decreasing = "#2F8F83", D2_late_increasing = "#B84A62"), labels = c("D1", "D2"), name = NULL) + labs(x = NULL, y = "Excluded genes") + theme_pub() + theme(axis.text.x = element_text(angle = 15, hjust = 1), legend.position = "top")
save_panel(p3b, file.path(f3, "Figure3B_clean_exclusion_categories"), 6.4, 4.2)

wdir <- file.path(root, "results/WGCNA")
w_sens <- read.csv(file.path(wdir, "WGCNA_selected_module_donor_sensitivity.csv"), stringsAsFactors = FALSE) %>%
  mutate(pair_label = paste0(ifelse(program == "clean_D1", "Clean D1", "Clean D2"), " / ", aligned_module),
         pair = factor(pair_label, levels = rev(pair_label)))
p3c <- ggplot(w_sens, aes(y = pair)) +
  geom_errorbar(aes(xmin = donor_equal_cluster_bootstrap_q025, xmax = donor_equal_cluster_bootstrap_q975),
                orientation = "y", width = .16, linewidth = .65, color = "#777777") +
  geom_point(aes(x = donor_equal_rho, fill = pair), shape = 21, size = 4.0, color = "white", stroke = .65) +
  geom_point(aes(x = donor_equal_lodo_min_rho), shape = 17, size = 2.8, color = "#222222") +
  scale_fill_manual(values = setNames(c("#2F8F83", "#B84A62"), w_sens$pair_label), guide = "none") +
  scale_x_continuous(limits = c(.65, 1), breaks = seq(.7, 1, .1)) +
  labs(x = "Module-program Spearman correlation", y = NULL) + theme_pub() +
  annotate("text", x = .665, y = 2.34, hjust = 0,
           label = "Circle: donor-equal estimate | Triangle: minimum leave-one-donor-out estimate\nLine: donor-cluster bootstrap 95% interval",
           size = 2.65, color = "#444444") +
  theme(axis.line.y = element_blank(), axis.ticks.y = element_blank())
save_panel(p3c, file.path(f3, "Figure3F_donor_aware_WGCNA_sensitivity"), 7.0, 4.2)

ov <- read.csv(file.path(wdir, "WGCNA_clean_D1D2_overlap.csv"), stringsAsFactors = FALSE)
ov <- ov %>% mutate(aligned_fraction = overlap_n / set_n_in_network,
                    label = paste0(overlap_n, "/", set_n_in_network, " network-represented genes\n",
                                   module, " module"))
p3d <- ggplot(ov, aes(program, aligned_fraction, fill = program)) +
  geom_col(width = .58, color = "#333333", linewidth = .25) +
  geom_text(aes(label = label), vjust = -.25, size = 2.8, lineheight = .95) +
  scale_fill_manual(values = c(clean_D1 = "#2F8F83", clean_D2 = "#B84A62"), guide = "none") +
  scale_x_discrete(labels = c(clean_D1 = "Clean D1", clean_D2 = "Clean D2")) +
  scale_y_continuous(labels = percent, limits = c(0, 1.12)) +
  labs(x = NULL, y = "Fraction assigned to aligned module") + theme_pub()
save_panel(p3d, file.path(f3, "Figure3G_WGCNA_dynamic_program_overlap"), 4.8, 4.2)

obj_branch <- audit %>% filter(!is.na(clean_module)) %>% count(clean_module, branch_pattern) %>% group_by(clean_module) %>% mutate(frac = n / sum(n)) %>% ungroup()
p3e <- ggplot(obj_branch, aes(clean_module, frac, fill = branch_pattern)) + geom_col(width = .62, color = "white") + scale_fill_manual(values = branch_cols, labels = branch_labels, name = "Branch pattern") + scale_y_continuous(labels = percent) + scale_x_discrete(labels = c(clean_D1 = "Clean D1", clean_D2 = "Clean D2")) + labs(x = NULL, y = "Gene fraction") + theme_pub()
save_panel(p3e, file.path(f3, "Figure3D_clean_program_branch_structure"), 6.2, 4.4)

go <- read.csv(file.path(root, "results/clean_programs/clean_D1D2_GO_KEGG_enrichment.csv"), stringsAsFactors = FALSE)
go_show <- go %>% filter(gene_set %in% c("clean_D1", "clean_D2"), is.finite(p.adjust)) %>% group_by(gene_set) %>% arrange(p.adjust, .by_group = TRUE) %>% slice_head(n = 6) %>% ungroup() %>% mutate(term = reorder(Description, -log10(p.adjust)))
p3f <- ggplot(go_show, aes(-log10(p.adjust), term, size = Count, color = gene_set)) + geom_point(alpha = .9) + scale_color_manual(values = c(clean_D1 = "#2F8F83", clean_D2 = "#B84A62"), labels = c("Clean D1", "Clean D2"), name = NULL) + labs(x = "-log10(FDR)", y = NULL, size = "Genes") + theme_pub() + theme(legend.position = "right")
save_panel(p3f, file.path(f3, "Figure3C_clean_program_GO_enrichment"), 7.4, 5.0)

# Figure 4: marker-annotated external evaluation and matched-unit network preservation.
f4 <- file.path(main_root, "Figure4"); dir.create(f4, recursive = TRUE, showWarnings = FALSE)
unlink(list.files(f4, pattern = "^Figure4.*\\.(png|pdf)$", full.names = TRUE))
std_dir <- file.path(root, "results/GSE191067_standard_validation")
hm <- read.csv(file.path(std_dir, "standard_validation_cell_metadata.csv"), stringsAsFactors = FALSE)
prog_cols <- c("Homeostatic/regulatory fibroblasts" = "#2F8F83", "Matrix-remodeling fibroblasts" = "#B68D40", "ECM-producing fibroblasts" = "#B84A62")
program_order <- c("Homeostatic/regulatory fibroblasts", "Matrix-remodeling fibroblasts", "ECM-producing fibroblasts")
short_labels <- c("Homeostatic/regulatory fibroblasts" = "Homeostatic", "Matrix-remodeling fibroblasts" = "Matrix-remodeling", "ECM-producing fibroblasts" = "ECM-producing")
hm$standard_program <- factor(hm$standard_program, levels = program_order)
hm$program_short <- factor(unname(short_labels[as.character(hm$standard_program)]), levels = unname(short_labels[program_order]))
p4a <- ggplot(hm, aes(UMAP_1, UMAP_2, color = standard_program)) + geom_point(size = .24, alpha = .82, stroke = 0) + scale_color_manual(values = prog_cols, name = NULL) + coord_equal() + labs(x = "Standard UMAP 1", y = "Standard UMAP 2") + theme_pub() + theme(legend.position = "right")
save_panel(p4a, file.path(f4, "Figure4A_standard_GSE191067_state_UMAP"), 6.4, 5.0)

marker_stats <- read.csv(file.path(std_dir, "standard_marker_dotplot_statistics.csv"), stringsAsFactors = FALSE)
marker_show <- c("WIF1", "PI16", "ADH1B", "CFD", "APOD", "CTSC", "NGFR", "FOXS1", "CHST1", "HAS2", "ASPN", "CTHRC1", "SPARC", "COL1A1", "POSTN")
marker_stats <- marker_stats %>% filter(gene %in% marker_show) %>% mutate(gene = factor(gene, levels = rev(marker_show)), standard_program = factor(standard_program, levels = program_order))
p4b <- ggplot(marker_stats, aes(standard_program, gene, size = pct_expressing, color = mean_gene_z)) + geom_point() + scale_size(range = c(.5, 6), name = "% expressing") + scale_color_gradient2(low = "#3A6EA5", mid = "white", high = "#B84A62", midpoint = 0, name = "Mean gene z") + scale_x_discrete(labels = short_labels) + labs(x = NULL, y = NULL) + theme_pub() + theme(axis.text.x = element_text(angle = 18, hjust = 1), legend.position = "right")
save_panel(p4b, file.path(f4, "Figure4B_standard_state_marker_dotplot"), 7.2, 5.8)

donor_state <- read.csv(file.path(std_dir, "standard_donor_state_scores.csv"), stringsAsFactors = FALSE) %>%
  mutate(program_short = factor(unname(short_labels[standard_program]), levels = unname(short_labels[program_order])),
         condition = factor(condition, levels = c("Normal", "Keloid")))
make_state_score_plot <- function(score_col, stem, ylab) {
  summary_data <- donor_state %>%
    group_by(program_short) %>%
    summarise(donor_median = median(.data[[score_col]]),
              q1 = quantile(.data[[score_col]], .25), q3 = quantile(.data[[score_col]], .75), .groups = "drop") %>%
    mutate(label = sprintf("Median %.2f", donor_median))
  label_offset <- diff(range(donor_state[[score_col]])) * .09
  p <- ggplot(donor_state, aes(program_short, .data[[score_col]])) +
    geom_errorbar(data = summary_data, aes(x = program_short, ymin = q1, ymax = q3), width = .12, linewidth = .55, color = "#222222", inherit.aes = FALSE) +
    geom_point(data = summary_data, aes(x = program_short, y = donor_median), shape = 23, size = 3.4, fill = "white", color = "#111111", stroke = .8, inherit.aes = FALSE) +
    geom_point(aes(fill = program_short, size = n_cells), shape = 21, color = "white", stroke = .35,
               position = position_jitter(width = .10, height = 0, seed = 42), alpha = .92) +
    geom_text(data = summary_data, aes(x = program_short, y = q3 + label_offset, label = label), size = 2.8, vjust = 0, inherit.aes = FALSE) +
    scale_fill_manual(values = setNames(unname(prog_cols[program_order]), unname(short_labels[program_order])), guide = "none") +
    scale_size_continuous(range = c(2.0, 4.8), name = "Cells") +
    scale_y_continuous(expand = expansion(mult = c(.08, .20))) +
    labs(x = NULL, y = ylab) + theme_pub() +
    theme(axis.text.x = element_text(angle = 15, hjust = 1), legend.position = "right")
  save_panel(p, file.path(f4, stem), 5.8, 4.6)
}
make_state_score_plot("median_D1", "Figure4C_clean_D1_donor_state_summary", "Donor median clean D1 score")
make_state_score_plot("median_D2", "Figure4D_clean_D2_donor_state_summary", "Donor median clean D2 score")

frac <- read.csv(file.path(std_dir, "standard_program_fractions_by_sample.csv"), stringsAsFactors = FALSE)
# The external evaluation includes three keloid and three normal-skin samples.
frac$sample <- factor(frac$sample, levels = c("HNS1", "HNS2", "HNS3", "HK1", "HK2", "HK3"))
frac$standard_program <- factor(frac$standard_program, levels = rev(program_order))
p4e <- ggplot(frac, aes(sample, fraction, fill = standard_program)) + geom_col(width = .72, color = "white", linewidth = .2) + scale_fill_manual(values = prog_cols[rev(program_order)], name = NULL) + scale_y_continuous(labels = percent) + labs(x = NULL, y = "Fibroblast state fraction") + theme_pub() + theme(axis.text.x = element_text(angle = 20, hjust = 1), legend.position = "right")
save_panel(p4e, file.path(f4, "Figure4E_state_composition_by_sample"), 6.5, 4.5)

pres <- read.csv(file.path(std_dir, "standard_module_preservation_2000_permutations.csv"), stringsAsFactors = FALSE)
ps <- pres[!is.na(pres$program_alignment) & pres$program_alignment != "NA", , drop = FALSE]
zcol <- grep("Zsummary", names(ps), value = TRUE)[1]; rankcol <- grep("medianRank", names(ps), value = TRUE)[1]
ps$Zsummary <- ps[[zcol]]; ps$medianRank <- ps[[rankcol]]
ps$program <- factor(ifelse(ps$program_alignment == "clean_D1-aligned", "Clean D1-aligned", "Clean D2-aligned"),
                     levels = c("Clean D2-aligned", "Clean D1-aligned"))
ps$label <- paste0("Zsummary = ", sprintf("%.2f", ps$Zsummary),
                   " | medianRank = ", ps$medianRank, "\nThreshold met")
ps$label_x <- ifelse(ps$Zsummary < 8, ps$Zsummary + .45, ps$Zsummary - .45)
ps$label_hjust <- ifelse(ps$Zsummary < 8, 0, 1)
p4f <- ggplot(ps, aes(Zsummary, program, color = program)) +
  annotate("rect", xmin = 2, xmax = 15.5, ymin = -Inf, ymax = Inf, fill = "#E8F2ED", alpha = .75) +
  geom_vline(xintercept = 2, linetype = "dashed", color = "#3F6E62", linewidth = .65) +
  geom_point(size = 5.0, shape = 21, aes(fill = program), color = "white", stroke = .7) +
  geom_text(aes(x = label_x, label = label, hjust = label_hjust), nudge_y = .22,
            lineheight = .95, size = 2.8, color = "#222222") +
  annotate("text", x = 2.25, y = 2.47, label = "Preservation evidence threshold (Zsummary > 2)",
           hjust = 0, size = 2.9, color = "#315B52") +
  scale_color_manual(values = c("Clean D1-aligned" = "#2F8F83", "Clean D2-aligned" = "#B84A62"), guide = "none") +
  scale_fill_manual(values = c("Clean D1-aligned" = "#2F8F83", "Clean D2-aligned" = "#B84A62"), guide = "none") +
  scale_x_continuous(breaks = c(0, 2, 5, 10, 15), limits = c(0, 15.5), expand = expansion(mult = c(.01, .02))) +
  labs(x = "Module preservation Zsummary", y = NULL) + theme_pub() +
  theme(panel.grid.major.x = element_line(color = "#E5E5E5", linewidth = .3),
        axis.line.y = element_blank(), axis.ticks.y = element_blank())
save_panel(p4f, file.path(f4, "Figure4F_module_preservation_threshold"), 7.0, 4.1)

# Figure 5: exploratory recurrence translation with discrimination and calibration-aware robustness.
f5 <- file.path(main_root, "Figure5"); dir.create(f5, recursive = TRUE, showWarnings = FALSE)
scores <- read.csv(file.path(root, "results/bulk_recurrence/bulk_primary_scores.csv"), stringsAsFactors = FALSE)
stats <- read.csv(file.path(root, "results/bulk_recurrence/bulk_ROC_AP_permutation_statistics.csv"), stringsAsFactors = FALSE)
scores$status <- factor(ifelse(scores$recurrence == 1, "Recurrence", "Non-recurrence"), levels = c("Non-recurrence", "Recurrence"))
status_cols <- c("Non-recurrence" = "#4C78A8", "Recurrence" = "#C44E52")
cohort <- as.data.frame(table(scores$status)); names(cohort) <- c("status", "n")
p5a <- ggplot(cohort, aes(status, n, fill = status)) + geom_col(width = .58, color = "#333333", linewidth = .25) + geom_text(aes(label = n), vjust = -.35, size = 3.2) + scale_fill_manual(values = status_cols, guide = "none") + labs(x = NULL, y = "Samples") + theme_pub()
save_panel(p5a, file.path(f5, "Figure5A_recurrence_cohort_composition"), 4.2, 4.0)

cov <- read.csv(file.path(root, "results/bulk_recurrence/bulk_gene_set_coverage.csv"), stringsAsFactors = FALSE) %>% mutate(label = recode(gene_set, clean_D1 = "Clean D1", clean_D2 = "Clean D2", five_hub = "Five-hub"), fraction = detected_n / defined_n)
p5b <- ggplot(cov, aes(label, fraction, fill = label)) + geom_col(width = .6, color = "#333333", linewidth = .25) + geom_text(aes(label = paste0(detected_n, "/", defined_n)), vjust = -.35, size = 3) + scale_y_continuous(labels = percent, limits = c(0, 1.08)) + scale_fill_manual(values = c("Clean D1" = "#2F8F83", "Clean D2" = "#B84A62", "Five-hub" = "#6F5AA8"), guide = "none") + labs(x = NULL, y = "Bulk gene coverage") + theme_pub() + theme(axis.text.x = element_text(angle = 15, hjust = 1))
save_panel(p5b, file.path(f5, "Figure5B_bulk_gene_coverage"), 4.8, 4.0)

make_bulk_box <- function(score_col, stem, ylab) {
  p <- ggplot(scores, aes(status, .data[[score_col]], fill = status, color = status)) + geom_boxplot(width = .48, outlier.shape = NA, alpha = .22, color = "#333333", linewidth = .3) + geom_point(position = position_jitter(width = .13, seed = 21), size = 1.4, alpha = .78) + scale_fill_manual(values = status_cols, guide = "none") + scale_color_manual(values = status_cols, guide = "none") + labs(x = NULL, y = ylab) + theme_pub()
  save_panel(p, file.path(f5, stem), 4.5, 4.2)
}
make_bulk_box("clean_D1_loss", "Figure5C_clean_D1_loss_by_recurrence", "Clean D1 loss score")
make_bulk_box("clean_D2", "Figure5D_clean_D2_by_recurrence", "Clean D2 score")

roc_scores <- c("clean_D1_loss" = "Clean D1 loss", "clean_D2" = "Clean D2", "five_hub" = "Five-hub")
roc_df <- bind_rows(lapply(names(roc_scores), function(nm) {
  r <- pROC::roc(scores$recurrence, scores[[nm]], direction = "<", quiet = TRUE)
  data.frame(score = roc_scores[[nm]], fpr = 1 - r$specificities, sensitivity = r$sensitivities)
}))
curve_cols <- c("Clean D1 loss" = "#2F8F83", "Clean D2" = "#B84A62", "Five-hub" = "#6F5AA8")
p5e <- ggplot(roc_df, aes(fpr, sensitivity, color = score, group = score)) + geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "#999999") + geom_path(linewidth = .95, lineend = "butt", linejoin = "mitre") + scale_color_manual(values = curve_cols, name = NULL) + coord_equal() + labs(x = "1 - specificity", y = "Sensitivity") + theme_pub() + theme(legend.position = "bottom")
save_panel(p5e, file.path(f5, "Figure5E_ROC_curves"), 5.2, 4.8)

# Average-precision step representation: each horizontal interval ends when one
# additional recurrent sample is retrieved, so the plotted area equals AP.
pr_df <- bind_rows(lapply(names(roc_scores), function(nm) {
  o <- order(scores[[nm]], decreasing = TRUE)
  yy <- scores$recurrence[o]
  positive_ranks <- which(yy == 1)
  positive_n <- length(positive_ranks)
  positive_precision <- seq_len(positive_n) / positive_ranks
  curve <- data.frame(recall = 0, precision = 1)
  curve <- rbind(curve, data.frame(recall = 0, precision = positive_precision[1]))
  for (j in seq_len(positive_n)) {
    curve <- rbind(curve, data.frame(recall = j / positive_n, precision = positive_precision[j]))
    if (j < positive_n) {
      curve <- rbind(curve, data.frame(recall = j / positive_n, precision = positive_precision[j + 1]))
    }
  }
  transform(curve, score = roc_scores[[nm]])
}))
p5f <- ggplot(pr_df, aes(recall, precision, color = score, group = score)) + geom_hline(yintercept = mean(scores$recurrence), linetype = "dashed", color = "#999999") + geom_path(linewidth = .95, lineend = "butt", linejoin = "mitre") + scale_color_manual(values = curve_cols, name = NULL) + coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) + labs(x = "Recall", y = "Precision") + theme_pub() + theme(legend.position = "bottom")
save_panel(p5f, file.path(f5, "Figure5F_precision_recall_curves"), 5.2, 4.8)

primary_names <- c(clean_D1_loss = "Clean D1 loss", clean_D2 = "Clean D2", five_hub = "Five-hub")
primary_auc <- stats %>% filter(score %in% names(primary_names)) %>% mutate(label = factor(unname(primary_names[score]), levels = rev(unname(primary_names))))
p5g <- ggplot(primary_auc, aes(auc, label, color = label)) + geom_errorbarh(aes(xmin = auc_ci_low, xmax = auc_ci_high), height = .18, color = "#666666") + geom_point(size = 2.8) + geom_vline(xintercept = .5, linetype = "dashed", color = "#888888") + scale_color_manual(values = curve_cols, guide = "none") + scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, .25)) + labs(x = "ROC AUC (95% bootstrap CI)", y = NULL) + theme_pub()
save_panel(p5g, file.path(f5, "Figure5G_primary_AUC_estimates"), 6.0, 3.8)

loo <- read.csv(file.path(root, "results/bulk_recurrence/bulk_leave_one_recurrence_case_out.csv"), stringsAsFactors = FALSE) %>% filter(score %in% c("clean_D1_loss", "five_hub")) %>% mutate(label = recode(score, clean_D1_loss = "Clean D1 loss", five_hub = "Five-hub"))
p5h <- ggplot(loo, aes(label, auc, color = label)) + geom_boxplot(width = .45, outlier.shape = NA, color = "#333333") + geom_point(position = position_jitter(width = .08, seed = 22), size = 2) + geom_hline(yintercept = .5, linetype = "dashed", color = "#888888") + scale_color_manual(values = c("Clean D1 loss" = "#2F8F83", "Five-hub" = "#6F5AA8"), guide = "none") + labs(x = NULL, y = "Leave-one-recurrence-out AUC") + theme_pub()
save_panel(p5h, file.path(f5, "Figure5H_leave_one_recurrence_case_out"), 4.8, 4.2)

# Supplementary Figure S1: QC and integration sensitivity.
s1 <- file.path(supp_root, "Supplementary_Figure_S1"); dir.create(s1, recursive = TRUE, showWarnings = FALSE)
hmeta <- read.csv(file.path(root, "results/QC_Harmony/Harmony_before_after_UMAP_metadata.csv"), stringsAsFactors = FALSE)
pS1c <- ggplot(hmeta, aes(UMAP_1, UMAP_2, color = sample)) + geom_point(size = .25, alpha = .75) + facet_wrap(~method) + coord_equal() + labs(x = "UMAP 1", y = "UMAP 2", color = "Sample") + theme_pub(9.5)
save_panel(pS1c, file.path(s1, "Supplementary_Figure_S1C_Harmony_before_after_by_sample"), 8.0, 4.2)
dbl <- read.csv(file.path(root, "results/QC_Harmony/candidate_fibroblast_scDblFinder_summary.csv"), stringsAsFactors = FALSE) %>% filter(included_in_trajectory_analysis)
pS1d <- ggplot(dbl, aes(sample, predicted_doublet_fraction, fill = sample)) + geom_col(width = .65, color = "#333333", linewidth = .2) + scale_y_continuous(labels = percent) + guides(fill = "none") + labs(x = NULL, y = "Predicted doublet fraction") + theme_pub() + theme(axis.text.x = element_text(angle = 25, hjust = 1))
save_panel(pS1d, file.path(s1, "Supplementary_Figure_S1D_scDblFinder_audit"), 5.8, 4.0)
mtx <- read.csv(file.path(root, "results/QC_Harmony/fibroblast_percent_mt_threshold_sensitivity.csv"), stringsAsFactors = FALSE)
pS1e <- ggplot(mtx, aes(threshold_percent_mt, retained_fraction, color = sample)) + geom_line() + geom_point() + scale_y_continuous(labels = percent) + scale_x_continuous(breaks = c(10,15,20,25)) + labs(x = "Mitochondrial threshold (%)", y = "Trajectory-cell retention", color = "Sample") + theme_pub()
save_panel(pS1e, file.path(s1, "Supplementary_Figure_S1E_percent_mt_sensitivity"), 6.0, 4.2)

# Supplementary Figure S2: trajectory robustness.
s2 <- file.path(supp_root, "Supplementary_Figure_S2"); dir.create(s2, recursive = TRUE, showWarnings = FALSE)
os <- read.csv(file.path(root, "results/ordering_gene_sensitivity/ordering_gene_trajectory_sensitivity_summary.csv"), stringsAsFactors = FALSE)
pS2b <- ggplot(os, aes(ordering_set, spearman_with_current_graph_pseudotime, fill = ordering_set)) + geom_col(width = .65, color = "#333333", linewidth = .2) + geom_text(aes(label = sprintf("%.2f", spearman_with_current_graph_pseudotime)), vjust = -.35, size = 3) + ylim(0,1.05) + guides(fill = "none") + labs(x = NULL, y = "Spearman correlation with primary ordering") + theme_pub() + theme(axis.text.x = element_text(angle = 25, hjust = 1))
save_panel(pS2b, file.path(s2, "Supplementary_Figure_S2B_ordering_gene_sensitivity"), 6.5, 4.3)

donor_trend <- read.csv(file.path(root, "results/donor_aware_trajectory/module_score_within_donor_trends.csv"), stringsAsFactors = FALSE) %>%
  dplyr::select(donor, clean_D1_rho, clean_D2_rho) %>%
  pivot_longer(cols = c(clean_D1_rho, clean_D2_rho), names_to = "program", values_to = "rho") %>%
  mutate(program = recode(program, clean_D1_rho = "Clean D1", clean_D2_rho = "Clean D2"))
pS2d <- ggplot(donor_trend, aes(rho, donor, color = program)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "#777777", linewidth = .45) +
  geom_point(size = 2.6, position = position_dodge(width = .45)) +
  scale_color_manual(values = c("Clean D1" = "#2F8F83", "Clean D2" = "#B84A62"), name = NULL) +
  scale_x_continuous(limits = c(-1, 1), breaks = c(-1, -.5, 0, .5, 1)) +
  labs(x = "Within-donor Spearman correlation with pseudotime", y = "Donor") + theme_pub() +
  theme(legend.position = "top")
save_panel(pS2d, file.path(s2, "Supplementary_Figure_S2D_within_donor_dynamic_program_trends"), 6.2, 4.6)

robust <- read.csv(file.path(root, "results/donor_aware_trajectory/gene_set_robustness_summary.csv"), stringsAsFactors = FALSE) %>%
  filter(gene_set %in% c("clean_D1", "clean_D2"),
         analysis %in% c("Per-donor majority direction among genes evaluable in >=3 donors",
                         "Fixed-trajectory LODO retained in >=5/6",
                         "Equal-donor downsampling: 6 donors x 11 cells; expected direction in >=80/100 repeats",
                         "GAM: donor_fixed, FDR<0.05 and expected endpoint direction",
                         "GAM: donor_random, FDR<0.05 and expected endpoint direction",
                         "GAM: donor_equal_weight, FDR<0.05 and expected endpoint direction",
                         "tradeSeq with donor covariates, association FDR<0.05")) %>%
  mutate(method = recode(analysis,
    `Per-donor majority direction among genes evaluable in >=3 donors` = "Per-donor majority",
    `Fixed-trajectory LODO retained in >=5/6` = "Leave-one-donor-out",
    `Equal-donor downsampling: 6 donors x 11 cells; expected direction in >=80/100 repeats` = "Equal-donor subsampling (>=80/100)",
    `GAM: donor_fixed, FDR<0.05 and expected endpoint direction` = "Donor-fixed GAM",
    `GAM: donor_random, FDR<0.05 and expected endpoint direction` = "Donor-random GAM",
    `GAM: donor_equal_weight, FDR<0.05 and expected endpoint direction` = "Donor-equal GAM",
    `tradeSeq with donor covariates, association FDR<0.05` = "Donor-adjusted tradeSeq"),
    program = recode(gene_set, clean_D1 = "Clean D1", clean_D2 = "Clean D2"),
    support = ifelse(grepl("Per-donor majority", analysis), supported_fraction_evaluable, supported_fraction_total),
    method = factor(method, levels = rev(c("Per-donor majority", "Leave-one-donor-out", "Equal-donor subsampling (>=80/100)",
                                          "Donor-fixed GAM", "Donor-random GAM", "Donor-equal GAM", "Donor-adjusted tradeSeq"))))
pS2e <- ggplot(robust, aes(support, method, color = program)) +
  geom_point(size = 2.8, position = position_dodge(width = .48)) +
  scale_color_manual(values = c("Clean D1" = "#2F8F83", "Clean D2" = "#B84A62"), name = NULL) +
  scale_x_continuous(labels = percent, limits = c(0, 1.02)) +
  labs(x = "Fraction of genes meeting the method-specific criterion", y = NULL) + theme_pub() +
  theme(legend.position = "top")
save_panel(pS2e, file.path(f3, "Figure3E_donor_aware_gene_support"), 7.0, 5.0)

pca_sens <- read.csv(file.path(root, "results/WGCNA_PCA_metacell_sensitivity/PCA_metacell_primary_module_concordance.csv"), stringsAsFactors = FALSE) %>%
  transmute(
    program = recode(program, clean_D1 = "Clean D1", clean_D2 = "Clean D2"),
    `Primary pseudotime-bin module` = ifelse(program == "Clean D1", w_sens$donor_equal_rho[w_sens$program == "clean_D1"], w_sens$donor_equal_rho[w_sens$program == "clean_D2"]),
    `PCA-metacell score-aligned module` = pca_metacell_aligned_rho,
    `PCA-metacell best-overlap module` = best_overlap_module_rho
  ) %>% pivot_longer(-program, names_to = "analysis", values_to = "rho")
pS2f <- ggplot(pca_sens, aes(rho, program, shape = analysis, color = program)) +
  geom_point(size = 3.0, position = position_dodge(width = .45)) +
  scale_color_manual(values = c("Clean D1" = "#2F8F83", "Clean D2" = "#B84A62"), guide = "none") +
  scale_shape_manual(values = c("Primary pseudotime-bin module" = 16, "PCA-metacell score-aligned module" = 17, "PCA-metacell best-overlap module" = 15), name = NULL) +
  scale_x_continuous(limits = c(.65, 1), breaks = seq(.7, 1, .1)) +
  labs(x = "Donor-equal module-program Spearman correlation", y = NULL) + theme_pub() +
  theme(legend.position = "bottom")
save_panel(pS2f, file.path(s2, "Supplementary_Figure_S2E_pseudotime_independent_WGCNA_sensitivity"), 7.2, 4.3)

# Supplementary Figure S3: external single-cell reclustering audit.
s3 <- file.path(supp_root, "Supplementary_Figure_S3"); dir.create(s3, recursive = TRUE, showWarnings = FALSE)
pS3a <- ggplot(hm, aes(UMAP_1, UMAP_2, color = sample)) + geom_point(size = .22, alpha = .78, stroke = 0) + coord_equal() + labs(x = "Standard UMAP 1", y = "Standard UMAP 2", color = "Sample") + theme_pub() + theme(legend.position = "right")
save_panel(pS3a, file.path(s3, "Supplementary_Figure_S3A_standard_UMAP_by_sample"), 6.3, 5.0)
pS3b <- ggplot(hm, aes(UMAP_1, UMAP_2, color = factor(standard_cluster))) + geom_point(size = .22, alpha = .78, stroke = 0) + coord_equal() + labs(x = "Standard UMAP 1", y = "Standard UMAP 2", color = "Cluster") + theme_pub() + theme(legend.position = "right")
save_panel(pS3b, file.path(s3, "Supplementary_Figure_S3B_standard_UMAP_by_cluster"), 6.3, 5.0)
ra <- read.csv(file.path(std_dir, "standard_reclustering_resolution_audit.csv"), stringsAsFactors = FALSE)
pS3c <- ggplot(ra, aes(resolution, mean_silhouette)) + geom_line(color = "#4C78A8", linewidth = .8) + geom_point(aes(fill = resolution == .2), shape = 21, size = 2.8, color = "#333333") + scale_fill_manual(values = c(`FALSE` = "white", `TRUE` = "#B84A62"), guide = "none") + scale_x_continuous(breaks = ra$resolution) + labs(x = "Clustering resolution", y = "Mean silhouette width") + theme_pub()
save_panel(pS3c, file.path(s3, "Supplementary_Figure_S3C_standard_resolution_audit"), 5.5, 4.2)

conf_frac <- read.csv(file.path(std_dir, "standard_program_fractions_confident_clusters.csv"), stringsAsFactors = FALSE)
conf_frac$sample <- factor(conf_frac$sample, levels = c("HNS1", "HNS2", "HNS3", "HK1", "HK2", "HK3"))
conf_frac$standard_program <- factor(conf_frac$standard_program, levels = rev(program_order))
pS3d <- ggplot(conf_frac, aes(sample, proportion, fill = standard_program)) +
  geom_col(width = .72, color = "white", linewidth = .2) +
  scale_fill_manual(values = prog_cols[rev(program_order)], name = NULL) +
  scale_y_continuous(labels = percent) +
  labs(x = NULL, y = "Fibroblast state fraction") + theme_pub() +
  theme(axis.text.x = element_text(angle = 20, hjust = 1), legend.position = "right")
save_panel(pS3d, file.path(s3, "Supplementary_Figure_S3D_confident_cluster_state_composition_sensitivity"), 6.5, 4.5)

pres_loo <- read.csv(file.path(std_dir, "standard_module_preservation_leave_one_donor_out.csv"), stringsAsFactors = FALSE) %>%
  filter(module %in% ps$module) %>%
  mutate(program = factor(ifelse(module == ps$module[ps$program == "Clean D1-aligned"], "Clean D1-aligned", "Clean D2-aligned"),
                          levels = c("Clean D1-aligned", "Clean D2-aligned")))
pS3e <- ggplot(pres_loo, aes(omitted_donor, Zsummary.pres, color = program, group = program)) +
  geom_hline(yintercept = 2, linetype = "dashed", color = "#52796F", linewidth = .55) +
  geom_line(linewidth = .7) + geom_point(size = 2.4) +
  facet_wrap(~program, scales = "free_y", ncol = 1) +
  scale_color_manual(values = c("Clean D1-aligned" = "#2F8F83", "Clean D2-aligned" = "#B84A62"), guide = "none") +
  labs(x = "Omitted validation donor", y = "Module preservation Zsummary") + theme_pub() +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))
save_panel(pS3e, file.path(s3, "Supplementary_Figure_S3E_module_preservation_leave_one_donor_out"), 6.6, 5.2)

# Supplementary Figure S4: primary bulk-score comparisons.
s4 <- file.path(supp_root, "Supplementary_Figure_S4"); dir.create(s4, recursive = TRUE, showWarnings = FALSE)
cormat <- cor(scores[, c("clean_D1_loss", "clean_D2", "five_hub")], method = "spearman")
score_display <- c(clean_D1_loss = "Clean D1 loss", clean_D2 = "Clean D2", five_hub = "Five-hub")
rownames(cormat) <- unname(score_display[rownames(cormat)])
colnames(cormat) <- unname(score_display[colnames(cormat)])
cdf <- as.data.frame(cormat) %>% tibble::rownames_to_column("score") %>% pivot_longer(-score, names_to = "score2", values_to = "rho")
pS4a <- ggplot(cdf, aes(score2, score, fill = rho)) + geom_tile(color = "white") + geom_text(aes(label = sprintf("%.2f", rho)), size = 2.5) + scale_fill_gradient2(low = "#3A6EA5", mid = "white", high = "#B84A62", limits = c(-1,1), name = "rho") + theme_minimal(base_size = 8.5) + theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 35, hjust = 1), plot.title = element_blank()) + labs(x = NULL, y = NULL)
save_panel(pS4a, file.path(s4, "Supplementary_Figure_S4A_bulk_score_correlations"), 7.0, 5.8)
pair <- read.csv(file.path(root, "results/bulk_recurrence/bulk_paired_bootstrap_AUC_comparisons.csv"), stringsAsFactors = FALSE)
pair$comparison <- recode(pair$comparison, clean_D1loss_vs_fivehub = "Clean D1 loss - Five-hub", clean_D1loss_vs_clean_D2 = "Clean D1 loss - Clean D2")
pS4b <- ggplot(pair, aes(observed_auc_difference, reorder(comparison, observed_auc_difference))) + geom_errorbarh(aes(xmin = bootstrap_ci_low, xmax = bootstrap_ci_high), height = .18, color = "#666666") + geom_point(size = 2.5, color = "#2F8F83") + geom_vline(xintercept = 0, linetype = "dashed") + labs(x = "Paired bootstrap AUC difference (95% CI)", y = NULL) + theme_pub()
save_panel(pS4b, file.path(s4, "Supplementary_Figure_S4B_paired_AUC_differences"), 7.0, 4.5)

cat("submission figures built\n")
