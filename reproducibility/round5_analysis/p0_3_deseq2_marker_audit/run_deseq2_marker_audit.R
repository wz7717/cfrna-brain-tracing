#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) stop("usage: run_deseq2_marker_audit.R COUNTS INPUT_DIR OUTPUT_DIR")
counts_path <- normalizePath(args[[1]], mustWork = TRUE)
input_dir <- normalizePath(args[[2]], mustWork = TRUE)
output_dir <- normalizePath(args[[3]], mustWork = FALSE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages(library(DESeq2))

metadata <- read.delim(file.path(input_dir, "bo2023_819_deseq2_metadata.tsv"), check.names = FALSE, stringsAsFactors = FALSE)
rownames(metadata) <- metadata$sample_id
metadata$donor_id <- factor(metadata$donor_id)
metadata$network <- factor(metadata$network)

counts <- read.delim(counts_path, header = TRUE, row.names = 1, check.names = FALSE)
if (!identical(colnames(counts), rownames(metadata))) stop("count columns do not exactly match ordered metadata")
count_matrix <- as.matrix(counts)
storage.mode(count_matrix) <- "integer"
if (anyNA(count_matrix) || any(count_matrix < 0)) stop("counts contain missing or negative values")

keep <- rowSums(count_matrix >= 10L) >= 9L
filter_table <- data.frame(
  n_raw_genes = nrow(count_matrix),
  n_tested_genes = sum(keep),
  n_prefiltered_genes = sum(!keep),
  rule = "raw count >=10 in at least 9 samples"
)
write.csv(filter_table, file.path(output_dir, "gene_filter_summary.csv"), row.names = FALSE)
count_matrix <- count_matrix[keep, , drop = FALSE]

dds <- DESeqDataSetFromMatrix(countData = count_matrix, colData = metadata, design = ~ donor_id + network)
full_matrix <- model.matrix(~ donor_id + network, data = as.data.frame(colData(dds)))
reduced_matrix <- model.matrix(~ donor_id, data = as.data.frame(colData(dds)))
if (qr(full_matrix)$rank != ncol(full_matrix)) stop("full DESeq2 design matrix is rank deficient")

set.seed(20260803)
dds_lrt <- DESeq(dds, test = "LRT", reduced = ~ donor_id, minReplicatesForReplace = Inf, parallel = FALSE, quiet = TRUE)
lrt <- as.data.frame(results(dds_lrt, independentFiltering = TRUE))
lrt$gene_id <- rownames(lrt)
lrt <- lrt[, c("gene_id", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")]

annotation <- read.delim(file.path(input_dir, "bo2023_gene_annotation.tsv"), check.names = FALSE, stringsAsFactors = FALSE)
lrt <- merge(lrt, annotation, by = "gene_id", all.x = TRUE, sort = FALSE)
lrt <- lrt[order(lrt$padj, lrt$pvalue, na.last = TRUE), ]
write.csv(lrt, file.path(output_dir, "deseq2_network_lrt_all_genes.csv"), row.names = FALSE, na = "")

network_levels <- levels(metadata$network)

# The global LRT is the sole formal primary test. Effect direction is summarized
# descriptively from DESeq2 size-factor-normalized counts after donor centering,
# avoiding ten extra families of post-hoc gene-wise tests.
log_norm <- log2(counts(dds_lrt, normalized = TRUE) + 1)
donor_centered <- log_norm
for (donor in levels(metadata$donor_id)) {
  idx <- which(metadata$donor_id == donor)
  donor_centered[, idx] <- sweep(log_norm[, idx, drop = FALSE], 1, rowMeans(log_norm[, idx, drop = FALSE]), "-")
}
network_means <- sapply(network_levels, function(net) rowMeans(donor_centered[, metadata$network == net, drop = FALSE]))
target_idx <- max.col(network_means, ties.method = "first")
other_mean <- (rowSums(network_means) - network_means[cbind(seq_len(nrow(network_means)), target_idx)]) / (length(network_levels) - 1)
effect <- data.frame(
  gene_id = rownames(network_means),
  highest_network = network_levels[target_idx],
  donor_centered_log2_effect_vs_other_network_mean = network_means[cbind(seq_len(nrow(network_means)), target_idx)] - other_mean,
  stringsAsFactors = FALSE
)
effect <- merge(effect, annotation, by = "gene_id", all.x = TRUE, sort = FALSE)
write.csv(effect, file.path(output_dir, "donor_centered_network_effects_all_genes.csv"), row.names = FALSE, na = "")

locked <- read.delim(file.path(input_dir, "locked_network_top200.tsv"), check.names = FALSE, stringsAsFactors = FALSE)
symbol_lrt <- lrt[!is.na(lrt$gene_symbol) & nzchar(lrt$gene_symbol), ]
symbol_lrt <- symbol_lrt[order(symbol_lrt$padj, symbol_lrt$pvalue, na.last = TRUE), ]
symbol_lrt <- symbol_lrt[!duplicated(symbol_lrt$gene_symbol), ]
audit <- merge(locked, symbol_lrt, by = "gene_symbol", all.x = TRUE, sort = FALSE)
symbol_effect <- effect[!is.na(effect$gene_symbol) & nzchar(effect$gene_symbol), ]
symbol_effect <- symbol_effect[!duplicated(symbol_effect$gene_symbol), ]
audit <- merge(audit, symbol_effect[, c("gene_symbol", "highest_network", "donor_centered_log2_effect_vs_other_network_mean")], by = "gene_symbol", all.x = TRUE, sort = FALSE)
audit$lrt_tested <- !is.na(audit$pvalue)
audit$lrt_fdr_lt_0_05 <- !is.na(audit$padj) & audit$padj < 0.05

# Conditional donor-level intervals for the descriptively selected highest
# Network. These intervals are supporting uncertainty summaries, not a second
# multiplicity-controlled inferential family.
audit$effect_n_donors <- NA_integer_
audit$effect_se <- NA_real_
audit$effect_ci95_low <- NA_real_
audit$effect_ci95_high <- NA_real_
gene_rows <- match(audit$gene_id, rownames(log_norm))
for (i in seq_len(nrow(audit))) {
  row_i <- gene_rows[[i]]
  target <- audit$highest_network[[i]]
  if (is.na(row_i) || is.na(target) || !nzchar(target)) next
  donor_effects <- c()
  for (donor in levels(metadata$donor_id)) {
    target_samples <- metadata$donor_id == donor & metadata$network == target
    other_samples <- metadata$donor_id == donor & metadata$network != target
    if (any(target_samples) && any(other_samples)) {
      donor_effects <- c(donor_effects, mean(log_norm[row_i, target_samples]) - mean(log_norm[row_i, other_samples]))
    }
  }
  n_eff <- length(donor_effects)
  audit$effect_n_donors[[i]] <- n_eff
  if (n_eff >= 3) {
    se <- sd(donor_effects) / sqrt(n_eff)
    crit <- qt(0.975, df = n_eff - 1)
    audit$effect_se[[i]] <- se
    audit$effect_ci95_low[[i]] <- mean(donor_effects) - crit * se
    audit$effect_ci95_high[[i]] <- mean(donor_effects) + crit * se
  }
}
audit <- audit[order(audit$fisher_score, decreasing = TRUE), ]
write.csv(audit, file.path(output_dir, "locked_top200_deseq2_support.csv"), row.names = FALSE, na = "")

tested_locked <- audit[audit$lrt_tested, ]
summary <- data.frame(
  metric = c(
    "raw_genes", "tested_genes", "lrt_fdr_lt_0_05_all_genes", "locked_panel_genes",
    "locked_panel_mapped_and_tested", "locked_panel_lrt_fdr_lt_0_05", "locked_panel_support_fraction",
    "spearman_fisher_vs_lrt_stat", "spearman_fisher_vs_neglog10_padj"
  ),
  value = c(
    nrow(counts), nrow(lrt), sum(lrt$padj < 0.05, na.rm = TRUE), nrow(locked), nrow(tested_locked),
    sum(tested_locked$lrt_fdr_lt_0_05), mean(tested_locked$lrt_fdr_lt_0_05),
    suppressWarnings(cor(tested_locked$fisher_score, tested_locked$stat, method = "spearman", use = "complete.obs")),
    suppressWarnings(cor(tested_locked$fisher_score, -log10(tested_locked$padj), method = "spearman", use = "complete.obs"))
  )
)
write.csv(summary, file.path(output_dir, "deseq2_marker_audit_summary.csv"), row.names = FALSE)

top_var <- order(matrixStats::rowVars(log_norm), decreasing = TRUE)[seq_len(min(500, nrow(log_norm)))]
pca_fit <- prcomp(t(log_norm[top_var, , drop = FALSE]), center = TRUE, scale. = FALSE)
pca <- data.frame(sample_id = rownames(pca_fit$x), PC1 = pca_fit$x[, 1], PC2 = pca_fit$x[, 2], metadata[rownames(pca_fit$x), c("network", "donor_id")])
percent_var <- (pca_fit$sdev^2) / sum(pca_fit$sdev^2)
write.csv(pca, file.path(output_dir, "log2_normalized_pca_coordinates.csv"), row.names = FALSE)
png(file.path(output_dir, "log2_normalized_pca_network.png"), width = 1800, height = 1400, res = 180)
cols <- grDevices::hcl.colors(length(network_levels), "Dark 3")
plot(pca$PC1, pca$PC2, col = cols[as.integer(factor(pca$network, levels = network_levels))], pch = 16,
     xlab = paste0("PC1: ", round(percent_var[[1]] * 100, 1), "%"),
     ylab = paste0("PC2: ", round(percent_var[[2]] * 100, 1), "%"), main = "Bo2023 log2 normalized-count PCA by Network")
legend("topright", legend = network_levels, col = cols, pch = 16, cex = 0.65)
dev.off()

png(file.path(output_dir, "deseq2_dispersion_plot.png"), width = 1600, height = 1200, res = 180)
plotDispEsts(dds_lrt)
dev.off()

writeLines(capture.output(sessionInfo()), file.path(output_dir, "session_info.txt"))
writeLines(c(
  paste("counts source file:", basename(counts_path)),
  paste("samples:", ncol(count_matrix)),
  paste("raw genes:", nrow(counts)),
  paste("tested genes:", nrow(count_matrix)),
  "full design: ~ donor_id + network",
  "reduced design: ~ donor_id",
  "primary test: likelihood-ratio test for any Network effect",
  "effect summary: donor-centered log2 normalized-count highest-Network versus other-Network mean",
  "uncertainty: conditional donor-level t intervals for locked Top200 effects; descriptive, not a second test family",
  "production marker panel unchanged: this is an inferential audit"
), file.path(output_dir, "run_manifest.txt"))
print(summary)
