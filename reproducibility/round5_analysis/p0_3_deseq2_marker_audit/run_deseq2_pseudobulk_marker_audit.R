#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) stop("usage: run_deseq2_pseudobulk_marker_audit.R COUNTS INPUT_DIR OUTPUT_DIR")
counts_path <- normalizePath(args[[1]], mustWork = TRUE)
input_dir <- normalizePath(args[[2]], mustWork = TRUE)
output_dir <- normalizePath(args[[3]], mustWork = FALSE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
suppressPackageStartupMessages(library(DESeq2))

metadata <- read.delim(file.path(input_dir, "bo2023_819_deseq2_metadata.tsv"), check.names = FALSE, stringsAsFactors = FALSE)
rownames(metadata) <- metadata$sample_id
metadata$donor_id <- factor(metadata$donor_id)
metadata$network <- factor(metadata$network)
raw <- read.delim(counts_path, header = TRUE, row.names = 1, check.names = FALSE)
if (!identical(colnames(raw), rownames(metadata))) stop("count columns do not match metadata")
raw_matrix <- as.matrix(raw)
storage.mode(raw_matrix) <- "integer"

unit <- interaction(metadata$donor_id, metadata$network, drop = TRUE, sep = "__")
unit_levels <- levels(unit)
pb <- sapply(unit_levels, function(u) rowSums(raw_matrix[, unit == u, drop = FALSE]))
pb_meta <- do.call(rbind, lapply(unit_levels, function(u) {
  idx <- which(unit == u)
  data.frame(
    pseudobulk_id = u,
    donor_id = as.character(metadata$donor_id[idx[[1]]]),
    network = as.character(metadata$network[idx[[1]]]),
    n_tissue_samples = length(idx),
    stringsAsFactors = FALSE
  )
}))
rownames(pb_meta) <- pb_meta$pseudobulk_id
pb_meta$donor_id <- factor(pb_meta$donor_id)
pb_meta$network <- factor(pb_meta$network)
if (!identical(colnames(pb), rownames(pb_meta))) stop("pseudobulk metadata order mismatch")
write.csv(pb_meta, file.path(output_dir, "pseudobulk_unit_metadata.csv"), row.names = FALSE)

keep <- rowSums(pb >= 10L) >= 3L
write.csv(data.frame(
  n_raw_genes = nrow(pb), n_tested_genes = sum(keep), n_prefiltered_genes = sum(!keep),
  rule = "pseudobulk raw count >=10 in at least 3 donor-Network units"
), file.path(output_dir, "pseudobulk_gene_filter_summary.csv"), row.names = FALSE)
pb <- pb[keep, , drop = FALSE]

full_matrix <- model.matrix(~ donor_id + network, data = pb_meta)
reduced_matrix <- model.matrix(~ donor_id, data = pb_meta)
if (qr(full_matrix)$rank != ncol(full_matrix)) stop("pseudobulk full design is rank deficient")
write.csv(data.frame(
  n_pseudobulk_units = nrow(pb_meta), n_donors = nlevels(pb_meta$donor_id),
  n_networks = nlevels(pb_meta$network), full_columns = ncol(full_matrix), full_rank = qr(full_matrix)$rank,
  reduced_columns = ncol(reduced_matrix), reduced_rank = qr(reduced_matrix)$rank,
  min_tissues_per_unit = min(pb_meta$n_tissue_samples), max_tissues_per_unit = max(pb_meta$n_tissue_samples)
), file.path(output_dir, "pseudobulk_design_audit.csv"), row.names = FALSE)

dds <- DESeqDataSetFromMatrix(pb, pb_meta, design = ~ donor_id + network)
set.seed(20260803)
dds <- DESeq(dds, test = "LRT", reduced = ~ donor_id, minReplicatesForReplace = Inf, parallel = FALSE, quiet = TRUE)
res <- as.data.frame(results(dds, independentFiltering = TRUE))
res$gene_id <- rownames(res)
annotation <- read.delim(file.path(input_dir, "bo2023_gene_annotation.tsv"), check.names = FALSE, stringsAsFactors = FALSE)
res <- merge(res, annotation, by = "gene_id", all.x = TRUE, sort = FALSE)
res <- res[, c("gene_id", "gene_symbol", "gene_type", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")]
res <- res[order(res$padj, res$pvalue, na.last = TRUE), ]
write.csv(res, file.path(output_dir, "pseudobulk_deseq2_network_lrt_all_genes.csv"), row.names = FALSE, na = "")

norm_log <- log2(counts(dds, normalized = TRUE) + 1)
network_levels <- levels(pb_meta$network)
network_means <- sapply(network_levels, function(net) rowMeans(norm_log[, pb_meta$network == net, drop = FALSE]))
target_idx <- max.col(network_means, ties.method = "first")
target_effect <- network_means[cbind(seq_len(nrow(network_means)), target_idx)] -
  (rowSums(network_means) - network_means[cbind(seq_len(nrow(network_means)), target_idx)]) / (length(network_levels) - 1)
effect <- data.frame(gene_id = rownames(norm_log), highest_network = network_levels[target_idx], pseudobulk_log2_effect_vs_other_network_mean = target_effect)
effect <- merge(effect, annotation, by = "gene_id", all.x = TRUE, sort = FALSE)
write.csv(effect, file.path(output_dir, "pseudobulk_network_effects_all_genes.csv"), row.names = FALSE, na = "")

locked <- read.delim(file.path(input_dir, "locked_network_top200.tsv"), check.names = FALSE, stringsAsFactors = FALSE)
symbol_res <- res[!is.na(res$gene_symbol) & nzchar(res$gene_symbol), ]
symbol_res <- symbol_res[order(symbol_res$padj, symbol_res$pvalue, na.last = TRUE), ]
symbol_res <- symbol_res[!duplicated(symbol_res$gene_symbol), ]
symbol_effect <- effect[!is.na(effect$gene_symbol) & nzchar(effect$gene_symbol), ]
symbol_effect <- symbol_effect[!duplicated(symbol_effect$gene_symbol), ]
audit <- merge(locked, symbol_res, by = "gene_symbol", all.x = TRUE, sort = FALSE)
audit <- merge(audit, symbol_effect[, c("gene_symbol", "highest_network", "pseudobulk_log2_effect_vs_other_network_mean")], by = "gene_symbol", all.x = TRUE, sort = FALSE)
audit$lrt_tested <- !is.na(audit$pvalue)
audit$lrt_fdr_lt_0_05 <- !is.na(audit$padj) & audit$padj < 0.05
audit <- audit[order(audit$fisher_score, decreasing = TRUE), ]
write.csv(audit, file.path(output_dir, "pseudobulk_locked_top200_deseq2_support.csv"), row.names = FALSE, na = "")

tested_locked <- audit[audit$lrt_tested, ]
summary <- data.frame(
  metric = c("raw_genes", "tested_genes", "pseudobulk_units", "lrt_fdr_lt_0_05_all_genes", "locked_panel_genes",
             "locked_panel_mapped_and_tested", "locked_panel_lrt_fdr_lt_0_05", "locked_panel_support_fraction",
             "spearman_fisher_vs_lrt_stat", "spearman_fisher_vs_neglog10_padj"),
  value = c(nrow(raw), nrow(res), nrow(pb_meta), sum(res$padj < 0.05, na.rm = TRUE), nrow(locked), nrow(tested_locked),
            sum(tested_locked$lrt_fdr_lt_0_05), mean(tested_locked$lrt_fdr_lt_0_05),
            suppressWarnings(cor(tested_locked$fisher_score, tested_locked$stat, method = "spearman", use = "complete.obs")),
            suppressWarnings(cor(tested_locked$fisher_score, -log10(tested_locked$padj), method = "spearman", use = "complete.obs")))
)
write.csv(summary, file.path(output_dir, "pseudobulk_deseq2_marker_audit_summary.csv"), row.names = FALSE)

top_var <- order(matrixStats::rowVars(norm_log), decreasing = TRUE)[seq_len(min(500, nrow(norm_log)))]
pca_fit <- prcomp(t(norm_log[top_var, , drop = FALSE]), center = TRUE, scale. = FALSE)
pca <- data.frame(pseudobulk_id = rownames(pca_fit$x), PC1 = pca_fit$x[, 1], PC2 = pca_fit$x[, 2], pb_meta[rownames(pca_fit$x), c("network", "donor_id", "n_tissue_samples")])
write.csv(pca, file.path(output_dir, "pseudobulk_pca_coordinates.csv"), row.names = FALSE)
percent_var <- pca_fit$sdev^2 / sum(pca_fit$sdev^2)
png(file.path(output_dir, "pseudobulk_pca_network.png"), width = 1800, height = 1400, res = 180)
cols <- grDevices::hcl.colors(length(network_levels), "Dark 3")
plot(pca$PC1, pca$PC2, col = cols[as.integer(factor(pca$network, levels = network_levels))], pch = 16,
     xlab = paste0("PC1: ", round(percent_var[[1]] * 100, 1), "%"), ylab = paste0("PC2: ", round(percent_var[[2]] * 100, 1), "%"),
     main = "Bo2023 donor-Network pseudobulk PCA")
legend("topright", legend = network_levels, col = cols, pch = 16, cex = 0.65)
dev.off()
png(file.path(output_dir, "pseudobulk_deseq2_dispersion_plot.png"), width = 1600, height = 1200, res = 180)
plotDispEsts(dds)
dev.off()
writeLines(capture.output(sessionInfo()), file.path(output_dir, "pseudobulk_session_info.txt"))
writeLines(c(
  paste("counts source file:", basename(counts_path)), paste("raw tissue samples:", ncol(raw_matrix)), paste("donor-Network pseudobulk units:", ncol(pb)),
  paste("raw genes:", nrow(raw)), paste("tested genes:", nrow(pb)), "aggregation: raw counts summed within each observed donor-Network cell",
  "full design: ~ donor_id + network", "reduced design: ~ donor_id", "primary test: DESeq2 LRT for any Network effect",
  "multiplicity: BH across tested genes", "production marker panel unchanged: inferential audit only"
), file.path(output_dir, "pseudobulk_run_manifest.txt"))
print(summary)
