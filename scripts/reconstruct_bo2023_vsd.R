suppressPackageStartupMessages({
  library(DESeq2)
  library(limma)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) {
  stop(
    "Usage: Rscript reconstruct_bo2023_vsd.R ",
    "<counts.tsv> <author_vsd.tsv> <metadata.csv> <output_dir>"
  )
}

counts_path <- args[[1]]
target_path <- args[[2]]
metadata_path <- args[[3]]
output_dir <- args[[4]]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

message("Reading Bo2023 matrices...")
counts <- as.matrix(read.delim(
  counts_path,
  row.names = 1,
  check.names = FALSE,
  stringsAsFactors = FALSE
))
storage.mode(counts) <- "integer"

target <- as.matrix(read.delim(
  target_path,
  row.names = 1,
  check.names = FALSE,
  stringsAsFactors = FALSE
))
storage.mode(target) <- "double"

metadata <- read.csv(
  metadata_path,
  row.names = 1,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

if (!identical(colnames(counts), colnames(target))) {
  stop("Raw count and author VSD sample columns are not identical.")
}
if (!all(colnames(counts) %in% rownames(metadata))) {
  stop("Metadata is missing one or more count-matrix samples.")
}
if (!all(rownames(target) %in% rownames(counts))) {
  stop("Author VSD genes are not a subset of raw count genes.")
}

metadata <- metadata[colnames(counts), , drop = FALSE]
counts <- counts[rownames(target), , drop = FALSE]
stopifnot(identical(colnames(counts), rownames(metadata)))
stopifnot(identical(rownames(counts), rownames(target)))

metadata$Batch <- factor(metadata$Batch)
metadata$batch2 <- factor(metadata$batch2)
metadata$Region <- factor(metadata$Region)
metadata$SaleemNetworks <- factor(metadata$SaleemNetworks)

safe_cor <- function(x, y) {
  keep <- is.finite(x) & is.finite(y)
  if (sum(keep) < 3L || sd(x[keep]) == 0 || sd(y[keep]) == 0) {
    return(NA_real_)
  }
  cor(x[keep], y[keep])
}

matrix_metrics <- function(candidate, target_matrix, transform_name, correction_name) {
  delta <- candidate - target_matrix
  sample_cor <- vapply(
    seq_len(ncol(candidate)),
    function(i) safe_cor(candidate[, i], target_matrix[, i]),
    numeric(1)
  )
  gene_cor <- vapply(
    seq_len(nrow(candidate)),
    function(i) safe_cor(candidate[i, ], target_matrix[i, ]),
    numeric(1)
  )
  flat_candidate <- as.numeric(candidate)
  flat_target <- as.numeric(target_matrix)
  affine <- lm(flat_target ~ flat_candidate)
  affine_residual <- residuals(affine)

  data.frame(
    transform = transform_name,
    correction = correction_name,
    global_cor = safe_cor(flat_candidate, flat_target),
    median_sample_cor = median(sample_cor, na.rm = TRUE),
    p10_sample_cor = unname(quantile(sample_cor, 0.10, na.rm = TRUE)),
    median_gene_cor = median(gene_cor, na.rm = TRUE),
    p10_gene_cor = unname(quantile(gene_cor, 0.10, na.rm = TRUE)),
    rmse = sqrt(mean(delta^2)),
    mae = mean(abs(delta)),
    affine_rmse = sqrt(mean(affine_residual^2)),
    affine_intercept = unname(coef(affine)[[1]]),
    affine_slope = unname(coef(affine)[[2]]),
    stringsAsFactors = FALSE
  )
}

make_design <- function(column) {
  model.matrix(reformulate(column), data = metadata)
}

correction_specs <- list(
  none = list(batch = NULL, batch2 = NULL, design = NULL),
  platform = list(batch = metadata$batch2, batch2 = NULL, design = NULL),
  sequencing_batch = list(batch = metadata$Batch, batch2 = NULL, design = NULL),
  platform_preserve_region = list(
    batch = metadata$batch2,
    batch2 = NULL,
    design = make_design("Region")
  ),
  sequencing_batch_preserve_region = list(
    batch = metadata$Batch,
    batch2 = NULL,
    design = make_design("Region")
  ),
  platform_preserve_network = list(
    batch = metadata$batch2,
    batch2 = NULL,
    design = make_design("SaleemNetworks")
  ),
  sequencing_batch_preserve_network = list(
    batch = metadata$Batch,
    batch2 = NULL,
    design = make_design("SaleemNetworks")
  )
)

apply_correction <- function(matrix_values, spec) {
  if (is.null(spec$batch) && is.null(spec$batch2)) {
    return(matrix_values)
  }
  correction_args <- list(x = matrix_values)
  if (!is.null(spec$batch)) {
    correction_args$batch <- spec$batch
  }
  if (!is.null(spec$batch2)) {
    correction_args$batch2 <- spec$batch2
  }
  if (!is.null(spec$design)) {
    correction_args$design <- spec$design
  }
  do.call(removeBatchEffect, correction_args)
}

fit_transform <- function(method) {
  dds <- DESeqDataSetFromMatrix(
    countData = counts,
    colData = metadata,
    design = ~1
  )
  dds <- estimateSizeFactors(dds)
  if (method == "exact_vst") {
    dds <- estimateDispersions(dds, fitType = "parametric", quiet = TRUE)
    transformed <- varianceStabilizingTransformation(dds, blind = FALSE)
  } else if (method == "fast_vst") {
    transformed <- vst(dds, blind = TRUE, fitType = "parametric")
    dds <- estimateDispersions(dds, fitType = "parametric", quiet = TRUE)
  } else {
    stop("Unknown transform method: ", method)
  }
  list(
    matrix = assay(transformed),
    dds = dds,
    size_factors = sizeFactors(dds)
  )
}

all_metrics <- list()
best <- NULL
best_rmse <- Inf
best_name <- NULL
best_fit <- NULL

for (transform_name in c("exact_vst", "fast_vst")) {
  message("Fitting ", transform_name, "...")
  fit <- fit_transform(transform_name)

  for (correction_name in names(correction_specs)) {
    message("Evaluating ", transform_name, " + ", correction_name)
    candidate <- apply_correction(
      fit$matrix,
      correction_specs[[correction_name]]
    )
    metric <- matrix_metrics(
      candidate,
      target,
      transform_name,
      correction_name
    )
    all_metrics[[length(all_metrics) + 1L]] <- metric

    if (is.finite(metric$rmse) && metric$rmse < best_rmse) {
      best_rmse <- metric$rmse
      best <- candidate
      best_name <- paste(transform_name, correction_name, sep = "__")
      best_fit <- fit
    }
    rm(candidate)
    gc(verbose = FALSE)
  }
}

metrics <- do.call(rbind, all_metrics)
metrics <- metrics[order(metrics$rmse, -metrics$global_cor), ]
write.csv(
  metrics,
  file.path(output_dir, "candidate_reconstruction_metrics.csv"),
  row.names = FALSE
)

best_output <- gzfile(
  file.path(output_dir, "best_reconstructed_vsd.tsv.gz"),
  open = "wt"
)
write.table(
  best,
  best_output,
  sep = "\t",
  quote = FALSE,
  col.names = NA
)
close(best_output)

ratio_geomeans <- exp(rowMeans(log(counts)))
positive_geomeans <- exp(
  rowSums(log(pmax(counts, 1))) / rowSums(counts > 0)
)
positive_geomeans[rowSums(counts > 0) == 0] <- NA_real_

reference_fit <- list(
  selected_candidate = best_name,
  genes = rownames(counts),
  samples = colnames(counts),
  size_factors = best_fit$size_factors,
  ratio_geomeans = ratio_geomeans,
  positive_geomeans = positive_geomeans,
  dispersion_function = dispersionFunction(best_fit$dds),
  session_info = capture.output(sessionInfo())
)
saveRDS(
  reference_fit,
  file.path(output_dir, "bo2023_frozen_vst_reference.rds")
)

writeLines(
  c(
    paste("selected_candidate:", best_name),
    paste("selected_rmse:", format(best_rmse, digits = 8)),
    paste("genes:", nrow(counts)),
    paste("samples:", ncol(counts))
  ),
  file.path(output_dir, "reconstruction_summary.txt")
)

message("Completed. Best candidate: ", best_name, "; RMSE=", best_rmse)
