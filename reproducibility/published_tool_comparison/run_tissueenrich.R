args <- commandArgs(trailingOnly = TRUE)
repo <- if (length(args) >= 1) args[[1]] else "."
out_dir <- if (length(args) >= 2) args[[2]] else file.path(repo, "reproducibility", "published_tool_comparison", "outputs")
lib <- normalizePath(file.path(dirname(repo), "external_tools", "R-library"), winslash = "/", mustWork = TRUE)
.libPaths(c(lib, .libPaths()))

suppressPackageStartupMessages(library(TissueEnrich))
suppressPackageStartupMessages(library(GSEABase))
suppressPackageStartupMessages(library(SummarizedExperiment))

panel_file <- file.path(repo, "data", "models", "bo2023_saleem_network_top200_model_genes.csv")
background_file <- file.path(dirname(repo), "code", "reproduction_validation_workspace_20260802", "reproducibility", "p0_bio3_projector", "projector_gene_parameters.csv")
stopifnot(file.exists(panel_file), file.exists(background_file))

panel <- unique(na.omit(read.csv(panel_file, check.names = FALSE)$gene_symbol))
background <- unique(na.omit(read.csv(background_file, check.names = FALSE)$gene_symbol))
stopifnot(length(panel) == 200L, all(panel %in% background))

panel_set <- GeneSet(geneIds = panel, organism = "Homo Sapiens", geneIdType = SymbolIdentifier())
background_set <- GeneSet(geneIds = background, organism = "Homo Sapiens", geneIdType = SymbolIdentifier())
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

extract_result <- function(x, dataset_name) {
  se <- x[[1]]
  tab <- setNames(data.frame(assay(se), check.names = FALSE), colData(se)[, 1])
  tab$tissue <- rowData(se)[, 1]
  tab$dataset <- dataset_name
  tab$adjusted_p <- 10^(-tab$Log10PValue)
  tab <- tab[, c("dataset", "tissue", "Log10PValue", "adjusted_p", "Tissue.Specific.Genes", "fold.change", "samples")]
  tab[order(tab$adjusted_p, -tab$fold.change, tab$tissue), ]
}

all_results <- list()
all_unmapped <- list()
for (dataset_id in c(1L, 2L)) {
  dataset_name <- if (dataset_id == 1L) "HPA" else "GTEx"
  result <- teEnrichment(
    inputGenes = panel_set,
    rnaSeqDataset = dataset_id,
    tissueSpecificGeneType = 1L,
    multiHypoCorrection = TRUE,
    backgroundGenes = background_set
  )
  all_results[[dataset_name]] <- extract_result(result, dataset_name)
  all_unmapped[[dataset_name]] <- data.frame(dataset = dataset_name, gene_symbol = geneIds(result[[4]]))
}

write.csv(do.call(rbind, all_results), file.path(out_dir, "tissueenrich_hpa_gtex_results.csv"), row.names = FALSE)
write.csv(do.call(rbind, all_unmapped), file.path(out_dir, "tissueenrich_unmapped_panel_genes.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(out_dir, "tissueenrich_session_info.txt"))

manifest <- c(
  paste0("TissueEnrich=", as.character(packageVersion("TissueEnrich"))),
  paste0("R=", R.version.string),
  paste0("panel_n=", length(panel)),
  paste0("background_n=", length(background)),
  "datasets=HPA,GTEx",
  "tissueSpecificGeneType=All",
  "multiple_testing=BH within each dataset",
  "interpretation=tissue enrichment only; not anatomical localization accuracy"
)
writeLines(manifest, file.path(out_dir, "tissueenrich_run_manifest.txt"))
