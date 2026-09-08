#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(argparse)
  library(TCA)
})

args0 <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args0, value = TRUE)
run_dir <- if (length(file_arg) > 0) dirname(normalizePath(sub("^--file=", "", file_arg[[1]]))) else getwd()
source(file.path(dirname(run_dir), "modules", "config_helpers.R"))
source(file.path(dirname(run_dir), "modules", "runner_helpers.R"))
source(file.path(dirname(run_dir), "modules", "marker_helpers.R"))
source(file.path(dirname(run_dir), "modules", "mapping_helpers.R"))
parser <- ArgumentParser(add_help = TRUE)
parser$add_argument("--dataset", type = "character", required = TRUE)
parser$add_argument("--config_id", type = "character", required = TRUE)
parser$add_argument("--n_core", type = "integer", required = FALSE, default = 15)
parser$add_argument("--extra_args", type = "character", required = FALSE, default = "")

args <- parser$parse_args()
repo_root <- find_repo_root(run_dir)
extra <- apply_method_default_extra_args(parse_extra_args(args$extra_args), "TCA", repo_root = repo_root)
validate_extra_args(
  extra,
  c("map_cell_types", "use_test_samples", "use_limma_top_genes", "top_n"),
  "TCA"
)

dataset <- args$dataset
config_id <- args$config_id
n_core <- args$n_core
map_cell_types <- extra_arg(extra, "map_cell_types", TRUE)
use_test_samples <- extra_arg(extra, "use_test_samples", TRUE)
use_limma_top_genes <- extra_arg(extra, "use_limma_top_genes", FALSE)
top_n <- extra_arg(extra, "top_n", 100)
validate_top_n_arg(extra, use_limma_top_genes, "TCA")

paths <- resolve_deconv_paths(dataset, config_id, repo_root = repo_root)
bulk_input <- paths$config$bulk_input[[1]]
bulk_scale <- paths$config$bulk_scale[[1]]
bulk_normalization <- paths$config$bulk_normalization[[1]]
frac_input <- paths$config$frac_input[[1]]
refType <- paths$refType
message_log_path <- method_message_log_path(paths, "TCA")
cleanup_message_log <- start_message_log(message_log_path)
on.exit(cleanup_message_log(), add = TRUE)

bulk_expr <- read_bulk(paths)
bulk_prep <- prepare_bulk_for_deconv(bulk_expr, paths, method = "TCA")
bulk_expr <- bulk_prep$bulk_expr
frac <- read_frac_for_method(paths, repo_root = repo_root)
sample_restriction <- restrict_to_test_samples(
  bulk_expr,
  paths,
  frac = frac,
  use_test_samples = use_test_samples
)
bulk_expr <- sample_restriction$bulk
frac <- sample_restriction$frac
frac <- validate_bulk_frac(bulk_expr, frac)
validate_fraction_cell_types(frac, "TCA")
bulk_input_dim <- format_deconv_dim(bulk_expr)

log_bulk <- prepare_log_bulk_input(bulk_expr, bulk_input, config_id, method = "TCA")
X <- log_bulk$bulk_expr

# Preserve archived TCA behavior: exclude genes with logged variance <= 1e-8.
gene_var <- matrixStats::rowVars(X)
keep_genes <- is.finite(gene_var) & gene_var > 1e-8
X <- X[keep_genes, , drop = FALSE]
if (nrow(X) == 0) {
  stop("No genes available for TCA after filtering near-constant genes")
}

limma_gene_count <- NA_integer_
if (use_limma_top_genes) {
  limma_genes <- read_limma_top_gene_union(paths, top_n = top_n)
  limma_gene_count <- length(limma_genes)
  genes <- intersect(rownames(X), limma_genes)
  if (length(genes) == 0) {
    stop("No genes available for TCA after applying limma top gene filter")
  }
  X <- X[genes, , drop = FALSE]
}

deconv_input_dim <- format_deconv_dim(X)

message("Running TCA")
message("  dataset: ", dataset)
message("  config_id: ", config_id)
message("  bulk_input: ", bulk_input, " (", bulk_input_dim, ")")
message("  bulk_scale: ", bulk_scale)
message("  bulk_normalization: ", bulk_normalization)
message("  bulk_preparation: ", bulk_prep$action)
message("  refType: ", refType)
message("  frac_input: ", frac_input, " (", nrow(frac), " samples x ", ncol(frac), " cell types)")
message("  use_test_samples: ", use_test_samples)
message("  n_core: ", n_core)
message("  map_cell_types: ", map_cell_types)
message("  use_limma_top_genes: ", use_limma_top_genes)
if (use_limma_top_genes) {
  message("  top_n: ", top_n)
  message("  limma_top_gene_union: ", limma_gene_count, " genes before intersecting bulk")
}
message("  note: ", log_bulk$note)
message("Starting deconvolution")
message("  deconv_input: ", deconv_input_dim)

with_runtime_log(paths, "TCA", n_core = n_core, deconv_input = deconv_input_dim, expr = {
  TCA_res <- tca(X = X, W = frac, parallel = TRUE, num_cores = n_core)
  Z_hat <- TCA::tensor(X = X, tca.mdl = TCA_res)

  if (!isTRUE(all.equal(colnames(TCA_res[["mus_hat"]]), colnames(frac)))) {
    stop("TCA model cell types do not match fraction input cell types")
  }

  Z_array <- array(NA_real_, dim = c(nrow(Z_hat[[1]]), ncol(Z_hat[[1]]), length(Z_hat)))
  dimnames(Z_array) <- list(rownames(Z_hat[[1]]), colnames(Z_hat[[1]]), colnames(frac))
  for (i in seq_along(Z_hat)) {
    Z_array[, , i] <- Z_hat[[i]]
  }

  if (!identical(dim(Z_array), c(nrow(X), ncol(X), ncol(frac)))) {
    stop(
      "Unexpected TCA CTSE dimensions: ", paste(dim(Z_array), collapse = " x "),
      "; expected ", nrow(X), " x ", ncol(X), " x ", ncol(frac)
    )
  }

  if (file.exists("TCA.log")) {
    unlink("TCA.log")
  }

  Z_export <- Z_array
  if (map_cell_types) {
    Z_export <- apply_indep_ref_cell_type_mapping(Z_export, paths, repo_root = repo_root)
  }

  message_existing_method_output(paths, "TCA")
  write_ctse_array(Z_export, paths$output_dir, "TCA")
})

message("Output saved")
message("  ctse: ", file.path(paths$output_dir, "TCA", "<cell_type>.txt.gz"))
if (map_cell_types && refType == "indep") {
  message("  ctse naming note: mapped cell types use benchmark names; unmapped reference cell types keep reference names")
} else {
  message("  ctse naming note: exported cell types use reference cell type names")
}
message("done")
