#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(argparse)
  library(Unico)
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
extra <- apply_method_default_extra_args(parse_extra_args(args$extra_args), "Unico", repo_root = repo_root)
validate_extra_args(
  extra,
  c("map_cell_types", "use_test_samples", "use_limma_top_genes", "top_n"),
  "Unico"
)

dataset <- args$dataset
config_id <- args$config_id
n_core <- args$n_core
map_cell_types <- extra_arg(extra, "map_cell_types", TRUE)
use_test_samples <- extra_arg(extra, "use_test_samples", TRUE)
use_limma_top_genes <- extra_arg(extra, "use_limma_top_genes", FALSE)
top_n <- extra_arg(extra, "top_n", 100)
validate_top_n_arg(extra, use_limma_top_genes, "Unico")

paths <- resolve_deconv_paths(dataset, config_id, repo_root = repo_root)
bulk_input <- paths$config$bulk_input[[1]]
bulk_scale <- paths$config$bulk_scale[[1]]
bulk_normalization <- paths$config$bulk_normalization[[1]]
frac_input <- paths$config$frac_input[[1]]
refType <- paths$refType
message_log_path <- method_message_log_path(paths, "Unico")
cleanup_message_log <- start_message_log(message_log_path)
on.exit(cleanup_message_log(), add = TRUE)

bulk_expr <- read_bulk(paths)
bulk_prep <- prepare_bulk_for_deconv(bulk_expr, paths, method = "Unico")
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
validate_fraction_cell_types(frac, "Unico")
bulk_input_dim <- format_deconv_dim(bulk_expr)

# X must not include features with standard deviation less than 1e-04.
gene_sd <- matrixStats::rowSds(bulk_expr)
keep_genes <- is.finite(gene_sd) & gene_sd > 1e-4
bulk_expr <- bulk_expr[keep_genes, , drop = FALSE]
if (nrow(bulk_expr) == 0) {
  stop("No genes available for Unico after filtering near-constant genes")
}

limma_gene_count <- NA_integer_
if (use_limma_top_genes) {
  limma_genes <- read_limma_top_gene_union(paths, top_n = top_n)
  limma_gene_count <- length(limma_genes)
  genes <- intersect(rownames(bulk_expr), limma_genes)
  if (length(genes) == 0) {
    stop("No genes available for Unico after applying limma top gene filter")
  }
  bulk_expr <- bulk_expr[genes, , drop = FALSE]
}

deconv_input_dim <- format_deconv_dim(bulk_expr)

message("Running Unico")
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

message("Starting deconvolution")
message("  deconv_input: ", deconv_input_dim)

with_runtime_log(paths, "Unico", n_core = n_core, deconv_input = deconv_input_dim, expr = {
  params_hat <- Unico(bulk_expr, frac, C1 = NULL, C2 = NULL, num_cores = n_core)
  Z_hat <- Unico::tensor(bulk_expr, frac, C1 = NULL, C2 = NULL, params_hat, num_cores = n_core)
  Z_hat <- aperm(Z_hat, c(2, 3, 1))

  if (!identical(dim(Z_hat), c(nrow(bulk_expr), ncol(bulk_expr), ncol(frac)))) {
    stop(
      "Unexpected Unico CTSE dimensions: ", paste(dim(Z_hat), collapse = " x "),
      "; expected ", nrow(bulk_expr), " x ", ncol(bulk_expr), " x ", ncol(frac)
    )
  }
  z_cell_types <- dimnames(Z_hat)[[3]]
  if (!is.null(z_cell_types) && !identical(z_cell_types, colnames(frac))) {
    stop(
      "Unico CTSE cell types do not match fraction input cell types. ",
      "CTSE: ", paste(z_cell_types, collapse = ", "),
      "; frac: ", paste(colnames(frac), collapse = ", ")
    )
  }
  dimnames(Z_hat) <- list(rownames(bulk_expr), colnames(bulk_expr), colnames(frac))

  Z_export <- Z_hat
  if (map_cell_types) {
    Z_export <- apply_indep_ref_cell_type_mapping(Z_export, paths, repo_root = repo_root)
  }

  message_existing_method_output(paths, "Unico")
  write_ctse_array(Z_export, paths$output_dir, "Unico")
})

message("Output saved")
message("  ctse: ", file.path(paths$output_dir, "Unico", "<cell_type>.txt.gz"))
if (map_cell_types && refType == "indep") {
  message("  ctse naming note: mapped cell types use benchmark names; unmapped reference cell types keep reference names")
} else {
  message("  ctse naming note: exported cell types use reference cell type names")
}
message("done")
