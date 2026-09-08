#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(argparse)
  library(ENIGMA)
})

args0 <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args0, value = TRUE)
run_dir <- if (length(file_arg) > 0) dirname(normalizePath(sub("^--file=", "", file_arg[[1]]))) else getwd()
source(file.path(dirname(run_dir), "modules", "config_helpers.R"))
source(file.path(dirname(run_dir), "modules", "runner_helpers.R"))
source(file.path(dirname(run_dir), "modules", "mapping_helpers.R"))
normalize_enigma_mode <- function(mode) {
  if (!mode %in% c("L2", "trace")) {
    stop("ENIGMAmode must be one of: L2, trace. Received: ", mode)
  }
  mode
}

enigma_method_name <- function(mode) {
  if (mode == "L2") {
    return("ENIGMAL2")
  }
  if (mode == "trace") {
    return("ENIGMAtrace")
  }
  stop("Unsupported ENIGMAmode: ", mode)
}

validate_enigma_array <- function(z_array, genes, samples, cell_types, label) {
  if (length(dim(z_array)) != 3) {
    stop(label, " must be a 3D array")
  }
  expected_dim <- c(length(genes), length(samples), length(cell_types))
  if (!identical(dim(z_array), expected_dim)) {
    stop(
      label, " has unexpected dimensions: ", paste(dim(z_array), collapse = " x "),
      "; expected ", paste(expected_dim, collapse = " x ")
    )
  }
  dimnames(z_array) <- list(genes, samples, cell_types)
  z_array
}

extract_enigma_array <- function(egm, norm_output, model_tracker) {
  if (isTRUE(model_tracker)) {
    return(sce2array(egm, norm_output = norm_output, model_name = "log"))
  }
  sce2array(egm, norm_output = norm_output)
}

parser <- ArgumentParser(add_help = TRUE)
parser$add_argument("--dataset", type = "character", required = TRUE)
parser$add_argument("--config_id", type = "character", required = TRUE)
parser$add_argument("--n_core", type = "integer", required = FALSE, default = NULL)
parser$add_argument("--extra_args", type = "character", required = FALSE, default = "")

args <- parser$parse_args()
repo_root <- find_repo_root(run_dir)
extra <- apply_method_default_extra_args(parse_extra_args(args$extra_args), "ENIGMA", repo_root = repo_root)
validate_extra_args(
  extra,
  c("map_cell_types", "use_test_samples", "ENIGMAmode", "model_tracker", "export_raw"),
  "ENIGMA"
)

dataset <- args$dataset
config_id <- args$config_id
n_core <- args$n_core
map_cell_types <- extra_arg(extra, "map_cell_types", TRUE)
use_test_samples <- extra_arg(extra, "use_test_samples", TRUE)
ENIGMAmode <- normalize_enigma_mode(extra_arg(extra, "ENIGMAmode", "L2"))
model_tracker <- extra_arg(extra, "model_tracker", FALSE)
export_raw <- extra_arg(extra, "export_raw", FALSE)
method_name <- enigma_method_name(ENIGMAmode)
raw_method_name <- paste0(method_name, "_raw")

paths <- resolve_deconv_paths(dataset, config_id, repo_root = repo_root)
bulk_input <- paths$config$bulk_input[[1]]
bulk_scale <- paths$config$bulk_scale[[1]]
bulk_normalization <- paths$config$bulk_normalization[[1]]
frac_input <- paths$config$frac_input[[1]]
refType <- paths$refType
message_log_path <- method_message_log_path(paths, method_name)
cleanup_message_log <- start_message_log(message_log_path)
on.exit(cleanup_message_log(), add = TRUE)

bulk_expr <- read_bulk(paths)
bulk_prep <- prepare_bulk_for_deconv(bulk_expr, paths, method = "ENIGMA")
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
bulk_input_dim <- format_deconv_dim(bulk_expr)

validate_non_log_bulk_input(bulk_expr, bulk_input, config_id, method = "ENIGMA")

ENIGMA_ref <- read_reference_matrix(paths$ref_dir, c("rowMeans_sig.csv", "rowMeans_sig.txt"), "ENIGMA reference")
aligned <- align_fraction_and_reference_cell_types(
  frac,
  colnames(ENIGMA_ref),
  paths,
  repo_root
)
frac <- aligned$frac
ENIGMA_ref <- ENIGMA_ref[, aligned$reference_cell_types, drop = FALSE]

genes <- intersect(rownames(bulk_expr), rownames(ENIGMA_ref))
if (length(genes) == 0) {
  stop("No genes available for ENIGMA after intersecting bulk and reference")
}
bulk_expr <- bulk_expr[genes, , drop = FALSE]
ENIGMA_ref <- ENIGMA_ref[genes, , drop = FALSE]
deconv_input_dim <- format_deconv_dim(bulk_expr)

if (!identical(rownames(bulk_expr), rownames(ENIGMA_ref))) {
  stop("ENIGMA bulk and reference gene names are not aligned after subsetting")
}
if (!identical(colnames(ENIGMA_ref), colnames(frac))) {
  stop("ENIGMA reference and fraction cell types are not aligned after subsetting")
}

message("Running ENIGMA")
message("  dataset: ", dataset)
message("  config_id: ", config_id)
message("  bulk_input: ", bulk_input, " (", bulk_input_dim, ")")
message("  bulk_scale: ", bulk_scale)
message("  bulk_normalization: ", bulk_normalization)
message("  bulk_preparation: ", bulk_prep$action)
message("  refType: ", refType)
message("  frac_input: ", frac_input, " (", nrow(frac), " samples x ", ncol(frac), " cell types)")
message("  ENIGMAmode: ", ENIGMAmode)
message("  method: ", method_name)
message("  use_test_samples: ", use_test_samples)
if (!is.null(n_core) && !is.na(n_core)) {
  message("  note: ENIGMA does not use n_core; value is ignored")
}
message("  map_cell_types: ", map_cell_types)
message("  model_tracker: ", model_tracker)
message("  export_raw: ", export_raw)
message("  reference_dim: ", nrow(ENIGMA_ref), " genes x ", ncol(ENIGMA_ref), " cell types")
message("Starting deconvolution")
message("  deconv_input: ", deconv_input_dim)

with_runtime_log(paths, method_name, deconv_input = deconv_input_dim, expr = {
  egm <- create_ENIGMA(bulk = bulk_expr, ref = ENIGMA_ref, ref_type = "aggre")
  egm@result_cell_proportion <- frac

  if (ENIGMAmode == "L2") {
    egm <- ENIGMA_L2_max_norm(
      egm,
      alpha = 0.1,
      model_tracker = model_tracker,
      model_name = "log",
      preprocess = "log"
    )
  } else if (ENIGMAmode == "trace") {
    egm <- ENIGMA_trace_norm(
      egm,
      alpha = 0.1,
      model_tracker = model_tracker,
      model_name = "log",
      preprocess = "log"
    )
  }

  cse_normalized <- extract_enigma_array(egm, norm_output = TRUE, model_tracker = model_tracker)
  cse_normalized <- validate_enigma_array(
    cse_normalized,
    rownames(bulk_expr),
    colnames(bulk_expr),
    colnames(frac),
    paste0(method_name, " normalized output")
  )

  Z_export <- cse_normalized
  if (map_cell_types) {
    Z_export <- apply_indep_ref_cell_type_mapping(Z_export, paths, repo_root = repo_root)
  }

  message_existing_method_output(paths, method_name)
  write_ctse_array(Z_export, paths$output_dir, method_name)

  if (export_raw) {
    cse_raw <- extract_enigma_array(egm, norm_output = FALSE, model_tracker = model_tracker)
    cse_raw <- validate_enigma_array(
      cse_raw,
      rownames(bulk_expr),
      colnames(bulk_expr),
      colnames(frac),
      paste0(method_name, " raw output")
    )
    if (map_cell_types) {
      cse_raw <- apply_indep_ref_cell_type_mapping(cse_raw, paths, repo_root = repo_root)
    }
    message_existing_method_output(paths, raw_method_name)
    write_ctse_array(cse_raw, paths$output_dir, raw_method_name)
  }
})

message("Output saved")
message("  ctse: ", file.path(paths$output_dir, method_name, "<cell_type>.txt.gz"))
if (export_raw) {
  message("  raw ctse: ", file.path(paths$output_dir, raw_method_name, "<cell_type>.txt.gz"))
}
if (map_cell_types && refType == "indep") {
  message("  ctse naming note: mapped cell types use benchmark names; unmapped reference cell types keep reference names")
} else {
  message("  ctse naming note: exported cell types use reference cell type names")
}
message("done")
