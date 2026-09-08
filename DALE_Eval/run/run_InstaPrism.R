#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(argparse)
  library(InstaPrism)
})

args0 <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args0, value = TRUE)
run_dir <- if (length(file_arg) > 0) dirname(normalizePath(sub("^--file=", "", file_arg[[1]]))) else getwd()
source(file.path(dirname(run_dir), "modules", "config_helpers.R"))
source(file.path(dirname(run_dir), "modules", "runner_helpers.R"))
source(file.path(dirname(run_dir), "modules", "mapping_helpers.R"))

normalize_lapCorrected <- function(slice, libSize, Nconstant, Dconstant) {
  division_value <- (colSums(slice) + Nconstant) / (libSize + Dconstant)
  sweep(slice, 2, division_value, "/")
}

InstaPrism_normalize <- function(data, normalize_method = "lapCorrected") {
  stopifnot(normalize_method %in% c("lapCorrected", "fracDiv"))
  scaled_array_3d <- array(dim = dim(data))
  dimnames(scaled_array_3d) <- dimnames(data)
  libSize <- apply(data, MARGIN = 2, sum)

  if ((range(libSize)[2] - range(libSize)[1]) > 0.01 * mean(libSize)) {
    warning("inconsistent libsize between reconstructed bulk profiles")
  }

  k <- dim(data)[[3]]
  Nconstant <- 0.1 * libSize / k
  Dconstant <- 0.1 * libSize / k

  for (cell_type in dimnames(data)[[3]]) {
    slice <- data[, , cell_type]
    if (normalize_method == "lapCorrected") {
      scaled_slice <- normalize_lapCorrected(slice, libSize, Nconstant, Dconstant)
    } else {
      scaled_slice <- normalize_lapCorrected(slice, libSize, Nconstant = 0, Dconstant = 0)
    }
    scaled_array_3d[, , cell_type] <- scaled_slice
  }
  scaled_array_3d
}

instaprism_key <- function(dataset, refType, ref_dir) {
  key_table <- data.frame(
    dataset = c(
      "BRCA_Bassez2021", "BRCA_Wu2021",
      "CRC_Lee2020", "CRC_Pelka2021",
      "LUAD_Kim2020", "NSCLC_Wu2021", "LUAD_Laughney2020"
    ),
    key = c(
      "Cancer_cell", "Cancer_Epithelial",
      "Malignant", "Malignant",
      "Malignant", "Malignant", "Malignant"
    ),
    stringsAsFactors = FALSE
  )

  if (grepl("PBMC", dataset) || grepl("ROSMAP", dataset)) {
    return(NA)
  }
  key_dataset <- if (refType == "indep") basename(ref_dir) else dataset
  key <- key_table$key[key_table$dataset == key_dataset]
  if (length(key) == 0) NA else key[[1]]
}

parser <- ArgumentParser(add_help = TRUE)
parser$add_argument("--dataset", type = "character", required = TRUE)
parser$add_argument("--config_id", type = "character", required = TRUE)
parser$add_argument("--n_core", type = "integer", required = FALSE, default = 15)
parser$add_argument("--extra_args", type = "character", required = FALSE, default = "")

args <- parser$parse_args()
repo_root <- find_repo_root(run_dir)
extra <- apply_method_default_extra_args(parse_extra_args(args$extra_args), "InstaPrism", repo_root = repo_root)
validate_extra_args(
  extra,
  c("n_iter", "run_update", "map_cell_types", "use_test_samples", "export_raw"),
  "InstaPrism"
)

dataset <- args$dataset
config_id <- args$config_id
n_core <- args$n_core
n_iter <- extra_arg(extra, "n_iter", 1000)
run_update <- extra_arg(extra, "run_update", FALSE)
map_cell_types <- extra_arg(extra, "map_cell_types", TRUE)
use_test_samples <- extra_arg(extra, "use_test_samples", TRUE)
export_raw <- extra_arg(extra, "export_raw", FALSE)

paths <- resolve_deconv_paths(dataset, config_id, repo_root = repo_root)
bulk_input <- paths$config$bulk_input[[1]]
bulk_scale <- paths$config$bulk_scale[[1]]
bulk_normalization <- paths$config$bulk_normalization[[1]]
frac_input <- paths$config$frac_input[[1]]
refType <- paths$refType
message_log_path <- method_message_log_path(paths, "InstaPrism")
cleanup_message_log <- start_message_log(message_log_path)
on.exit(cleanup_message_log(), add = TRUE)

if (frac_input != "InstaPrismfrac") {
  warning(
    "run_InstaPrism.R does not use frac_input. ",
    "Received frac_input=", frac_input, " for config_id=", config_id, ". ",
    "Results will be generated under this config, but they are expected to be ",
    "equivalent to running the same bulk_input/refType with frac_input=InstaPrismfrac."
  )
}

refPhi_path <- ref_file(paths$ref_dir, "refPhi.RDS", "refPhi")
bulk_expr <- read_bulk(paths)
bulk_prep <- prepare_bulk_for_deconv(bulk_expr, paths, method = "InstaPrism")
bulk_expr <- bulk_prep$bulk_expr
sample_restriction <- restrict_to_test_samples(
  bulk_expr,
  paths,
  use_test_samples = use_test_samples
)
bulk_expr <- sample_restriction$bulk
bulk_input_dim <- format_deconv_dim(bulk_expr)
deconv_input_dim <- bulk_input_dim
refPhi <- readRDS(refPhi_path)
key <- instaprism_key(dataset, refType, paths$ref_dir)

validate_non_log_bulk_input(bulk_expr, bulk_input, config_id, method = "InstaPrism")

message("Running InstaPrism")
message("  dataset: ", dataset)
message("  config_id: ", config_id)
message("  bulk_input: ", bulk_input, " (", bulk_input_dim, ")")
message("  bulk_scale: ", bulk_scale)
message("  bulk_normalization: ", bulk_normalization)
message("  bulk_preparation: ", bulk_prep$action)
message("  refType: ", refType)
message("  config_frac_input: ", frac_input)
message("  n_iter: ", n_iter)
message("  run_update: ", run_update)
message("  map_cell_types: ", map_cell_types)
message("  use_test_samples: ", use_test_samples)
message("  export_raw: ", export_raw)
message("  note: base InstaPrism does not read frac_input")

message("Starting deconvolution")
message("  deconv_input: ", deconv_input_dim)

instaprism_runtime <- with_runtime_log(paths, "InstaPrism", n_core = n_core, deconv_input = deconv_input_dim, expr = {
  InstaPrism_res <- InstaPrism(
    bulk_Expr = bulk_expr,
    refPhi_cs = refPhi,
    n.iter = n_iter,
    filter = FALSE,
    n.core = n_core
  )
  inferred_frac <- t(InstaPrism_res@Post.ini.ct@theta)
  Z <- get_Z_array(InstaPrism_res, resolution = "ct", n.core = n_core)
  Z <- aperm(Z, c(2, 1, 3))
  Z_normalized <- InstaPrism_normalize(Z, normalize_method = "lapCorrected")
  if (map_cell_types) {
    Z_normalized <- apply_indep_ref_cell_type_mapping(Z_normalized, paths, repo_root = repo_root)
  }

  message_existing_method_output(paths, "InstaPrism")
  write_ctse_array(Z_normalized, paths$output_dir, "InstaPrism")
  write_method_fraction(inferred_frac, paths, "InstaPrism", "InstaPrismfrac")
  if (export_raw) {
    Z_raw <- Z
    if (map_cell_types) {
      Z_raw <- apply_indep_ref_cell_type_mapping(Z_raw, paths, repo_root = repo_root)
    }
    message_existing_method_output(paths, "InstaPrism_raw")
    write_ctse_array(Z_raw, paths$output_dir, "InstaPrism_raw")
  }

  list(InstaPrism_res = InstaPrism_res, bulk_expr = bulk_expr, key = key)
})
InstaPrism_res <- instaprism_runtime$value$InstaPrism_res
bulk_expr <- instaprism_runtime$value$bulk_expr
key <- instaprism_runtime$value$key

if (run_update) {
  with_runtime_log(paths, "InstaPrismUpdated", n_core = n_core, deconv_input = deconv_input_dim, expr = {
    InstaPrism_res_updated <- InstaPrism_update(
      InstaPrism_res,
      bulk_expr,
      key = key,
      n.iter = n_iter,
      n.core = n_core
    )
    Z_updated <- get_Z_array(InstaPrism_res_updated, resolution = "ct", n.core = n_core)
    Z_updated <- aperm(Z_updated, c(2, 1, 3))
    updated_frac <- t(InstaPrism_res_updated@theta)
    Z_updated_normalized <- InstaPrism_normalize(Z_updated, normalize_method = "lapCorrected")
    if (map_cell_types) {
      Z_updated_normalized <- apply_indep_ref_cell_type_mapping(Z_updated_normalized, paths, repo_root = repo_root)
    }

    message_existing_method_output(paths, "InstaPrismUpdated")
    write_ctse_array(Z_updated_normalized, paths$output_dir, "InstaPrismUpdated")
    write_method_fraction(updated_frac, paths, "InstaPrismUpdated", "InstaPrismfrac")
  })
}

message("Output saved")
message(
  "  ctse: ",
  file.path(paths$output_dir, "InstaPrism", "<cell_type>.txt.gz")
)

if (run_update) {
  message(
    "  ctse updated: ",
    file.path(paths$output_dir, "InstaPrismUpdated", "<cell_type>.txt.gz")
  )
}

if (export_raw) {
  message(
    "  raw ctse: ",
    file.path(paths$output_dir, "InstaPrism_raw", "<cell_type>.txt.gz")
  )
}

message(
  "  frac: ",
  file.path(paths$output_dir, "InstaPrism", "InstaPrismfrac.txt")
)


if (run_update) {
  message(
    "  updated frac: ",
    file.path(paths$output_dir, "InstaPrismUpdated", "InstaPrismfrac.txt")
  )
}

if (map_cell_types && refType == "indep") {
  message("  ctse naming note: mapped cell types use benchmark names; unmapped reference cell types keep reference names")
} else {
  message("  ctse naming note: exported cell types use reference cell type names")
}

message("done")
