#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(argparse)
})

args0 <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args0, value = TRUE)
run_dir <- if (length(file_arg) > 0) dirname(normalizePath(sub("^--file=", "", file_arg[[1]]))) else getwd()
source(file.path(dirname(run_dir), "modules", "config_helpers.R"))
source(file.path(dirname(run_dir), "modules", "runner_helpers.R"))
source(file.path(dirname(run_dir), "modules", "marker_helpers.R"))
source(file.path(dirname(run_dir), "modules", "mapping_helpers.R"))
read_optional_covariance <- function(ref_dir) {
  covariance_path <- file.path(ref_dir, "bMIND_covariance.RDS")
  if (!file.exists(covariance_path)) {
    return(NULL)
  }
  readRDS(covariance_path)
}

subset_covariance <- function(covariance, genes, cell_types) {
  if (is.null(covariance)) {
    return(NULL)
  }
  if (length(dim(covariance)) != 3) {
    stop("bMIND_covariance.RDS must be a 3D array")
  }

  dn <- dimnames(covariance)
  if (!is.null(dn[[1]])) {
    missing_genes <- setdiff(genes, dn[[1]])
    if (length(missing_genes) > 0) {
      stop("Covariance missing genes used for deconvolution: ", paste(head(missing_genes, 10), collapse = ", "))
    }
    covariance <- covariance[genes, , , drop = FALSE]
    dn <- dimnames(covariance)
  }

  if (!is.null(dn[[2]])) {
    missing_cell_types <- setdiff(cell_types, dn[[2]])
    if (length(missing_cell_types) > 0) {
      stop("Covariance missing cell types used for deconvolution: ", paste(missing_cell_types, collapse = ", "))
    }
    covariance <- covariance[, cell_types, , drop = FALSE]
    dn <- dimnames(covariance)
  }

  if (!is.null(dn[[3]])) {
    missing_cell_types <- setdiff(cell_types, dn[[3]])
    if (length(missing_cell_types) > 0) {
      stop("Covariance missing cell types used for deconvolution: ", paste(missing_cell_types, collapse = ", "))
    }
    covariance <- covariance[, , cell_types, drop = FALSE]
  }

  covariance
}

parser <- ArgumentParser(add_help = TRUE)
parser$add_argument("--dataset", type = "character", required = TRUE)
parser$add_argument("--config_id", type = "character", required = TRUE)
parser$add_argument("--n_core", type = "integer", required = FALSE, default = 15)
parser$add_argument("--extra_args", type = "character", required = FALSE, default = "")

args <- parser$parse_args()
repo_root <- find_repo_root(run_dir)
extra <- apply_method_default_extra_args(parse_extra_args(args$extra_args), "bMIND", repo_root = repo_root)
validate_extra_args(
  extra,
  c("map_cell_types", "use_test_samples", "export_posterior", "run_epicunmix", "use_limma_top_genes", "top_n"),
  "bMIND"
)

dataset <- args$dataset
config_id <- args$config_id
n_core <- args$n_core
map_cell_types <- extra_arg(extra, "map_cell_types", TRUE)
use_test_samples <- extra_arg(extra, "use_test_samples", TRUE)
export_posterior <- extra_arg(extra, "export_posterior", FALSE)
run_epicunmix <- extra_arg(extra, "run_epicunmix", TRUE)
use_limma_top_genes <- extra_arg(extra, "use_limma_top_genes", FALSE)
top_n <- extra_arg(extra, "top_n", 100)
validate_top_n_arg(extra, use_limma_top_genes, "bMIND")

paths <- resolve_deconv_paths(dataset, config_id, repo_root = repo_root)
bulk_input <- paths$config$bulk_input[[1]]
bulk_scale <- paths$config$bulk_scale[[1]]
bulk_normalization <- paths$config$bulk_normalization[[1]]
frac_input <- paths$config$frac_input[[1]]
refType <- paths$refType
message_log_path <- method_message_log_path(paths, "bMIND")
cleanup_bmind_log <- start_message_log(message_log_path)
on.exit(cleanup_bmind_log(), add = TRUE)

bulk_expr <- read_bulk(paths)
bulk_prep <- prepare_bulk_for_deconv(bulk_expr, paths, method = "bMIND")
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
profile <- read_reference_matrix(paths$ref_dir, c("bMIND_profile.csv", "bMIND_profile.txt"), "bMIND profile")
bMIND_covariance <- read_optional_covariance(paths$ref_dir)

aligned <- align_fraction_and_reference_cell_types(
  frac,
  colnames(profile),
  paths,
  repo_root
)
frac <- aligned$frac
profile <- profile[, aligned$reference_cell_types, drop = FALSE]

log_bulk <- prepare_log_bulk_input(bulk_expr, bulk_input, config_id, method = "bMIND")
bulk_expr_logged <- log_bulk$bulk_expr

limma_gene_count <- NA_integer_
if (use_limma_top_genes) {
  limma_genes <- read_limma_top_gene_union(paths, top_n = top_n)
  limma_gene_count <- length(limma_genes)
} else {
  limma_genes <- NULL
}

genes <- intersect(rownames(bulk_expr_logged), rownames(profile))
if (use_limma_top_genes) {
  genes <- intersect(genes, limma_genes)
}
if (length(genes) == 0) {
  stop("No genes available for bMIND after intersecting bulk, profile, and configured gene filters")
}

bulk_expr_logged <- bulk_expr_logged[genes, , drop = FALSE]
profile <- profile[genes, , drop = FALSE]
bMIND_covariance <- subset_covariance(bMIND_covariance, genes, colnames(profile))
deconv_input_dim <- format_deconv_dim(bulk_expr_logged)

if (!identical(rownames(bulk_expr_logged), rownames(profile))) {
  stop("bMIND bulk and profile gene names are not aligned after subsetting")
}
if (!identical(colnames(profile), colnames(frac))) {
  stop("bMIND profile and fraction cell types are not aligned after subsetting")
}

message(if (run_epicunmix) "Running bMIND and EPICunmix" else "Running bMIND")
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
message("  export_posterior: ", export_posterior)
message("  run_epicunmix: ", run_epicunmix)
if (!run_epicunmix && !export_posterior) {
  warning(
    "run_epicunmix=false and export_posterior=false: this run will write bMIND CTSE only. ",
    "No bMIND_posterior.RDS will be available for run_EPICunmix_with_posterior.R."
  )
}
message("  use_limma_top_genes: ", use_limma_top_genes)
if (use_limma_top_genes) {
  message("  top_n: ", top_n)
}
message("  profile_dim: ", nrow(profile), " genes x ", ncol(profile), " cell types")
message("  note: ", log_bulk$note)

message("Starting deconvolution")
message("  deconv_input: ", deconv_input_dim)


bmind_runtime <- with_runtime_log(paths, "bMIND", n_core = n_core, deconv_input = deconv_input_dim, expr = {
  posterior <- MIND::bMIND(
    bulk_expr_logged,
    frac,
    profile = profile,
    covariance = bMIND_covariance,
    ncore = n_core
  )
  Z <- aperm(posterior$A, c(1, 3, 2))
  dimnames(Z)[[2]] <- colnames(bulk_expr_logged)
  dimnames(Z)[[3]] <- colnames(frac)

  Z_export <- Z
  if (map_cell_types) {
    Z_export <- apply_indep_ref_cell_type_mapping(Z_export, paths, repo_root = repo_root)
  }
  message_existing_method_output(paths, "bMIND")
  write_ctse_array(Z_export, paths$output_dir, "bMIND")

  if (export_posterior) {
    posterior_dir <- file.path(paths$output_dir, "bMIND")
    dir.create(posterior_dir, recursive = TRUE, showWarnings = FALSE)
    saveRDS(posterior, file = file.path(posterior_dir, "bMIND_posterior.RDS"))
  }

  list(
    posterior = posterior,
    Z = Z,
    bulk_expr_logged = bulk_expr_logged,
    frac = frac
  )
})
posterior <- bmind_runtime$value$posterior
Z <- bmind_runtime$value$Z
bulk_expr_logged <- bmind_runtime$value$bulk_expr_logged
frac <- bmind_runtime$value$frac
cleanup_bmind_log()
options(ctse.message_log_path = "")

if (run_epicunmix) {
  message_log_path <- method_message_log_path(paths, "EPICunmix")
  cleanup_epicunmix_log <- start_message_log(message_log_path)
  on.exit(cleanup_epicunmix_log(), add = TRUE)
  message("Running EPICunmix after bMIND")
  message("  dataset: ", dataset)
  message("  config_id: ", config_id)
  message("  bulk_input: ", bulk_input, " (", bulk_input_dim, ")")
message("  bulk_scale: ", bulk_scale)
message("  bulk_normalization: ", bulk_normalization)
message("  bulk_preparation: ", bulk_prep$action)
  message("  refType: ", refType)
  message("  frac_input: ", frac_input, " (", nrow(frac), " samples x ", ncol(frac), " cell types)")
  message("  n_core: ", n_core)
  message("  deconv_input: ", deconv_input_dim)

  with_runtime_log(paths, "EPICunmix", n_core = n_core, deconv_input = deconv_input_dim, expr = {
    epic_unmix <- EPICunmix::run_epic_unmix(
      bulk_expr_logged,
      frac,
      posterior,
      outf = FALSE,
      ncore = n_core
    )
    Z2 <- aperm(epic_unmix$A, c(1, 3, 2))
    dimnames(Z2)[[2]] <- colnames(bulk_expr_logged)
    dimnames(Z2)[[3]] <- dimnames(Z)[[3]]

    Z2_export <- Z2
    if (map_cell_types) {
      Z2_export <- apply_indep_ref_cell_type_mapping(Z2_export, paths, repo_root = repo_root)
    }
    message_existing_method_output(paths, "EPICunmix")
    write_ctse_array(Z2_export, paths$output_dir, "EPICunmix")
  })
}

message("Output saved")
message("  bMIND ctse: ", file.path(paths$output_dir, "bMIND", "<cell_type>.txt.gz"))
if (run_epicunmix) {
  message("  EPICunmix ctse: ", file.path(paths$output_dir, "EPICunmix", "<cell_type>.txt.gz"))
}
if (export_posterior) {
  message("  bMIND posterior: ", file.path(paths$output_dir, "bMIND", "bMIND_posterior.RDS"))
} else if (!run_epicunmix) {
  message("  bMIND posterior: not saved; run_EPICunmix_with_posterior.R cannot use this run")
}
if (map_cell_types && refType == "indep") {
  message("  ctse naming note: mapped cell types use benchmark names; unmapped reference cell types keep reference names")
} else {
  message("  ctse naming note: exported cell types use reference cell type names")
}
message("done")
