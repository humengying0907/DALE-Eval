#!/usr/bin/env Rscript

args0 <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args0, value = TRUE)
script_path <- if (length(file_arg) > 0L) {
  normalizePath(sub("^--file=", "", file_arg[[1L]]), mustWork = TRUE)
} else {
  normalizePath(
    "DALE_Eval/eval/evaluate_truth_sample_mean_nv.R",
    mustWork = FALSE
  )
}

find_repo_root_local <- function(start) {
  cur <- normalizePath(start, mustWork = TRUE)
  repeat {
    if (
      file.exists(file.path(cur, "DALE_Eval", "configs", "deconv_configs.txt")) &&
        dir.exists(file.path(cur, "DALE_Eval"))
    ) {
      return(cur)
    }
    parent <- dirname(cur)
    if (identical(parent, cur)) {
      stop("Could not find repository root from: ", start)
    }
    cur <- parent
  }
}

repo_root <- find_repo_root_local(dirname(script_path))
source(file.path(repo_root, "DALE_Eval", "modules", "config_helpers.R"))
source(file.path(repo_root, "DALE_Eval", "modules", "evalu.R"))
source(file.path(repo_root, "DALE_Eval", "modules", "scITD_helpers.R"))
source(file.path(
  repo_root,
  "DALE_Eval",
  "modules",
  "celltype_pairwise_cor_helpers.R"
))
source(file.path(
  repo_root,
  "DALE_Eval",
  "modules",
  "sample_mean_nv_helpers.R"
))

usage <- function() {
  cat(
    "Usage:\n",
    "  Rscript DALE_Eval/eval/evaluate_truth_sample_mean_nv.R \\\n",
    "    --dataset PBMC_Perez2022 --truth_type sumcount \\\n",
    "    --digits 6 [--dry_run]\n\n",
    "Options:\n",
    "  --dataset     Benchmark dataset name. Required.\n",
    "  --truth_type  sumcount, meancpm, or sumcount_cpm. Required.\n",
    "  --digits      Decimal places written to disk. Default: 6.\n",
    "  --dry_run     Resolve inputs and outputs without reading CTSE matrices or writing results.\n",
    "  --help, -h    Show this help.\n",
    sep = ""
  )
}

parse_cli_args <- function(args) {
  out <- list(digits = "6", dry_run = FALSE)
  allowed <- c("dataset", "truth_type", "digits", "dry_run")
  i <- 1L
  while (i <= length(args)) {
    raw_key <- args[[i]]
    if (raw_key %in% c("--help", "-h")) {
      usage()
      quit(status = 0L)
    }
    if (!startsWith(raw_key, "--")) {
      stop("Unexpected positional argument: ", raw_key)
    }
    key <- gsub("-", "_", sub("^--", "", raw_key))
    if (!key %in% allowed) {
      stop("Unsupported argument: ", raw_key)
    }
    if (identical(key, "dry_run")) {
      out[[key]] <- TRUE
      i <- i + 1L
      next
    }
    if (i == length(args) || startsWith(args[[i + 1L]], "--")) {
      stop("Missing value for ", raw_key)
    }
    out[[key]] <- args[[i + 1L]]
    i <- i + 2L
  }
  required <- c("dataset", "truth_type")
  missing <- required[!vapply(
    required,
    function(x) !is.null(out[[x]]) && nzchar(out[[x]]),
    logical(1)
  )]
  if (length(missing) > 0L) {
    usage()
    stop("Missing required arguments: ", paste(missing, collapse = ", "))
  }
  supported_truth_types <- c("sumcount", "meancpm", "sumcount_cpm")
  if (!out$truth_type %in% supported_truth_types) {
    stop(
      "Unsupported --truth_type=", out$truth_type,
      ". Supported values: ",
      paste(supported_truth_types, collapse = ", ")
    )
  }
  out$digits <- suppressWarnings(as.integer(out$digits))
  if (is.na(out$digits) || out$digits < 0L) {
    stop("--digits must be a non-negative integer")
  }
  out
}

args <- parse_cli_args(commandArgs(trailingOnly = TRUE))
obj_dir <- file.path(repo_root, "Benchmarking_obj", args$dataset)
if (!dir.exists(obj_dir)) {
  stop("Benchmark object not found: ", obj_dir)
}
scitd_config <- read_scitd_config(repo_root = repo_root)
truth_files <- scitd_truth_files(
  args$dataset,
  truth_type = args$truth_type,
  repo_root = repo_root
)
cell_types <- sort(names(truth_files))
validate_pairwise_cell_types(cell_types, label = "truth cell types")
test_info <- read_pairwise_test_samples(obj_dir)
output_paths <- sample_stats_truth_output_paths(
  args$dataset,
  truth_type = args$truth_type,
  repo_root = repo_root
)

message("Evaluating truth sample mean and normalized variance")
message("  dataset: ", args$dataset)
message("  truth_type: ", args$truth_type)
message("  truth_input: ctse_truth/", args$truth_type)
message("  selected_cell_types: ", paste(cell_types, collapse = ", "))
message("  sample_policy: all testing samples")
message("  n_test_samples_defined: ", length(test_info$samples))
message("  gene_policy: all genes shared across truth cell types")
message("  zero_profile_filter: none; zero profiles remain zero in Y")
message("  variable_gene_filter: none")
message("  normalized_variance_qc: none")
if (identical(args$truth_type, "sumcount")) {
  message(
    "  preprocessing: TMM normalization to ",
    scitd_config$scale_factor,
    " then log1p"
  )
} else {
  message("  preprocessing: log1p only")
}
message("  digits: ", args$digits)
message("  sample_mean_output: ", output_paths$sample_mean)
message("  normalized_variance_output: ", output_paths$normalized_variance)
message("  dry_run: ", args$dry_run)

if (!isTRUE(args$dry_run)) {
  sample_stats_enable_vendored_scitd(repo_root)
  prepared <- prepare_truth_sample_stats_input(
    truth_files,
    cell_types,
    test_samples = test_info$samples,
    truth_type = args$truth_type
  )
  processed <- preprocess_truth_sample_stats(
    prepared$matrices,
    truth_type = args$truth_type,
    scale_factor = scitd_config$scale_factor,
    norm_method = scitd_config$truth_norm_method,
    zero_tolerance = scitd_config$zero_profile_tolerance
  )
  result <- calculate_sample_mean_nv(
    processed$y,
    genes = prepared$genes,
    cell_types = cell_types
  )
  write_sample_stats_matrix(
    result$sample_mean,
    output_paths$sample_mean,
    digits = args$digits
  )
  write_sample_stats_matrix(
    result$normalized_variance,
    output_paths$normalized_variance,
    digits = args$digits
  )
  message(
    "wrote ", length(prepared$genes), " genes x ", length(cell_types),
    " cell types using ", length(prepared$samples), " testing samples"
  )
  message(
    "retained all-zero sample/cell-type profiles: ",
    nrow(processed$zero_profiles)
  )
  message(
    "non-finite NV values written as NA by cell type: ",
    paste(
      paste0(names(result$invalid_nv_by_cell_type), "=", result$invalid_nv_by_cell_type),
      collapse = ", "
    )
  )
  if (nrow(result$nv_gam_skipped) > 0L) {
    message("normalized-variance GAM skipped by cell type:")
    for (i in seq_len(nrow(result$nv_gam_skipped))) {
      skipped <- result$nv_gam_skipped[i, , drop = FALSE]
      message(
        "  ", skipped$cell_type,
        ": ", skipped$reason,
        " (usable_genes=", skipped$n_usable_genes,
        ", unique_log_means=", skipped$n_unique_log_means,
        ", k=", skipped$gam_k,
        "); normalized variance written as NA"
      )
    }
  }
}

message(if (isTRUE(args$dry_run)) "dry run complete" else "done")
