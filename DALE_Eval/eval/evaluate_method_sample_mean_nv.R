#!/usr/bin/env Rscript

args0 <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args0, value = TRUE)
script_path <- if (length(file_arg) > 0L) {
  normalizePath(sub("^--file=", "", file_arg[[1L]]), mustWork = TRUE)
} else {
  normalizePath(
    "DALE_Eval/eval/evaluate_method_sample_mean_nv.R",
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
source(file.path(repo_root, "DALE_Eval", "modules", "mapping_helpers.R"))
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
    "  Rscript DALE_Eval/eval/evaluate_method_sample_mean_nv.R \\\n",
    "    --dataset PBMC_Perez2022 --config_id config01 \\\n",
    "    --methods all --lib_norm false --digits 6 [--dry_run]\n\n",
    "Options:\n",
    "  --dataset     Benchmark dataset name. Required.\n",
    "  --config_id   Active evaluation config ID. Required.\n",
    "  --methods     Comma-separated method names, or all. Default: all.\n",
    "  --lib_norm    true or false. Default: false.\n",
    "                sample_logmean is written only when false.\n",
    "  --digits      Decimal places written to disk. Default: 6.\n",
    "  --dry_run     Resolve inputs and outputs without reading CTSE matrices or writing results.\n",
    "  --help, -h    Show this help.\n",
    sep = ""
  )
}

parse_boolean <- function(value, label) {
  value <- tolower(trimws(as.character(value)))
  if (!value %in% c("true", "false")) {
    stop(label, " must be true or false")
  }
  identical(value, "true")
}

parse_cli_args <- function(args) {
  out <- list(
    methods = "all",
    lib_norm = "false",
    digits = "6",
    dry_run = FALSE
  )
  allowed <- c(
    "dataset",
    "config_id",
    "methods",
    "lib_norm",
    "digits",
    "dry_run"
  )
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
  required <- c("dataset", "config_id")
  missing <- required[!vapply(
    required,
    function(x) !is.null(out[[x]]) && nzchar(out[[x]]),
    logical(1)
  )]
  if (length(missing) > 0L) {
    usage()
    stop("Missing required arguments: ", paste(missing, collapse = ", "))
  }
  out$lib_norm <- parse_boolean(out$lib_norm, "--lib_norm")
  out$digits <- suppressWarnings(as.integer(out$digits))
  if (is.na(out$digits) || out$digits < 0L) {
    stop("--digits must be a non-negative integer")
  }
  out
}

split_csv <- function(x) {
  values <- trimws(unlist(strsplit(x, ",", fixed = TRUE)))
  unique(values[nzchar(values)])
}

args <- parse_cli_args(commandArgs(trailingOnly = TRUE))
eval_config <- read_pairwise_eval_config(args$config_id, repo_root = repo_root)
paths <- resolve_deconv_paths(args$dataset, args$config_id, repo_root = repo_root)
if (!dir.exists(paths$output_dir)) {
  stop("Deconvolution result directory not found: ", paths$output_dir)
}
expected_truth_type <- trimws(as.character(
  eval_config$expected_truth_type[[1L]]
))
performance_slug <- sample_stats_method_performance_slug(
  paths,
  expected_truth_type
)

available_methods <- sort(list.dirs(
  paths$output_dir,
  recursive = FALSE,
  full.names = FALSE
))
methods <- if (identical(args$methods, "all")) {
  available_methods
} else {
  requested <- split_csv(args$methods)
  missing <- setdiff(requested, available_methods)
  if (length(missing) > 0L) {
    stop("Requested methods not found: ", paste(missing, collapse = ", "))
  }
  requested
}
if (length(methods) == 0L) {
  stop("No methods selected")
}

scitd_config <- read_scitd_config(repo_root = repo_root)
method_preprocessing <- read_scitd_method_preprocessing(repo_root = repo_root)
test_info <- read_pairwise_test_samples(paths$obj_dir)

message("Evaluating inferred CTSE sample mean and normalized variance")
message("  dataset: ", args$dataset)
message("  config_id: ", args$config_id)
message("  config_slug: ", config_slug(paths$config))
message(
  "  default_truth_type: ", expected_truth_type,
  " (output routing only; truth is not an input)"
)
message("  performance_slug: ", performance_slug)
message("  methods: ", paste(methods, collapse = ", "))
message("  lib_norm: ", args$lib_norm)
message("  sample_policy: all method samples; require testing-only and aligned")
message("  gene_policy: all genes shared across selected method cell types")
message("  zero_profile_filter: none")
message("  variable_gene_filter: none")
message("  normalized_variance_qc: none")
message("  digits: ", args$digits)
message("  dry_run: ", args$dry_run)

for (method in methods) {
  method_dir <- file.path(paths$output_dir, method)
  method_files <- ctse_cell_type_files(method_dir)
  if (length(method_files) == 0L) {
    stop("No CTSE .txt.gz files found for method: ", method)
  }
  selection <- resolve_pairwise_method_cell_types(
    paths,
    method_files,
    repo_root = repo_root
  )
  preprocess_mode <- scitd_method_preprocess_mode(
    method,
    method_preprocessing
  )
  if (args$lib_norm && !identical(preprocess_mode, "log1p")) {
    stop(
      "lib_norm=true is not supported for ",
      method,
      " (preprocess_mode=",
      preprocess_mode,
      ")"
    )
  }
  output_paths <- sample_stats_method_output_paths(
    paths,
    method,
    expected_truth_type = expected_truth_type,
    lib_norm = args$lib_norm
  )

  message("\nMethod: ", method)
  message("  preprocess_mode: ", preprocess_mode)
  message("  selected_cell_types: ", paste(selection$selected, collapse = ", "))
  if (length(selection$extra) > 0L) {
    message(
      "  ignored_unmatched_cell_types: ",
      paste(sort(selection$extra), collapse = ", ")
    )
  }
  message("  sample_mean_output: ", output_paths$sample_mean)
  if (args$lib_norm) {
    message("  sample_logmean_output: skipped for lib_norm=true")
  } else {
    message("  sample_logmean_output: ", output_paths$sample_logmean)
  }
  message("  normalized_variance_output: ", output_paths$normalized_variance)
  if (isTRUE(args$dry_run)) {
    next
  }

  sample_stats_enable_vendored_scitd(repo_root)
  prepared <- prepare_method_sample_stats_input(
    method_files,
    selection$selected,
    test_samples = test_info$samples,
    label = paste0("method ", method)
  )
  sample_logmean_result <- if (args$lib_norm) {
    NULL
  } else {
    calculate_sample_logmean(
      prepared$matrices,
      genes = prepared$genes,
      cell_types = selection$selected,
      preprocess_mode = preprocess_mode
    )
  }
  processed <- preprocess_method_sample_stats(
    prepared$matrices,
    preprocess_mode = preprocess_mode,
    scale_factor = scitd_config$scale_factor,
    lib_norm = args$lib_norm,
    zero_tolerance = scitd_config$zero_profile_tolerance
  )
  result <- calculate_sample_mean_nv(
    processed$y,
    genes = prepared$genes,
    cell_types = selection$selected
  )
  write_sample_stats_matrix(
    result$sample_mean,
    output_paths$sample_mean,
    digits = args$digits
  )
  if (!args$lib_norm) {
    write_sample_stats_matrix(
      sample_logmean_result$sample_logmean,
      output_paths$sample_logmean,
      digits = args$digits
    )
  }
  write_sample_stats_matrix(
    result$normalized_variance,
    output_paths$normalized_variance,
    digits = args$digits
  )
  message(
    "  wrote ", length(prepared$genes), " genes x ",
    length(selection$selected), " cell types using ",
    length(prepared$samples), " samples"
  )
  message("  negative gene/cell-type shifts: ", nrow(processed$negative_shifts))
  if (!args$lib_norm) {
    message(
      "  negative raw-mean gene shifts for sample_logmean: ",
      nrow(sample_logmean_result$negative_mean_shifts)
    )
  }
  message("  retained all-zero sample/cell-type profiles: ", nrow(processed$zero_profiles))
  message(
    "  non-finite NV values written as NA by cell type: ",
    paste(
      paste0(names(result$invalid_nv_by_cell_type), "=", result$invalid_nv_by_cell_type),
      collapse = ", "
    )
  )
  if (nrow(result$nv_gam_skipped) > 0L) {
    message("  normalized-variance GAM skipped by cell type:")
    for (i in seq_len(nrow(result$nv_gam_skipped))) {
      skipped <- result$nv_gam_skipped[i, , drop = FALSE]
      message(
        "    ", skipped$cell_type,
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
