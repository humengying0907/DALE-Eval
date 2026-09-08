#!/usr/bin/env Rscript

args0 <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args0, value = TRUE)
script_path <- if (length(file_arg) > 0L) {
  normalizePath(sub("^--file=", "", file_arg[[1L]]), mustWork = TRUE)
} else {
  normalizePath(
    "DALE_Eval/eval/build_gene_prioritization_scores.R",
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
source(file.path(repo_root, "DALE_Eval", "modules", "runner_helpers.R"))
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
source(file.path(
  repo_root,
  "DALE_Eval",
  "modules",
  "gene_prioritization_helpers.R"
))

usage <- function() {
  cat(
    "Usage:\n",
    "  Rscript DALE_Eval/eval/build_gene_prioritization_scores.R \\\n",
    "    --dataset LUAD_Kim2020 --config_id config01 \\\n",
    "    --scores all --methods all --lib_norm false --n_core 15 --digits 6 [--dry_run]\n\n",
    "Options:\n",
    "  --dataset     Benchmark dataset name. Required.\n",
    "  --config_id   Deconvolution config ID. Required.\n",
    "  --scores      Comma-separated score names, or all. Default: all.\n",
    "                Names: meanlog_margin, meanlog_mean_contrast,\n",
    "                logmean_margin, logmean_contrast, partial_r2,\n",
    "                genesigtest_adapted_fdr.\n",
    "  --methods     Comma-separated method names, or all. Default: all.\n",
    "  --lib_norm    Select matching sample-mean preprocessing: true or false.\n",
    "                Default: false.\n",
    "  --n_core      Workers used for partial R2. Default: 1.\n",
    "  --digits      Decimal places for method-specific scores and significant\n",
    "                digits for GeneSigTest_adapted FDR. Default: 6.\n",
    "  --dry_run     Resolve inputs and outputs without reading matrices or writing results.\n",
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

split_csv <- function(x) {
  values <- trimws(unlist(strsplit(x, ",", fixed = TRUE)))
  unique(values[nzchar(values)])
}

parse_cli_args <- function(args) {
  out <- list(
    scores = "all",
    methods = "all",
    lib_norm = "false",
    n_core = "1",
    digits = "6",
    dry_run = FALSE
  )
  allowed <- c(
    "dataset",
    "config_id",
    "scores",
    "methods",
    "lib_norm",
    "n_core",
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
  out$n_core <- suppressWarnings(as.integer(out$n_core))
  if (is.na(out$n_core) || out$n_core < 1L) {
    stop("--n_core must be a positive integer")
  }
  out$digits <- suppressWarnings(as.integer(out$digits))
  if (is.na(out$digits) || out$digits < 1L) {
    stop("--digits must be a positive integer")
  }

  supported_scores <- gene_prioritization_score_names()
  logmean_scores <- c("logmean_margin", "logmean_contrast")
  requested_all <- identical(tolower(out$scores), "all")
  out$scores <- if (requested_all) {
    if (out$lib_norm) {
      setdiff(supported_scores, logmean_scores)
    } else {
      supported_scores
    }
  } else {
    requested <- split_csv(out$scores)
    unsupported <- setdiff(requested, supported_scores)
    if (length(unsupported) > 0L) {
      stop("Unsupported scores: ", paste(unsupported, collapse = ", "))
    }
    if (out$lib_norm && any(requested %in% logmean_scores)) {
      stop("logmean scores require --lib_norm false")
    }
    requested
  }
  if (length(out$scores) == 0L) {
    stop("No scores selected")
  }
  out
}

args <- parse_cli_args(commandArgs(trailingOnly = TRUE))
paths <- resolve_deconv_paths(args$dataset, args$config_id, repo_root = repo_root)
method_score_names <- c(
  "meanlog_margin",
  "meanlog_mean_contrast",
  "logmean_margin",
  "logmean_contrast",
  "partial_r2"
)
selected_method_scores <- intersect(args$scores, method_score_names)
needs_methods <- length(selected_method_scores) > 0L
needs_sample_mean <- any(c(
  "meanlog_margin",
  "meanlog_mean_contrast"
) %in% selected_method_scores)
needs_sample_logmean <- any(c(
  "logmean_margin",
  "logmean_contrast"
) %in% selected_method_scores)
needs_sample_stats <- needs_sample_mean || needs_sample_logmean
needs_partial_r2 <- "partial_r2" %in% selected_method_scores
needs_genesigtest <- "genesigtest_adapted_fdr" %in% args$scores

methods <- character(0)
if (needs_methods) {
  if (!dir.exists(paths$output_dir)) {
    stop("Deconvolution result directory not found: ", paths$output_dir)
  }
  available_methods <- sort(list.dirs(
    paths$output_dir,
    recursive = FALSE,
    full.names = FALSE
  ))
  methods <- if (identical(tolower(args$methods), "all")) {
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
}

expected_truth_type <- NULL
if (needs_sample_stats) {
  eval_config <- read_pairwise_eval_config(args$config_id, repo_root = repo_root)
  expected_truth_type <- trimws(as.character(
    eval_config$expected_truth_type[[1L]]
  ))
}

message("Building truth-independent gene-prioritization scores")
message("  dataset: ", args$dataset)
message("  config_id: ", args$config_id)
message("  config_slug: ", config_slug(paths$config))
message("  output_root: ", gene_prioritization_output_root(paths))
message("  scores: ", paste(args$scores, collapse = ", "))
if (needs_methods) {
  message("  methods: ", paste(methods, collapse = ", "))
  message("  lib_norm: ", args$lib_norm)
}
if (needs_sample_stats) {
  message(
    "  sample_stat_routing_truth_type: ", expected_truth_type,
    " (input path only; truth is not read)"
  )
}
message("  sample_policy: testing samples only")
if (needs_partial_r2) {
  message("  partial_r2_bulk_transform: log1p(config-prepared bulk)")
}
message("  n_core: ", args$n_core)
message("  digits: ", args$digits)
message("  dry_run: ", args$dry_run)

test_info <- NULL
method_preprocessing <- NULL
scitd_config <- NULL
if (needs_partial_r2) {
  test_info <- read_pairwise_test_samples(paths$obj_dir)
  method_preprocessing <- read_scitd_method_preprocessing(repo_root = repo_root)
  scitd_config <- read_scitd_config(repo_root = repo_root)
}

bulk_prepared <- NULL
if ((needs_partial_r2 || needs_genesigtest) && !isTRUE(args$dry_run)) {
  bulk_raw <- read_bulk(paths)
  bulk_prepared <- prepare_bulk_for_deconv(
    bulk_raw,
    paths,
    method = "Gene prioritization"
  )
  message("  bulk_preparation: ", bulk_prepared$action)
}

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
  output_paths <- gene_prioritization_method_output_paths(
    paths,
    method,
    lib_norm = args$lib_norm
  )

  sample_mean_path <- NULL
  if (needs_sample_mean) {
    sample_mean_path <- sample_stats_method_output_paths(
      paths,
      method,
      expected_truth_type = expected_truth_type,
      lib_norm = args$lib_norm
    )$sample_mean
  }

  sample_logmean_path <- NULL
  if (needs_sample_logmean) {
    sample_logmean_path <- sample_stats_method_output_paths(
      paths,
      method,
      expected_truth_type = expected_truth_type,
      lib_norm = FALSE
    )$sample_logmean
  }

  preprocess_mode <- NULL
  if (needs_partial_r2) {
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
  }

  message("\nMethod: ", method)
  message("  selected_cell_types: ", paste(selection$selected, collapse = ", "))
  if (length(selection$extra) > 0L) {
    message(
      "  ignored_unmatched_cell_types: ",
      paste(sort(selection$extra), collapse = ", ")
    )
  }
  if (needs_sample_mean) {
    message("  sample_mean_input: ", sample_mean_path)
  }
  if (needs_sample_logmean) {
    message("  sample_logmean_input: ", sample_logmean_path)
  }
  if (needs_partial_r2) {
    message("  preprocess_mode: ", preprocess_mode)
  }
  for (score_name in selected_method_scores) {
    message("  ", score_name, "_output: ", output_paths[[score_name]])
  }
  if (isTRUE(args$dry_run)) {
    next
  }

  if (needs_sample_mean) {
    sample_mean <- read_gene_prioritization_matrix(
      sample_mean_path,
      paste0("sample mean for method ", method)
    )
    if (!identical(colnames(sample_mean), selection$selected)) {
      stop(
        "Sample-mean cell types do not match selected method cell types for ",
        method
      )
    }
    specificity <- calculate_meanlog_specificity(sample_mean)
    for (score_name in intersect(
      selected_method_scores,
      c("meanlog_margin", "meanlog_mean_contrast")
    )) {
      write_gene_prioritization_matrix(
        specificity[[score_name]],
        output_paths[[score_name]],
        digits = args$digits
      )
    }
  }

  if (needs_sample_logmean) {
    sample_logmean <- read_gene_prioritization_matrix(
      sample_logmean_path,
      paste0("sample log-mean for method ", method)
    )
    if (!identical(colnames(sample_logmean), selection$selected)) {
      stop(
        "Sample-logmean cell types do not match selected method cell types for ",
        method
      )
    }
    specificity <- calculate_logmean_specificity(sample_logmean)
    for (score_name in intersect(
      selected_method_scores,
      c("logmean_margin", "logmean_contrast")
    )) {
      write_gene_prioritization_matrix(
        specificity[[score_name]],
        output_paths[[score_name]],
        digits = args$digits
      )
    }
  }

  if (needs_partial_r2) {
    prepared <- prepare_method_sample_stats_input(
      method_files,
      selection$selected,
      test_samples = test_info$samples,
      label = paste0("method ", method)
    )
    processed <- preprocess_method_sample_stats(
      prepared$matrices,
      preprocess_mode = preprocess_mode,
      scale_factor = scitd_config$scale_factor,
      lib_norm = args$lib_norm,
      zero_tolerance = scitd_config$zero_profile_tolerance
    )
    partial_r2 <- calculate_post_hoc_partial_r2(
      processed$y,
      bulk_prepared$bulk_expr,
      n_core = args$n_core,
      min_positive_samples = 2L
    )
    write_gene_prioritization_matrix(
      partial_r2,
      output_paths$partial_r2,
      digits = args$digits
    )
    message(
      "  wrote partial R2 for ",
      nrow(partial_r2),
      " genes x ",
      ncol(partial_r2),
      " cell types"
    )
  }
}

if (needs_genesigtest) {
  output_path <- gene_prioritization_genesigtest_output_path(paths)
  frac_path <- resolve_frac_input_path(paths, repo_root = repo_root)
  message("\nGeneSigTest_adapted")
  message("  frac_input: ", paths$config$frac_input[[1L]])
  message("  frac_path: ", frac_path)
  message("  output: ", output_path)
  if (!isTRUE(args$dry_run)) {
    frac <- read_frac_for_method(paths, repo_root = repo_root)
    restricted <- restrict_to_test_samples(
      bulk_prepared$bulk_expr,
      paths,
      frac = frac,
      use_test_samples = TRUE
    )
    bulk <- restricted$bulk
    frac <- validate_bulk_frac(bulk, restricted$frac)
    result <- gene_sig_test_adapted(bulk, frac)
    adjusted_fdr <- map_gene_prioritization_indep_columns(
      result$pval,
      paths,
      repo_root = repo_root
    )
    write_gene_prioritization_matrix(
      adjusted_fdr,
      output_path,
      digits = args$digits,
      digit_mode = "significant"
    )
    message(
      "  wrote adapted GeneSigTest FDR for ",
      nrow(adjusted_fdr),
      " genes x ",
      ncol(adjusted_fdr),
      " cell types using ",
      nrow(frac),
      " test samples"
    )
  }
}

message(if (isTRUE(args$dry_run)) "dry run complete" else "done")
