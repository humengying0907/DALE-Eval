#!/usr/bin/env Rscript

args0 <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args0, value = TRUE)
script_path <- if (length(file_arg) > 0L) {
  normalizePath(sub("^--file=", "", file_arg[[1L]]), mustWork = TRUE)
} else {
  normalizePath(
    "DALE_Eval/eval/evaluate_method_celltype_pairwise_cor.R",
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
source(file.path(
  repo_root,
  "DALE_Eval",
  "modules",
  "celltype_pairwise_cor_helpers.R"
))

usage <- function() {
  cat(
    "Usage:\n",
    "  Rscript DALE_Eval/eval/evaluate_method_celltype_pairwise_cor.R \\\n",
    "    --dataset BRCA_Bassez2021 --config_id config01 \\\n",
    "    --methods all --digits 6 [--dry_run]\n\n",
    "Options:\n",
    "  --dataset     Benchmark dataset name. Required.\n",
    "  --config_id   Config present in both deconv_configs.txt and eval_configs.txt. Required.\n",
    "  --methods     Comma-separated method names, or all. Default: all.\n",
    "  --digits      Decimal places written to disk. Default: 6.\n",
    "  --dry_run     Resolve and validate inputs/outputs without reading CTSE matrices or writing results.\n",
    "  --help, -h    Show this help.\n",
    sep = ""
  )
}

parse_cli_args <- function(args) {
  out <- list(methods = "all", digits = "6", dry_run = FALSE)
  allowed <- c("dataset", "config_id", "methods", "digits", "dry_run")
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
    if (key == "dry_run") {
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
result_config_dir <- paths$output_dir
if (!dir.exists(result_config_dir)) {
  stop("Deconvolution result directory not found: ", result_config_dir)
}

available_methods <- sort(list.dirs(
  result_config_dir,
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

test_info <- read_pairwise_test_samples(paths$obj_dir)
performance_slug <- paste0(
  config_slug(paths$config),
  "__truth-",
  eval_config$expected_truth_type[[1L]]
)
metric_dir <- file.path(
  paths$obj_dir,
  "deconv_performance",
  performance_slug,
  "celltype_pairwise_pearson_cor"
)

message("Evaluating within-method cell-type pairwise Pearson correlation")
message("  dataset: ", args$dataset)
message("  config_id: ", args$config_id)
message("  config_slug: ", config_slug(paths$config))
message("  refType: ", paths$refType)
message("  methods: ", paste(methods, collapse = ", "))
message("  test_sample_path: ", test_info$path)
message("  n_test_samples_defined: ", length(test_info$samples))
message("  fraction_filter: none")
message("  sample_mask: none")
message("  digits: ", args$digits)
message("  output_dir: ", metric_dir)
message("  dry_run: ", args$dry_run)

for (method in methods) {
  method_dir <- file.path(result_config_dir, method)
  method_files <- ctse_cell_type_files(method_dir)
  if (length(method_files) == 0L) {
    stop("No CTSE .txt.gz files found for method: ", method)
  }
  selection <- resolve_pairwise_method_cell_types(
    paths,
    method_files,
    repo_root = repo_root
  )
  out_path <- file.path(metric_dir, paste0(method, ".txt"))

  message("\nMethod: ", method)
  message("  selected_cell_types: ", paste(selection$selected, collapse = ", "))
  if (length(selection$extra) > 0L) {
    message("  ignored_unmatched_cell_types: ", paste(sort(selection$extra), collapse = ", "))
  }
  message("  n_pairs: ", choose(length(selection$selected), 2L))
  message("  output: ", out_path)

  if (isTRUE(args$dry_run)) {
    next
  }
  result <- compute_celltype_pairwise_pearson(
    files = method_files,
    cell_types = selection$selected,
    test_samples = test_info$samples,
    label = paste0("method ", method)
  )
  write_pairwise_cor_matrix(result$matrix, out_path, digits = args$digits)
  message(
    "  wrote ", result$n_genes, " genes x ", nrow(result$pairs),
    " pairs using ", result$n_samples, " testing samples"
  )
}

message(if (isTRUE(args$dry_run)) "dry run complete" else "done")
