#!/usr/bin/env Rscript

args0 <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args0, value = TRUE)
script_path <- if (length(file_arg) > 0L) {
  normalizePath(sub("^--file=", "", file_arg[[1L]]), mustWork = TRUE)
} else {
  normalizePath(
    "DALE_Eval/eval/evaluate_truth_celltype_pairwise_cor.R",
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
source(file.path(
  repo_root,
  "DALE_Eval",
  "modules",
  "celltype_pairwise_cor_helpers.R"
))

usage <- function() {
  cat(
    "Usage:\n",
    "  Rscript DALE_Eval/eval/evaluate_truth_celltype_pairwise_cor.R \\\n",
    "    --dataset BRCA_Bassez2021 --truth_type meancpm \\\n",
    "    --digits 6 [--dry_run]\n\n",
    "Options:\n",
    "  --dataset     Benchmark dataset name. Required.\n",
    "  --truth_type  meancpm or sumcount_cpm. Required.\n",
    "  --digits      Decimal places written to disk. Default: 6.\n",
    "  --dry_run     Resolve and validate inputs/outputs without reading CTSE matrices or writing results.\n",
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
  out$digits <- suppressWarnings(as.integer(out$digits))
  if (is.na(out$digits) || out$digits < 0L) {
    stop("--digits must be a non-negative integer")
  }
  pairwise_truth_output_label(out$truth_type)
  out
}

args <- parse_cli_args(commandArgs(trailingOnly = TRUE))
obj_dir <- file.path(repo_root, "Benchmarking_obj", args$dataset)
if (!dir.exists(obj_dir)) {
  stop("Benchmark object not found: ", obj_dir)
}
truth_dir <- file.path(obj_dir, "ctse_truth", args$truth_type)
truth_files <- ctse_cell_type_files(truth_dir)
if (length(truth_files) == 0L) {
  stop("No truth CTSE .txt.gz files found under: ", truth_dir)
}
cell_types <- sort(names(truth_files))
validate_pairwise_cell_types(cell_types, label = "truth cell types")
test_info <- read_pairwise_test_samples(obj_dir)
out_path <- file.path(
  obj_dir,
  "deconv_performance",
  pairwise_truth_output_label(args$truth_type),
  "celltype_pairwise_pearson_cor",
  "Z_truth.txt"
)

message("Evaluating truth cell-type pairwise Pearson correlation")
message("  dataset: ", args$dataset)
message("  truth_type: ", args$truth_type)
message("  selected_cell_types: ", paste(cell_types, collapse = ", "))
message("  n_pairs: ", choose(length(cell_types), 2L))
message("  test_sample_path: ", test_info$path)
message("  n_test_samples_defined: ", length(test_info$samples))
message("  fraction_filter: none")
message("  sample_mask: none")
message("  digits: ", args$digits)
message("  output: ", out_path)
message("  dry_run: ", args$dry_run)

if (!isTRUE(args$dry_run)) {
  result <- compute_celltype_pairwise_pearson(
    files = truth_files,
    cell_types = cell_types,
    test_samples = test_info$samples,
    label = paste0("truth ", args$truth_type)
  )
  write_pairwise_cor_matrix(result$matrix, out_path, digits = args$digits)
  message(
    "wrote ", result$n_genes, " genes x ", nrow(result$pairs),
    " pairs using ", result$n_samples, " testing samples"
  )
}

message(if (isTRUE(args$dry_run)) "dry run complete" else "done")
