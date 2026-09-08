#!/usr/bin/env Rscript

args0 <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args0, value = TRUE)
script_path <- if (length(file_arg) > 0) {
  normalizePath(sub("^--file=", "", file_arg[[1]]), mustWork = TRUE)
} else {
  normalizePath("DALE_Eval/eval/evaluate_ctse_sample_cor.R", mustWork = FALSE)
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
      stop("Could not find repo root from: ", start)
    }
    cur <- parent
  }
}
repo_root <- find_repo_root_local(dirname(script_path))

source(file.path(repo_root, "DALE_Eval", "modules", "config_helpers.R"))
source(file.path(repo_root, "DALE_Eval", "modules", "evalu.R"))

usage <- function() {
  cat(
    "Usage:\n",
    "  Rscript ../DALE_Eval/eval/evaluate_ctse_sample_cor.R \\\n",
    "    --dataset LUAD_Kim2020 --config_id config01 --truth_type meancpm \\\n",
    "    --methods all --digits 6\n\n",
    "Computes within-sample, within-cell-type Spearman correlation across genes.\n",
    "All truth/method-overlapping samples are reported without fraction filtering.\n\n",
    "Options:\n",
    "  --dataset             Benchmark dataset name. Required.\n",
    "  --config_id           Deconvolution config id, e.g. config01. Required.\n",
    "  --truth_type          CTSE truth directory under ctse_truth/. Required.\n",
    "  --methods             Comma-separated methods, or all. Default: all.\n",
    "  --truth_transform     none or log2p1. Default: none.\n",
    "  --estimate_transform  none or log2p1. Default: none.\n",
    "  --hallmark_gmt        Hallmark GMT path, absolute or relative to repo root.\n",
    "                        Default: other_source_data/h.all.v7.5.1.symbols.gmt.\n",
    "  --digits              Correlation decimal places written to disk. Default: 6.\n",
    "  --dry_run             Show selected inputs and output paths without writing.\n",
    sep = ""
  )
}

split_csv <- function(x) {
  x <- trimws(unlist(strsplit(x, ",", fixed = TRUE)))
  x[nzchar(x)]
}

parse_bool <- function(x, label) {
  if (is.logical(x) && length(x) == 1L && !is.na(x)) {
    return(x)
  }
  value <- tolower(trimws(as.character(x)))
  if (!value %in% c("true", "false")) {
    stop(label, " must be true or false")
  }
  identical(value, "true")
}

parse_cli_args <- function(args) {
  out <- list(
    methods = "all",
    truth_transform = "none",
    estimate_transform = "none",
    hallmark_gmt = file.path(
      "other_source_data",
      "h.all.v7.5.1.symbols.gmt"
    ),
    digits = "6",
    dry_run = FALSE
  )
  allowed <- c(
    "dataset",
    "config_id",
    "truth_type",
    "methods",
    "truth_transform",
    "estimate_transform",
    "hallmark_gmt",
    "digits",
    "dry_run"
  )

  i <- 1L
  while (i <= length(args)) {
    key_raw <- args[[i]]
    if (key_raw %in% c("--help", "-h")) {
      usage()
      quit(status = 0)
    }
    if (!startsWith(key_raw, "--")) {
      stop("Unexpected positional argument: ", key_raw)
    }

    key <- gsub("-", "_", sub("^--", "", key_raw))
    if (!key %in% allowed) {
      stop("Unsupported argument: ", key_raw)
    }

    if (i == length(args) || startsWith(args[[i + 1L]], "--")) {
      out[[key]] <- TRUE
      i <- i + 1L
    } else {
      out[[key]] <- args[[i + 1L]]
      i <- i + 2L
    }
  }

  required <- c("dataset", "config_id", "truth_type")
  missing <- required[!vapply(
    required,
    function(x) !is.null(out[[x]]) && nzchar(as.character(out[[x]])),
    logical(1)
  )]
  if (length(missing) > 0) {
    usage()
    stop("Missing required arguments: ", paste(missing, collapse = ", "))
  }

  out$digits <- suppressWarnings(as.integer(out$digits))
  if (is.na(out$digits) || out$digits < 0L) {
    stop("--digits must be a non-negative integer")
  }
  out$dry_run <- parse_bool(out$dry_run, "--dry_run")
  out
}

resolve_repo_path <- function(path, repo_root) {
  if (!is.character(path) || length(path) != 1L || !nzchar(path)) {
    stop("--hallmark_gmt must be a non-empty path")
  }
  candidate <- if (startsWith(path, "/")) path else file.path(repo_root, path)
  normalizePath(candidate, mustWork = TRUE)
}

validate_transform <- function(value, label) {
  if (!value %in% c("none", "log2p1")) {
    stop(label, " must be none or log2p1")
  }
}

write_sample_cor_metadata <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  if (file.exists(path)) {
    warning("Overwriting existing sample-correlation metadata file: ", path)
  }
  write.table(
    x,
    file = path,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    na = "NA"
  )
}

evaluate_method_sample_cor <- function(truth_files,
                                       method_files,
                                       hallmark_genes,
                                       truth_transform,
                                       estimate_transform) {
  cell_types <- sort(intersect(names(truth_files), names(method_files)))
  if (length(cell_types) == 0) {
    stop("No overlapping cell types between truth and method CTSE")
  }

  all_vectors <- list()
  hallmark_vectors <- list()
  metadata <- NULL

  for (cell_type in cell_types) {
    truth <- read_ctse_matrix(
      truth_files[[cell_type]],
      paste0("truth ", cell_type)
    )
    estimate <- read_ctse_matrix(
      method_files[[cell_type]],
      paste0("estimate ", cell_type)
    )
    truth <- apply_ctse_transform(truth, truth_transform)
    estimate <- apply_ctse_transform(estimate, estimate_transform)

    overlapping_genes <- intersect(rownames(truth), rownames(estimate))
    hallmark_genes_used <- intersect(overlapping_genes, hallmark_genes)
    common_samples <- intersect(colnames(truth), colnames(estimate))

    all_values <- compute_ctse_sample_cor(
      truth = truth,
      estimate = estimate,
      genes = overlapping_genes,
      samples = common_samples
    )
    hallmark_values <- compute_ctse_sample_cor(
      truth = truth,
      estimate = estimate,
      genes = hallmark_genes_used,
      samples = common_samples
    )

    all_vectors[[cell_type]] <- all_values
    hallmark_vectors[[cell_type]] <- hallmark_values

    message(
      "    ", cell_type,
      ": truth_genes=", nrow(truth),
      ", inferred_genes=", nrow(estimate),
      ", overlapping_genes=", length(overlapping_genes),
      ", hallmark_genes=", length(hallmark_genes_used),
      ", common_samples=", length(common_samples)
    )

    if (is.null(metadata)) {
      metadata <- data.frame(
        gene_scope = c("all_genes", "hallmark_genes"),
        n_truth_genes = rep(nrow(truth), 2L),
        n_inferred_genes = rep(nrow(estimate), 2L),
        n_overlapping_genes = rep(length(overlapping_genes), 2L),
        n_genes_used = c(
          length(overlapping_genes),
          length(hallmark_genes_used)
        ),
        n_common_samples = rep(length(common_samples), 2L),
        stringsAsFactors = FALSE
      )
    }
  }

  list(
    all_genes = merge_sample_metric_vectors(all_vectors),
    hallmark_genes = merge_sample_metric_vectors(hallmark_vectors),
    metadata = metadata
  )
}

args <- parse_cli_args(commandArgs(trailingOnly = TRUE))
validate_transform(args$truth_transform, "--truth_transform")
validate_transform(args$estimate_transform, "--estimate_transform")

hallmark_path <- resolve_repo_path(args$hallmark_gmt, repo_root)
hallmark_genes <- read_gmt_gene_union(hallmark_path)
if (length(hallmark_genes) == 0) {
  stop("No Hallmark genes found in: ", hallmark_path)
}

config_row <- read_deconv_config(args$config_id, repo_root = repo_root)
slug <- config_slug(config_row)
paths <- resolve_deconv_paths(args$dataset, args$config_id, repo_root = repo_root)
truth_dir <- file.path(paths$obj_dir, "ctse_truth", args$truth_type)
truth_files <- ctse_cell_type_files(truth_dir)
if (length(truth_files) == 0) {
  stop("No CTSE truth files found under: ", truth_dir)
}

result_config_dir <- file.path(paths$obj_dir, "deconv_res", slug)
if (!dir.exists(result_config_dir)) {
  stop("Deconvolution result config directory not found: ", result_config_dir)
}

available_methods <- list.dirs(
  result_config_dir,
  recursive = FALSE,
  full.names = FALSE
)
if (identical(args$methods, "all")) {
  methods <- available_methods
} else {
  methods <- split_csv(args$methods)
  missing_methods <- setdiff(methods, available_methods)
  if (length(missing_methods) > 0) {
    stop(
      "Requested methods not found for config: ",
      paste(missing_methods, collapse = ", ")
    )
  }
}
methods <- sort(methods)
if (length(methods) == 0) {
  stop("No methods selected for evaluation")
}

perf_slug <- paste0(slug, "__truth-", args$truth_type)
perf_dir <- file.path(paths$obj_dir, "deconv_performance", perf_slug)
metadata_path <- file.path(perf_dir, "sample_cor_metadata.txt")
metadata_results <- list()

message("Evaluating CTSE within-sample cross-gene Spearman correlation")
message("  dataset: ", args$dataset)
message("  config_id: ", args$config_id)
message("  config_slug: ", slug)
message("  truth_type: ", args$truth_type)
message("  truth_transform: ", args$truth_transform)
message("  estimate_transform: ", args$estimate_transform)
message("  sample_policy: all truth/method-overlapping samples")
message("  fraction_filter: none")
message("  hallmark_gmt: ", hallmark_path)
message("  hallmark_gene_universe: ", length(hallmark_genes))
message("  digits: ", args$digits)
message("  methods: ", paste(methods, collapse = ", "))
message("  output_dir: ", perf_dir)
message("  metadata_output: ", metadata_path)
message("  dry_run: ", args$dry_run)

for (method in methods) {
  method_dir <- file.path(result_config_dir, method)
  method_files <- ctse_cell_type_files(method_dir)
  if (length(method_files) == 0) {
    warning("Skipping method with no CTSE .txt.gz files: ", method)
    next
  }

  message("\nMethod: ", method)
  overlap <- sort(intersect(names(truth_files), names(method_files)))
  extra_truth <- sort(setdiff(names(truth_files), names(method_files)))
  extra_method <- sort(setdiff(names(method_files), names(truth_files)))
  message("  overlapping_cell_types: ", paste(overlap, collapse = ", "))
  if (length(extra_truth) > 0) {
    message("  truth_only_cell_types: ", paste(extra_truth, collapse = ", "))
  }
  if (length(extra_method) > 0) {
    message("  method_only_cell_types: ", paste(extra_method, collapse = ", "))
  }
  if (length(overlap) == 0) {
    warning("Skipping method with no truth-overlapping cell types: ", method)
    next
  }

  all_path <- file.path(
    perf_dir,
    "sample_cor_all_genes",
    paste0(method, ".txt")
  )
  hallmark_output_path <- file.path(
    perf_dir,
    "sample_cor_hallmark_genes",
    paste0(method, ".txt")
  )
  message("  sample_cor_all_genes_output: ", all_path)
  message("  sample_cor_hallmark_genes_output: ", hallmark_output_path)
  message("  metadata_counts_from_cell_type: ", overlap[[1]])

  if (args$dry_run) {
    next
  }

  result <- evaluate_method_sample_cor(
    truth_files = truth_files,
    method_files = method_files,
    hallmark_genes = hallmark_genes,
    truth_transform = args$truth_transform,
    estimate_transform = args$estimate_transform
  )

  write_metric_matrix(result$all_genes, all_path, digits = args$digits)
  write_metric_matrix(
    result$hallmark_genes,
    hallmark_output_path,
    digits = args$digits
  )
  metadata_results[[method]] <- data.frame(
    method = rep(method, nrow(result$metadata)),
    result$metadata,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  message(
    "  wrote all-gene correlations: ",
    nrow(result$all_genes), " samples x ",
    ncol(result$all_genes), " cell types"
  )
  message(
    "  wrote Hallmark correlations: ",
    nrow(result$hallmark_genes), " samples x ",
    ncol(result$hallmark_genes), " cell types"
  )
}

if (!args$dry_run && length(metadata_results) > 0) {
  metadata <- do.call(rbind, unname(metadata_results))
  rownames(metadata) <- NULL
  write_sample_cor_metadata(metadata, metadata_path)
  message(
    "wrote config-level metadata: ", metadata_path,
    " (", nrow(metadata), " method/scope rows)"
  )
}

message(if (args$dry_run) "dry run done" else "done")
