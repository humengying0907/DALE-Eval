#!/usr/bin/env Rscript

args0 <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args0, value = TRUE)
script_path <- if (length(file_arg) > 0) {
  normalizePath(sub("^--file=", "", file_arg[[1]]), mustWork = TRUE)
} else {
  normalizePath("DALE_Eval/eval/evaluate_ctse_performance.R", mustWork = FALSE)
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
    "  Rscript ../DALE_Eval/eval/evaluate_ctse_performance.R \\\n",
    "    --dataset BRCA_Bassez2021 --config_id config01 --truth_type meancpm \\\n",
    "    --metrics spearman_cor --filter_frac truth_cellfrac --min_frac 0.001\n\n",
    "Options:\n",
    "  --dataset             Benchmark dataset name. Required.\n",
    "  --config_id           Deconvolution config id, e.g. config01. Required.\n",
    "  --truth_type          CTSE truth directory under ctse_truth/. Required.\n",
    "  --methods             Comma-separated methods, or all. Default: all.\n",
    "  --metrics             Comma-separated: spearman_cor,pearson_cor. Default: spearman_cor.\n",
    "  --truth_transform     none or log2p1. Default: none.\n",
    "  --estimate_transform  none or log2p1. Default: none.\n",
    "  --filter_frac         none, truth_cellfrac, truth_transcriptfrac, method:InstaPrismfrac,\n",
    "                        or method:<method>:<file>. Default: truth_cellfrac.\n",
    "  --min_frac            Sample filter threshold for filter_frac. Default: 0.001.\n",
    "  --sample_mask         Optional sample-by-cell-type 0/1 mask path. Default: none.\n",
    "  --min_n_sample        Minimum samples per cell type. Default: 10.\n",
    "  --output_tag          Optional safe suffix for metric output folders. Default: empty.\n",
    "  --digits              Metric decimal places written to disk. Default: 2.\n",
    sep = ""
  )
}

parse_cli_args <- function(args) {
  out <- list(
    methods = "all",
    metrics = "spearman_cor",
    truth_transform = "none",
    estimate_transform = "none",
    filter_frac = "truth_cellfrac",
    min_frac = "0.001",
    sample_mask = "none",
    min_n_sample = "10",
    output_tag = "",
    digits = "2"
  )

  i <- 1L
  while (i <= length(args)) {
    key <- args[[i]]
    if (key %in% c("--help", "-h")) {
      usage()
      quit(status = 0)
    }
    if (!startsWith(key, "--")) {
      stop("Unexpected positional argument: ", key)
    }
    key <- gsub("-", "_", sub("^--", "", key))
    if (i == length(args) || startsWith(args[[i + 1L]], "--")) {
      out[[key]] <- TRUE
      i <- i + 1L
    } else {
      out[[key]] <- args[[i + 1L]]
      i <- i + 2L
    }
  }

  required <- c("dataset", "config_id", "truth_type")
  missing <- required[!vapply(required, function(x) !is.null(out[[x]]) && nzchar(out[[x]]), logical(1))]
  if (length(missing) > 0) {
    usage()
    stop("Missing required arguments: ", paste(missing, collapse = ", "))
  }

  out$min_frac <- as.numeric(out$min_frac)
  out$min_n_sample <- as.integer(out$min_n_sample)
  out$digits <- as.integer(out$digits)
  if (!is.finite(out$min_frac)) stop("--min_frac must be numeric")
  if (is.na(out$min_n_sample) || out$min_n_sample < 1) stop("--min_n_sample must be a positive integer")
  if (is.na(out$digits) || out$digits < 0) stop("--digits must be a non-negative integer")
  if (!is.character(out$sample_mask) || length(out$sample_mask) != 1L) {
    stop("--sample_mask must be a path or none")
  }
  if (!is.character(out$output_tag) || length(out$output_tag) != 1L) {
    stop("--output_tag must be a single safe identifier")
  }
  if (
    nzchar(out$output_tag) &&
      !grepl("^[A-Za-z0-9][A-Za-z0-9_-]*$", out$output_tag)
  ) {
    stop(
      "--output_tag must start with a letter or number and contain only ",
      "letters, numbers, underscores, or hyphens"
    )
  }

  out
}

split_csv <- function(x) {
  x <- trimws(unlist(strsplit(x, ",", fixed = TRUE)))
  x[nzchar(x)]
}

read_filter_matrix <- function(spec, paths) {
  if (identical(spec, "none")) {
    return(NULL)
  }

  if (spec %in% c("truth_cellfrac", "truth_transcriptfrac")) {
    path <- file.path(paths$obj_dir, "frac_input", paste0(spec, ".txt"))
    label <- spec
  } else if (startsWith(spec, "method:")) {
    token <- sub("^method:", "", spec)
    parts <- strsplit(token, ":", fixed = TRUE)[[1]]
    if (identical(token, "InstaPrismfrac")) {
      method <- "InstaPrism"
      filename <- "InstaPrismfrac.txt"
    } else if (length(parts) == 2) {
      method <- parts[[1]]
      filename <- parts[[2]]
      if (!grepl("\\.txt$", filename)) {
        filename <- paste0(filename, ".txt")
      }
    } else {
      stop(
        "Unsupported method filter_frac syntax: ", spec,
        ". Use method:InstaPrismfrac or method:<method>:<file>."
      )
    }
    path <- file.path(paths$output_dir, method, filename)
    label <- paste0("method fraction ", spec)
  } else if (file.exists(spec)) {
    path <- spec
    label <- spec
  } else {
    stop("Unsupported --filter_frac value: ", spec)
  }

  frac <- read_ctse_matrix(path, label = label)
  list(path = path, matrix = frac)
}

validate_choices <- function(values, allowed, label) {
  bad <- setdiff(values, allowed)
  if (length(bad) > 0) {
    stop(label, " contains unsupported values: ", paste(bad, collapse = ", "))
  }
}

read_sample_mask <- function(spec) {
  if (identical(spec, "none")) {
    return(NULL)
  }
  if (!file.exists(spec)) {
    stop("Sample mask not found: ", spec)
  }

  mask <- read.delim(
    spec,
    sep = "\t",
    check.names = FALSE,
    row.names = 1
  )
  mask <- as.matrix(mask)
  suppressWarnings(storage.mode(mask) <- "numeric")

  valid <- is.na(mask) | mask == 0 | mask == 1
  if (any(!valid)) {
    stop("Sample mask must contain only 0, 1, or NA: ", spec)
  }
  if (anyDuplicated(rownames(mask)) || anyDuplicated(colnames(mask))) {
    stop("Sample mask has duplicated sample or cell-type names: ", spec)
  }

  mask[is.na(mask)] <- 0
  list(path = normalizePath(spec), matrix = mask)
}

evaluate_method_metric <- function(truth_files,
                                   method_files,
                                   filter_frac,
                                   sample_mask,
                                   metric,
                                   truth_transform,
                                   estimate_transform,
                                   min_frac,
                                   min_n_sample) {
  cell_types <- sort(intersect(names(truth_files), names(method_files)))
  if (!is.null(filter_frac)) {
    cell_types <- intersect(cell_types, colnames(filter_frac$matrix))
  }
  if (!is.null(sample_mask)) {
    cell_types <- intersect(cell_types, colnames(sample_mask$matrix))
  }
  if (length(cell_types) == 0) {
    stop(
      "No overlapping cell types after applying truth/method/filter/mask ",
      "intersections"
    )
  }

  metric_vectors <- list()
  diagnostics <- data.frame()

  for (cell_type in cell_types) {
    truth <- read_ctse_matrix(truth_files[[cell_type]], paste0("truth ", cell_type))
    estimate <- read_ctse_matrix(method_files[[cell_type]], paste0("estimate ", cell_type))
    truth <- apply_ctse_transform(truth, truth_transform)
    estimate <- apply_ctse_transform(estimate, estimate_transform)

    filter_samples <- NULL
    if (!is.null(filter_frac)) {
      filter_values <- filter_frac$matrix[, cell_type]
      filter_samples <- rownames(filter_frac$matrix)[is.finite(filter_values) & filter_values > min_frac]
    }

    mask_samples <- NULL
    if (!is.null(sample_mask)) {
      mask_values <- sample_mask$matrix[, cell_type]
      mask_samples <- rownames(sample_mask$matrix)[mask_values == 1]
    }

    common_genes <- intersect(rownames(truth), rownames(estimate))
    common_samples <- intersect(colnames(truth), colnames(estimate))
    eval_samples <- common_samples
    if (!is.null(filter_samples)) {
      eval_samples <- intersect(eval_samples, filter_samples)
    }
    if (!is.null(mask_samples)) {
      eval_samples <- intersect(eval_samples, mask_samples)
    }

    message(
      "    ", cell_type,
      ": genes=", length(common_genes),
      ", common_samples=", length(common_samples),
      if (!is.null(filter_samples)) {
        paste0(", filter_pass=", length(intersect(common_samples, filter_samples)))
      } else {
        ""
      },
      if (!is.null(mask_samples)) {
        paste0(", mask_pass=", length(intersect(common_samples, mask_samples)))
      } else {
        ""
      },
      ", eval_samples=", length(eval_samples)
    )

    metric_vectors[[cell_type]] <- compute_ctse_metric_rows(
      truth = truth,
      estimate = estimate,
      metric = metric,
      samples = eval_samples,
      min_n_sample = min_n_sample
    )

    diagnostics <- rbind(
      diagnostics,
      data.frame(
        cell_type = cell_type,
        n_genes = length(common_genes),
        n_common_samples = length(common_samples),
        n_eval_samples = length(eval_samples),
        stringsAsFactors = FALSE
      )
    )
  }

  list(
    matrix = merge_metric_vectors(metric_vectors),
    diagnostics = diagnostics
  )
}

args <- parse_cli_args(commandArgs(trailingOnly = TRUE))
metrics <- split_csv(args$metrics)
validate_choices(metrics, c("spearman_cor", "pearson_cor"), "--metrics")
validate_choices(args$truth_transform, c("none", "log2p1"), "--truth_transform")
validate_choices(args$estimate_transform, c("none", "log2p1"), "--estimate_transform")

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

available_methods <- list.dirs(result_config_dir, recursive = FALSE, full.names = FALSE)
if (identical(args$methods, "all")) {
  methods <- available_methods
} else {
  methods <- split_csv(args$methods)
  missing_methods <- setdiff(methods, available_methods)
  if (length(missing_methods) > 0) {
    stop("Requested methods not found for config: ", paste(missing_methods, collapse = ", "))
  }
}
methods <- sort(methods)
if (length(methods) == 0) {
  stop("No methods selected for evaluation")
}

filter_frac <- read_filter_matrix(args$filter_frac, paths)
sample_mask <- read_sample_mask(args$sample_mask)
perf_slug <- paste0(slug, "__truth-", args$truth_type)
perf_dir <- file.path(paths$obj_dir, "deconv_performance", perf_slug)

message("Evaluating CTSE performance")
message("  dataset: ", args$dataset)
message("  config_id: ", args$config_id)
message("  config_slug: ", slug)
message("  truth_type: ", args$truth_type)
message("  truth_transform: ", args$truth_transform)
message("  estimate_transform: ", args$estimate_transform)
message("  filter_frac: ", args$filter_frac)
if (!is.null(filter_frac)) {
  message("  filter_frac_path: ", filter_frac$path)
}
message("  min_frac: ", args$min_frac)
message("  sample_mask: ", args$sample_mask)
if (!is.null(sample_mask)) {
  message("  sample_mask_path: ", sample_mask$path)
}
message("  min_n_sample: ", args$min_n_sample)
message(
  "  output_tag: ",
  if (nzchar(args$output_tag)) args$output_tag else "<default>"
)
message("  digits: ", args$digits)
message("  methods: ", paste(methods, collapse = ", "))
message("  metrics: ", paste(metrics, collapse = ", "))
message("  output_dir: ", perf_dir)

for (method in methods) {
  method_dir <- file.path(result_config_dir, method)
  method_files <- ctse_cell_type_files(method_dir)
  if (length(method_files) == 0) {
    warning("Skipping method with no CTSE .txt.gz files: ", method)
    next
  }

  message("\nMethod: ", method)
  overlap <- sort(intersect(names(truth_files), names(method_files)))
  extra_truth <- setdiff(names(truth_files), names(method_files))
  extra_method <- setdiff(names(method_files), names(truth_files))
  message("  overlapping_cell_types: ", paste(overlap, collapse = ", "))
  if (length(extra_truth) > 0) {
    message("  truth_only_cell_types: ", paste(sort(extra_truth), collapse = ", "))
  }
  if (length(extra_method) > 0) {
    message("  method_only_cell_types: ", paste(sort(extra_method), collapse = ", "))
  }

  for (metric in metrics) {
    message("  Metric: ", metric)
    result <- evaluate_method_metric(
      truth_files = truth_files,
      method_files = method_files,
      filter_frac = filter_frac,
      sample_mask = sample_mask,
      metric = metric,
      truth_transform = args$truth_transform,
      estimate_transform = args$estimate_transform,
      min_frac = args$min_frac,
      min_n_sample = args$min_n_sample
    )

    metric_output_folder <- if (nzchar(args$output_tag)) {
      paste0(metric, "_", args$output_tag)
    } else {
      metric
    }
    out_path <- file.path(
      perf_dir,
      metric_output_folder,
      paste0(method, ".txt")
    )
    write_metric_matrix(result$matrix, out_path, digits = args$digits)
    message(
      "    wrote ", out_path,
      " (", nrow(result$matrix), " genes x ", ncol(result$matrix), " cell types)"
    )
  }
}

message("done")
