#!/usr/bin/env Rscript

# Add default noisy bulk inputs to Benchmarking_obj/<dataset>/bulk_input/.
#
# Two explicit modes are supported:
#   copy_noise_sweep_results=true:
#     Copy the four default noisy bulk inputs, plus matching deconv_res and
#     deconv_performance artifacts when present, from scripts/noise_sweep_test/.
#     Missing noisy bulk inputs stop the script.
#   copy_noise_sweep_results=false:
#     Regenerate all four default noisy bulk inputs from clean wcpm/sumcount,
#     overwriting existing noisy inputs and writing bulk_input/noise_qc/noisy_bulk_qc.pdf.

args0 <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args0, value = TRUE)
script_path <- if (length(file_arg) > 0) {
  normalizePath(sub("^--file=", "", file_arg[[1]]), mustWork = TRUE)
} else {
  normalizePath("scripts/step4_add_bulk_noise.R", mustWork = FALSE)
}
script_dir <- dirname(script_path)
repo_root_guess <- normalizePath(file.path(script_dir, ".."), mustWork = TRUE)

source(file.path(repo_root_guess, "DALE_Eval", "modules", "config_helpers.R"))
repo_root <- find_repo_root(repo_root_guess)
source(file.path(repo_root, "DALE_Eval", "modules", "pseudobulk.R"))

parse_bool <- function(x, label) {
  x <- tolower(trimws(x))
  if (!x %in% c("true", "false")) {
    stop(label, " must be true or false")
  }
  identical(x, "true")
}

parse_cli_args <- function(args) {
  out <- list(dataset = NULL, seed = 123L, copy_noise_sweep_results = FALSE)
  i <- 1L
  while (i <= length(args)) {
    key <- args[[i]]
    if (key %in% c("--help", "-h")) {
      cat(
        "Usage:\n",
        "  Rscript scripts/step4_add_bulk_noise.R [--dataset DATASET] [--seed 123] [--copy_noise_sweep_results false]\n\n",
        "Behavior:\n",
        "  copy_noise_sweep_results=true:\n",
        "    Copy the four default noisy bulk inputs from scripts/noise_sweep_test/<dataset>/noise_inputs/.\n",
        "    Also copy matching deconv_res and deconv_performance artifacts when present.\n",
        "    Missing expected noisy bulk inputs stop the script.\n",
        "  copy_noise_sweep_results=false:\n",
        "    Regenerate all four noisy bulk inputs from clean wcpm/sumcount, overwriting existing noisy inputs.\n",
        "    Write bulk_input/noise_qc/noisy_bulk_qc.pdf.\n\n",
        "If --dataset is omitted, all Benchmarking_obj datasets are considered.\n",
        sep = ""
      )
      quit(status = 0)
    }
    if (!startsWith(key, "--")) {
      stop("Unexpected positional argument: ", key)
    }
    if (i == length(args) || startsWith(args[[i + 1L]], "--")) {
      stop("Missing value for argument: ", key)
    }
    value <- args[[i + 1L]]
    key <- sub("^--", "", key)
    if (key == "dataset") {
      out$dataset <- value
    } else if (key == "seed") {
      out$seed <- as.integer(value)
      if (is.na(out$seed)) stop("--seed must be an integer")
    } else if (key == "copy_noise_sweep_results") {
      out$copy_noise_sweep_results <- parse_bool(value, "--copy_noise_sweep_results")
    } else {
      stop("Unsupported argument: --", key)
    }
    i <- i + 2L
  }
  out
}

read_matrix <- function(path) {
  x <- read.delim(path, sep = "\t", check.names = FALSE, row.names = 1)
  x <- as.matrix(x)
  storage.mode(x) <- "numeric"
  x[is.na(x)] <- 0
  x
}

write_matrix <- function(x, path, digits = 2) {
  if (file.exists(path)) {
    warning("overwriting existing noisy bulk input: ", path, immediate. = TRUE)
  }
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  idx <- is.finite(x) & x != 0
  x[idx] <- round(x[idx], digits)
  write.table(
    x,
    file = path,
    sep = "\t",
    quote = FALSE,
    col.names = NA
  )
  message("  generated: ", path)
  invisible(TRUE)
}

copy_file_verbose <- function(from, to) {
  if (!file.exists(from)) {
    stop("Missing source file: ", from)
  }
  if (file.exists(to)) {
    warning("overwriting existing target file: ", to, immediate. = TRUE)
  }
  dir.create(dirname(to), recursive = TRUE, showWarnings = FALSE)
  ok <- file.copy(from, to, overwrite = TRUE, copy.date = TRUE)
  if (!ok) stop("Failed to copy file from ", from, " to ", to)
  message("  copied: ", from, " -> ", to)
  invisible(TRUE)
}

copy_dir_verbose <- function(from, to) {
  if (!dir.exists(from)) {
    message("  missing source directory; skipped: ", from)
    return(invisible(FALSE))
  }
  if (dir.exists(to)) {
    warning("overwriting existing target directory: ", to, immediate. = TRUE)
    unlink(to, recursive = TRUE, force = TRUE)
  }
  dir.create(dirname(to), recursive = TRUE, showWarnings = FALSE)
  ok <- file.copy(from, dirname(to), recursive = TRUE, overwrite = TRUE, copy.date = TRUE)
  if (!ok) stop("Failed to copy directory from ", from, " to ", to)
  copied_path <- file.path(dirname(to), basename(from))
  if (!identical(normalizePath(copied_path, mustWork = FALSE), normalizePath(to, mustWork = FALSE))) {
    if (dir.exists(to)) unlink(to, recursive = TRUE, force = TRUE)
    ok <- file.rename(copied_path, to)
    if (!ok) stop("Failed to rename copied directory from ", copied_path, " to ", to)
  }
  message("  copied: ", from, " -> ", to)
  invisible(TRUE)
}

select_dataset_dirs <- function(repo_root, dataset = NULL) {
  benchmarking_dir <- file.path(repo_root, "Benchmarking_obj")
  if (!dir.exists(benchmarking_dir)) {
    stop("Benchmarking_obj directory not found: ", benchmarking_dir)
  }

  if (!is.null(dataset)) {
    dataset_dir <- file.path(benchmarking_dir, dataset)
    if (!dir.exists(dataset_dir)) {
      stop("Dataset directory not found: ", dataset_dir)
    }
    return(dataset_dir)
  }

  list.dirs(benchmarking_dir, recursive = FALSE, full.names = TRUE)
}

noise_sweep_dataset_dir <- function(repo_root, dataset) {
  file.path(repo_root, "scripts", "noise_sweep_test", dataset)
}

read_noise_sweep_configs <- function(repo_root, dataset) {
  noise_dir <- noise_sweep_dataset_dir(repo_root, dataset)
  shared_config <- file.path(repo_root, "scripts", "noise_sweep_test", "configs.tsv")
  dataset_config <- file.path(noise_dir, "configs.tsv")
  config_path <- if (file.exists(shared_config)) shared_config else dataset_config
  if (!file.exists(config_path)) {
    return(NULL)
  }
  read.delim(config_path, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
}

local_config_slug <- function(config_row) {
  paste0(
    config_row$config_id,
    "_bulk-", config_row$bulk_input,
    "__frac-", config_row$frac_input,
    "__ref-", config_row$refType
  )
}

find_noise_sweep_config <- function(repo_root, dataset, config_id) {
  configs <- read_noise_sweep_configs(repo_root, dataset)
  if (is.null(configs)) {
    return(NULL)
  }
  row <- configs[configs$config_id == config_id, , drop = FALSE]
  if (nrow(row) == 0) {
    return(NULL)
  }
  if (nrow(row) != 1) {
    stop("Expected one noise_sweep_test config for ", config_id, "; found ", nrow(row))
  }
  row
}

# Default noise settings from the settled centered/mu-shift benchmark design.
sigma_log2 <- 0.15
candidate_expression_threshold_cpm <- 1
positive_fraction <- 0.4
positive_mu_shift_mean_log2 <- -0.2
positive_mu_shift_sd_log2 <- 0.4
expressed_background_mu_shift <- 0
low_mu_shift_mean <- -1
low_mu_shift_sd <- 0.4

noise_map <- data.frame(
  target_name = c(
    "wcpm_centered_noise",
    "wcpm_mushift_noise",
    "sumcount_centered_noise",
    "sumcount_mushift_noise"
  ),
  source_config_id = c(
    "wcpm_centered_noise_sigma015",
    "wcpm_mushift_frac040_meanm020_sd040_sigma015",
    "sumcount_centered_noise_sigma015",
    "sumcount_mushift_frac040_meanm020_sd040_sigma015"
  ),
  target_config_id = c("config03", "config04", "config05", "config06"),
  input_family = c("wcpm", "wcpm", "sumcount", "sumcount"),
  expected_bulk_scale = c("cpm", "cpm", "counts", "counts"),
  truth_type = c("meancpm", "meancpm", "sumcount_cpm", "sumcount_cpm"),
  stringsAsFactors = FALSE
)

validate_target_configs <- function(repo_root) {
  for (i in seq_len(nrow(noise_map))) {
    row <- noise_map[i, , drop = FALSE]
    config <- read_deconv_config(row$target_config_id, repo_root = repo_root)
    expected <- list(
      bulk_input = row$target_name,
      bulk_scale = row$expected_bulk_scale,
      bulk_normalization = "cpm",
      frac_input = "InstaPrismfrac",
      refType = "indep"
    )
    for (field in names(expected)) {
      observed <- as.character(config[[field]][[1]])
      if (!identical(observed, expected[[field]])) {
        stop(
          "Target config ", row$target_config_id, " has unexpected ", field,
          ": observed=", observed, ", expected=", expected[[field]]
        )
      }
    }
  }
}

build_default_noise_outputs <- function(dataset_dir, seed) {
  bulk_input_dir <- file.path(dataset_dir, "bulk_input")
  wcpm_file <- file.path(bulk_input_dir, "wcpm.txt")
  sumcount_file <- file.path(bulk_input_dir, "sumcount.txt")

  if (!file.exists(wcpm_file)) stop("Missing clean wcpm input for fallback generation: ", wcpm_file)
  if (!file.exists(sumcount_file)) stop("Missing clean sumcount input for fallback generation: ", sumcount_file)

  set.seed(seed)

  wcpm_clean <- read_matrix(wcpm_file)
  sumcount_clean <- read_matrix(sumcount_file)
  sumcount_cpm_clean <- counts_to_cpm(sumcount_clean)

  wcpm_library_sizes <- matrix_col_sums(wcpm_clean)
  sumcount_library_sizes <- matrix_col_sums(sumcount_clean)

  wcpm_mu_shift_log2 <- make_sparse_expression_mu_shift(
    pseudobulk_cpm = wcpm_clean,
    candidate_expression_threshold_cpm = candidate_expression_threshold_cpm,
    positive_fraction = positive_fraction,
    positive_mu_shift_mean = positive_mu_shift_mean_log2,
    positive_mu_shift_sd = positive_mu_shift_sd_log2,
    expressed_background_mu_shift = expressed_background_mu_shift,
    low_mu_shift_mean = low_mu_shift_mean,
    low_mu_shift_sd = low_mu_shift_sd
  )

  sumcount_mu_shift_log2 <- make_sparse_expression_mu_shift(
    pseudobulk_cpm = sumcount_cpm_clean,
    candidate_expression_threshold_cpm = candidate_expression_threshold_cpm,
    positive_fraction = positive_fraction,
    positive_mu_shift_mean = positive_mu_shift_mean_log2,
    positive_mu_shift_sd = positive_mu_shift_sd_log2,
    expressed_background_mu_shift = expressed_background_mu_shift,
    low_mu_shift_mean = low_mu_shift_mean,
    low_mu_shift_sd = low_mu_shift_sd
  )

  list(
    wcpm_centered_noise = simulate_cpm_direct_multiplicative_noise(
      clean_cpm = wcpm_clean,
      mu_shift_log2 = 0,
      sigma = sigma_log2,
      renormalize_cpm = TRUE
    )$bulk,
    wcpm_mushift_noise = simulate_cpm_probability_multinomial_noise(
      clean_cpm = wcpm_clean,
      mu_shift_log2 = wcpm_mu_shift_log2,
      sigma = sigma_log2,
      library_sizes = wcpm_library_sizes,
      output_scale = "cpm"
    )$bulk,
    sumcount_centered_noise = simulate_count_probability_multinomial_noise(
      clean_counts = sumcount_clean,
      mu_shift_log2 = 0,
      sigma = sigma_log2,
      library_sizes = sumcount_library_sizes,
      count_sampler = "multinomial",
      output_scale = "counts"
    )$bulk,
    sumcount_mushift_noise = simulate_count_probability_multinomial_noise(
      clean_counts = sumcount_clean,
      mu_shift_log2 = sumcount_mu_shift_log2,
      sigma = sigma_log2,
      library_sizes = sumcount_library_sizes,
      count_sampler = "multinomial",
      output_scale = "counts"
    )$bulk
  )
}

write_default_noise_qc <- function(dataset_dir) {
  if (!requireNamespace("ggplot2", quietly = TRUE) ||
      !requireNamespace("gridExtra", quietly = TRUE)) {
    stop("Packages ggplot2 and gridExtra are required to write noisy bulk QC")
  }

  bulk_input_dir <- file.path(dataset_dir, "bulk_input")
  qc_dir <- file.path(bulk_input_dir, "noise_qc")
  qc_path <- file.path(qc_dir, "noisy_bulk_qc.pdf")

  if (file.exists(qc_path)) {
    warning("overwriting existing noisy bulk QC: ", qc_path, immediate. = TRUE)
  }

  wcpm_clean <- read_matrix(file.path(bulk_input_dir, "wcpm.txt"))
  sumcount_clean <- read_matrix(file.path(bulk_input_dir, "sumcount.txt"))

  qc_specs <- list(
    list(
      target_name = "wcpm_centered_noise",
      clean_bulk = wcpm_clean,
      clean_scale = "cpm",
      label = paste0("wcpm_centered_noise", " | sigma=", sigma_log2, " | mu_shift=0")
    ),
    list(
      target_name = "wcpm_mushift_noise",
      clean_bulk = wcpm_clean,
      clean_scale = "cpm",
      label = paste0(
        "wcpm_mushift_noise",
        " | sigma=", sigma_log2,
        " | fraction=", positive_fraction,
        " | mu_mean=", positive_mu_shift_mean_log2,
        " | mu_sd=", positive_mu_shift_sd_log2,
        " | threshold_cpm=", candidate_expression_threshold_cpm
      )
    ),
    list(
      target_name = "sumcount_centered_noise",
      clean_bulk = sumcount_clean,
      clean_scale = "counts",
      label = paste0("sumcount_centered_noise", " | sigma=", sigma_log2, " | mu_shift=0")
    ),
    list(
      target_name = "sumcount_mushift_noise",
      clean_bulk = sumcount_clean,
      clean_scale = "counts",
      label = paste0(
        "sumcount_mushift_noise",
        " | sigma=", sigma_log2,
        " | fraction=", positive_fraction,
        " | mu_mean=", positive_mu_shift_mean_log2,
        " | mu_sd=", positive_mu_shift_sd_log2,
        " | threshold_cpm=", candidate_expression_threshold_cpm
      )
    )
  )

  qc_grobs <- lapply(qc_specs, function(spec) {
    noisy_path <- file.path(bulk_input_dir, paste0(spec$target_name, ".txt"))
    if (!file.exists(noisy_path)) {
      stop("Missing noisy bulk input for QC: ", noisy_path)
    }
    plot_pseudobulk_noise_qc_grid(
      clean_bulk = spec$clean_bulk,
      simulated_obj = read_matrix(noisy_path),
      clean_scale = spec$clean_scale,
      label = spec$label
    )
  })

  dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)
  qc_plot <- gridExtra::arrangeGrob(grobs = qc_grobs, ncol = 1)
  ggplot2::ggsave(qc_path, plot = qc_plot, width = 12, height = 13.2)
  message("  wrote QC: ", qc_path)
  invisible(TRUE)
}

copy_noise_sweep_artifacts <- function(repo_root, dataset_dir) {
  dataset <- basename(dataset_dir)
  noise_dir <- noise_sweep_dataset_dir(repo_root, dataset)
  bulk_input_dir <- file.path(dataset_dir, "bulk_input")

  message("Copying default noisy bulk inputs and matching noise-sweep artifacts")
  message("  dataset: ", dataset)
  message("  source: ", noise_dir)

  source_paths <- file.path(
    noise_dir,
    "noise_inputs",
    paste0(noise_map$source_config_id, ".txt")
  )
  missing_sources <- source_paths[!file.exists(source_paths)]
  if (length(missing_sources) > 0) {
    stop(
      "copy_noise_sweep_results=true but expected noisy bulk inputs are missing:\n  ",
      paste(missing_sources, collapse = "\n  "),
      "\nRun with --copy_noise_sweep_results false to regenerate from clean bulk instead."
    )
  }

  for (i in seq_len(nrow(noise_map))) {
    row <- noise_map[i, , drop = FALSE]
    source_path <- file.path(noise_dir, "noise_inputs", paste0(row$source_config_id, ".txt"))
    target_path <- file.path(bulk_input_dir, paste0(row$target_name, ".txt"))
    copy_file_verbose(source_path, target_path)
  }

  for (i in seq_len(nrow(noise_map))) {
    row <- noise_map[i, , drop = FALSE]
    local_config <- find_noise_sweep_config(repo_root, dataset, row$source_config_id)
    if (is.null(local_config)) {
      message("  missing noise_sweep config; skipped artifacts for: ", row$source_config_id)
      next
    }

    source_slug <- local_config_slug(local_config)
    target_config <- read_deconv_config(row$target_config_id, repo_root = repo_root)
    target_slug <- config_slug(target_config)

    source_res_dir <- file.path(noise_dir, "deconv_res", source_slug)
    target_res_dir <- file.path(dataset_dir, "deconv_res", target_slug)
    copy_dir_verbose(source_res_dir, target_res_dir)

    source_perf_dir <- file.path(
      noise_dir,
      "deconv_performance",
      paste0(source_slug, "__truth-", row$truth_type)
    )
    target_perf_dir <- file.path(
      dataset_dir,
      "deconv_performance",
      paste0(target_slug, "__truth-", row$truth_type)
    )
    copy_dir_verbose(source_perf_dir, target_perf_dir)
  }

  invisible(TRUE)
}

regenerate_noise_inputs <- function(dataset_dir, seed) {
  dataset <- basename(dataset_dir)
  bulk_input_dir <- file.path(dataset_dir, "bulk_input")

  message("Regenerating default noisy bulk inputs")
  message("  dataset: ", dataset)
  message("  seed: ", seed)
  message("  source: clean wcpm.txt and sumcount.txt")

  generated_outputs <- build_default_noise_outputs(dataset_dir, seed = seed)
  for (i in seq_len(nrow(noise_map))) {
    row <- noise_map[i, , drop = FALSE]
    target_path <- file.path(bulk_input_dir, paste0(row$target_name, ".txt"))
    write_matrix(generated_outputs[[row$target_name]], target_path)
  }
  write_default_noise_qc(dataset_dir)
  invisible(TRUE)
}

args <- parse_cli_args(commandArgs(trailingOnly = TRUE))
validate_target_configs(repo_root)
dataset_dirs <- select_dataset_dirs(repo_root, dataset = args$dataset)

for (dataset_dir in dataset_dirs) {
  if (isTRUE(args$copy_noise_sweep_results)) {
    copy_noise_sweep_artifacts(repo_root, dataset_dir)
  } else {
    regenerate_noise_inputs(dataset_dir, seed = args$seed)
  }
}

message("Done")
