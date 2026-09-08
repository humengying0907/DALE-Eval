# Build the complete dataset-level bulk baseline correlation results for config01-config08.
#
# Notebook-style workflow:
#   1. Edit the settings block.
#   2. Source or run from repo root, benchmark_summary/, or its preparation/ folder.
#   3. Inspect messages and output files.
#
# Output:
#   Benchmarking_obj/<dataset>/
#     regressed_bulk/
#       bulk-<bulk_input>[__norm-cpm]__regress-<frac_source>.txt
#     deconv_performance/bulk_baseline/
#       spearman_cor/
#         truth-<truth_type>__bulk-<bulk_input>[__norm-cpm][__regress-<frac_source>].txt
#       manifest.tsv
#
# Baselines:
#   - raw/prepared bulk expression repeated as a per-cell-type estimate
#   - prepared bulk expression after regressing out the matching truth fraction
#   - prepared bulk expression after regressing out the matching config01-config08 InstaPrism fractions
#
# Correlations use declared test samples that pass the truth-fraction threshold.
# Regression fits retain all samples shared by bulk and the regression fractions;
# neither test-sample selection nor a high-fraction mask changes the fitted model.

## ---------------------------------------------------------------------------
## Settings
## ---------------------------------------------------------------------------

repo_root_candidates <- unique(normalizePath(
  c(".", "..", "../.."),
  mustWork = FALSE
))
repo_root_matches <- repo_root_candidates[
  file.exists(file.path(repo_root_candidates, "Benchmarking_obj"))
]
if (length(repo_root_matches) == 0L) {
  stop("Run from repo root, benchmark_summary/, or benchmark_summary/preparation/")
}
repo_root <- normalizePath(repo_root_matches[[1]])

dataset_include <- NULL

# The standard section uses all eligible test samples. The high-fraction
# section additionally applies its sample mask. Set either flag to FALSE
# to skip that section entirely.
run_all_sample_baselines <- TRUE
run_high_frac_baselines <- TRUE

# Edit this to select which reusable bulk baselines to build.
# Regression baselines write the regressed bulk matrices under
# Benchmarking_obj/<dataset>/regressed_bulk/ and their correlations under
# Benchmarking_obj/<dataset>/deconv_performance/bulk_baseline/spearman_cor/.
bulk_recipe_include <- c(
  "bulk-wcpm",
  "bulk-sumcount__norm-cpm",
  "bulk-wcpm_centered_noise",
  "bulk-wcpm_mushift_noise",
  "bulk-sumcount_centered_noise__norm-cpm",
  "bulk-sumcount_mushift_noise__norm-cpm"
)

metrics <- c("spearman_cor")
min_frac <- 0.001
min_n_sample <- 10L
digits <- 6L

output_folder_name <- "bulk_baseline"

bulk_recipes <- data.frame(
  recipe_id = c(
    "bulk-wcpm",
    "bulk-sumcount__norm-cpm",
    "bulk-wcpm_centered_noise",
    "bulk-wcpm_mushift_noise",
    "bulk-sumcount_centered_noise__norm-cpm",
    "bulk-sumcount_mushift_noise__norm-cpm"
  ),
  bulk_input = c(
    "wcpm",
    "sumcount",
    "wcpm_centered_noise",
    "wcpm_mushift_noise",
    "sumcount_centered_noise",
    "sumcount_mushift_noise"
  ),
  prep_config_id = c(
    "config01",
    "config02",
    "config03",
    "config04",
    "config05",
    "config06"
  ),
  truth_type = c(
    "meancpm",
    "sumcount_cpm",
    "meancpm",
    "meancpm",
    "sumcount_cpm",
    "sumcount_cpm"
  ),
  filter_frac = c(
    "truth_cellfrac",
    "truth_transcriptfrac",
    "truth_cellfrac",
    "truth_cellfrac",
    "truth_transcriptfrac",
    "truth_transcriptfrac"
  ),
  truth_regress_frac = c(
    "truth_cellfrac",
    "truth_transcriptfrac",
    "truth_cellfrac",
    "truth_cellfrac",
    "truth_transcriptfrac",
    "truth_transcriptfrac"
  ),
  instaprism_regress_config_ids = c(
    "config01,config07",
    "config02,config08",
    "config03",
    "config04",
    "config05",
    "config06"
  ),
  stringsAsFactors = FALSE
)

## ---------------------------------------------------------------------------
## Setup
## ---------------------------------------------------------------------------

source(file.path(repo_root, "DALE_Eval", "modules", "config_helpers.R"))
source(file.path(repo_root, "DALE_Eval", "modules", "runner_helpers.R"))
source(file.path(repo_root, "DALE_Eval", "modules", "evalu.R"))

validate_values <- function(values, allowed, label) {
  bad <- setdiff(values, allowed)
  if (length(bad) > 0) {
    stop(label, " contains unsupported values: ", paste(bad, collapse = ", "))
  }
}

validate_values(metrics, c("spearman_cor", "pearson_cor"), "metrics")

bulk_recipes <- bulk_recipes[bulk_recipes$recipe_id %in% bulk_recipe_include, , drop = FALSE]
missing_recipes <- setdiff(bulk_recipe_include, bulk_recipes$recipe_id)
if (length(missing_recipes) > 0) {
  stop("Unknown bulk_recipe_include values: ", paste(missing_recipes, collapse = ", "))
}
if (nrow(bulk_recipes) == 0) {
  stop("No bulk recipes selected")
}

benchmark_root <- file.path(repo_root, "Benchmarking_obj")
datasets <- list.dirs(benchmark_root, full.names = FALSE, recursive = FALSE)
datasets <- datasets[file.exists(file.path(benchmark_root, datasets, "deconv_performance"))]

if (!is.null(dataset_include)) {
  datasets <- intersect(datasets, dataset_include)
}

if (length(datasets) == 0) {
  stop("No datasets selected")
}

## ---------------------------------------------------------------------------
## Small Helpers
## ---------------------------------------------------------------------------

read_numeric_matrix <- function(path, label = "matrix") {
  if (!file.exists(path)) {
    stop(label, " not found: ", path)
  }
  x <- read.delim(path, sep = "\t", check.names = FALSE, row.names = 1)
  x <- as.matrix(x)
  storage.mode(x) <- "numeric"
  x[is.na(x)] <- 0
  x
}

write_numeric_matrix <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  write.table(
    x,
    file = path,
    sep = "\t",
    quote = FALSE,
    col.names = NA
  )
}

read_config_fraction <- function(dataset_root, config_id) {
  config_row <- read_deconv_config(config_id, repo_root = repo_root)
  slug <- config_slug(config_row)
  frac_path <- file.path(
    dataset_root,
    "deconv_res",
    slug,
    "InstaPrism",
    "InstaPrismfrac.txt"
  )
  list(path = frac_path, matrix = read_numeric_matrix(frac_path, "InstaPrism fraction"))
}

regress_out_fraction <- function(bulk_expr, frac) {
  bulk_expr <- as.matrix(bulk_expr)
  frac <- as.matrix(frac)

  common_samples <- intersect(colnames(bulk_expr), rownames(frac))
  if (length(common_samples) < 2L) {
    stop("Too few common samples for fraction regression")
  }

  y <- bulk_expr[, common_samples, drop = FALSE]
  mod <- frac[common_samples, , drop = FALSE]
  storage.mode(mod) <- "numeric"

  keep_cols <- colSums(is.finite(mod)) == nrow(mod) & colSums(abs(mod), na.rm = TRUE) > 0
  mod <- mod[, keep_cols, drop = FALSE]
  if (ncol(mod) == 0L) {
    stop("No usable fraction columns for regression")
  }

  mod_qr <- qr(mod)
  if (mod_qr$rank < ncol(mod)) {
    keep <- sort(mod_qr$pivot[seq_len(mod_qr$rank)])
    message("    fraction regression: dropping rank-dependent fraction columns")
    mod <- mod[, keep, drop = FALSE]
    mod_qr <- qr(mod)
  }

  beta <- qr.coef(mod_qr, t(y))
  fitted <- t(mod %*% beta)
  residuals <- y - fitted
  residuals + rowMeans(y)
}

evaluate_baseline_metric <- function(truth_files,
                                     estimate,
                                     filter_frac,
                                     metric,
                                     min_frac,
                                     min_n_sample,
                                     test_samples,
                                     sample_mask = NULL) {
  cell_types <- sort(intersect(names(truth_files), colnames(filter_frac)))
  if (!is.null(sample_mask)) {
    cell_types <- intersect(cell_types, colnames(sample_mask))
  }
  if (length(cell_types) == 0L) {
    stop("No overlapping cell types among truth, fraction filter, and sample mask")
  }

  metric_vectors <- list()

  for (cell_type in cell_types) {
    truth <- read_ctse_matrix(truth_files[[cell_type]], paste0("truth ", cell_type))
    filter_values <- filter_frac[, cell_type]
    filter_samples <- rownames(filter_frac)[is.finite(filter_values) & filter_values > min_frac]
    eval_samples <- intersect(filter_samples, test_samples)

    if (!is.null(sample_mask)) {
      mask_values <- sample_mask[, cell_type]
      mask_samples <- rownames(sample_mask)[is.finite(mask_values) & mask_values == 1]
      eval_samples <- intersect(eval_samples, mask_samples)
    }

    metric_vectors[[cell_type]] <- compute_ctse_metric_rows(
      truth = truth,
      estimate = estimate,
      metric = metric,
      samples = eval_samples,
      min_n_sample = min_n_sample
    )
  }

  merge_metric_vectors(metric_vectors)
}

write_manifest <- function(manifest_rows, path) {
  if (length(manifest_rows) == 0L) {
    return(invisible(NULL))
  }
  manifest <- do.call(rbind, manifest_rows)
  write.table(
    manifest,
    file = path,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
}

## ---------------------------------------------------------------------------
## Build High-Fraction Baseline Correlations
## ---------------------------------------------------------------------------

if (run_high_frac_baselines) {
  high_frac_recipes <- data.frame(
    recipe_id = c(
      "bulk-wcpm",
      "bulk-sumcount__norm-cpm"
    ),
    prep_config_id = c(
      "config01",
      "config02"
    ),
    truth_type = c(
      "meancpm",
      "sumcount_cpm"
    ),
    filter_frac = c(
      "truth_cellfrac",
      "truth_transcriptfrac"
    ),
    truth_regress_frac = c(
      "truth_cellfrac",
      "truth_transcriptfrac"
    ),
    sample_mask_file = c(
      "InstaPrismfrac_config01_gt0.1.txt",
      "InstaPrismfrac_config02_gt0.1.txt"
    ),
    stringsAsFactors = FALSE
  )

  for (dataset_name in datasets) {
    message("High-fraction dataset: ", dataset_name)

    dataset_root <- file.path(benchmark_root, dataset_name)
    out_dir <- file.path(
      dataset_root,
      "deconv_performance",
      output_folder_name,
      "spearman_cor_highFrac_sample"
    )

    for (recipe_idx in seq_len(nrow(high_frac_recipes))) {
      recipe <- high_frac_recipes[recipe_idx, , drop = FALSE]
      recipe_id <- recipe$recipe_id[[1]]
      config_id <- recipe$prep_config_id[[1]]
      truth_type <- recipe$truth_type[[1]]
      regress_frac_name <- recipe$truth_regress_frac[[1]]

      message("  High-fraction recipe: ", recipe_id)

      truth_dir <- file.path(dataset_root, "ctse_truth", truth_type)
      filter_frac_path <- file.path(
        dataset_root,
        "frac_input",
        paste0(recipe$filter_frac[[1]], ".txt")
      )
      sample_mask_path <- file.path(
        dataset_root,
        "evaluation_metadata",
        "sample_mask",
        recipe$sample_mask_file[[1]]
      )

      required_paths <- c(
        "truth directory" = truth_dir,
        "fraction filter" = filter_frac_path,
        "sample mask" = sample_mask_path
      )
      missing_required <- names(required_paths)[!file.exists(required_paths)]
      if (length(missing_required) > 0L) {
        message(
          "    skip: missing ",
          paste(missing_required, collapse = ", ")
        )
        next
      }

      paths <- tryCatch(
        resolve_deconv_paths(dataset_name, config_id, repo_root = repo_root),
        error = function(e) e
      )
      if (inherits(paths, "error")) {
        message("    skip: cannot resolve prep config: ", conditionMessage(paths))
        next
      }
      if (!file.exists(paths$bulk_path)) {
        message("    skip: missing bulk input ", paths$bulk_path)
        next
      }

      truth_paths <- list.files(
        truth_dir,
        pattern = "[.]txt([.]gz)?$",
        full.names = TRUE
      )
      if (length(truth_paths) == 0L) {
        message("    skip: no CTSE truth files in ", truth_dir)
        next
      }
      truth_files <- stats::setNames(
        truth_paths,
        sub("[.]txt([.]gz)?$", "", basename(truth_paths))
      )

      bulk_expr <- read_numeric_matrix(paths$bulk_path, "bulk input")
      bulk_prep <- prepare_bulk_for_deconv(
        bulk_expr,
        paths,
        method = "Bulk baseline"
      )
      bulk_estimate <- bulk_prep$bulk_expr
      filter_frac <- read_numeric_matrix(filter_frac_path, "fraction filter")
      sample_mask <- read_numeric_matrix(sample_mask_path, "sample mask")
      test_samples <- restrict_to_test_samples(
        bulk = bulk_estimate,
        paths = paths,
        frac = filter_frac,
        use_test_samples = TRUE
      )$selected_samples
      message("    evaluation: ", length(test_samples), " test samples before cell-type filtering")

      regressed_path <- file.path(
        dataset_root,
        "regressed_bulk",
        paste0(recipe_id, "__regress-", regress_frac_name, ".txt")
      )
      if (file.exists(regressed_path)) {
        regressed_estimate <- read_numeric_matrix(
          regressed_path,
          "truth-fraction-regressed bulk"
        )
      } else {
        message("    building missing regressed bulk: ", regressed_path)
        regress_frac_path <- file.path(
          dataset_root,
          "frac_input",
          paste0(regress_frac_name, ".txt")
        )
        regress_frac <- read_numeric_matrix(
          regress_frac_path,
          "truth regression fraction"
        )
        regressed_estimate <- regress_out_fraction(bulk_estimate, regress_frac)
        write_numeric_matrix(regressed_estimate, regressed_path)
      }

      estimates <- list(
        raw = bulk_estimate,
        truth_regressed = regressed_estimate
      )
      output_names <- c(
        raw = paste0("truth-", truth_type, "__", recipe_id, ".txt"),
        truth_regressed = paste0(
          "truth-", truth_type,
          "__", recipe_id,
          "__regress-", regress_frac_name,
          ".txt"
        )
      )

      for (estimate_name in names(estimates)) {
        metric_table <- evaluate_baseline_metric(
          truth_files = truth_files,
          estimate = estimates[[estimate_name]],
          filter_frac = filter_frac,
          metric = "spearman_cor",
          min_frac = min_frac,
          min_n_sample = min_n_sample,
          test_samples = test_samples,
          sample_mask = sample_mask
        )
        write_numeric_matrix(
          round(metric_table, digits = digits),
          file.path(out_dir, output_names[[estimate_name]])
        )
      }
    }
  }
}

## ---------------------------------------------------------------------------
## Build Standard Test-Sample Baseline Correlations
## ---------------------------------------------------------------------------

for (dataset_name in datasets) {
  if (!run_all_sample_baselines) {
    next
  }

  message("Dataset: ", dataset_name)

  dataset_root <- file.path(benchmark_root, dataset_name)
  out_root <- file.path(dataset_root, "deconv_performance", output_folder_name)
  dir.create(out_root, recursive = TRUE, showWarnings = FALSE)

  manifest_rows <- list()

  for (recipe_idx in seq_len(nrow(bulk_recipes))) {
    recipe <- bulk_recipes[recipe_idx, , drop = FALSE]
    recipe_id <- recipe$recipe_id[[1]]
    truth_type <- recipe$truth_type[[1]]

    message("  Recipe: ", recipe_id)

    truth_dir <- file.path(dataset_root, "ctse_truth", truth_type)
    filter_frac_path <- file.path(dataset_root, "frac_input", paste0(recipe$filter_frac[[1]], ".txt"))

    if (!dir.exists(truth_dir)) {
      message("    skip: missing truth dir ", truth_dir)
      next
    }
    if (!file.exists(filter_frac_path)) {
      message("    skip: missing filter fraction ", filter_frac_path)
      next
    }

    paths <- tryCatch(
      resolve_deconv_paths(dataset_name, recipe$prep_config_id[[1]], repo_root = repo_root),
      error = function(e) e
    )
    if (inherits(paths, "error")) {
      message("    skip: cannot resolve prep config: ", conditionMessage(paths))
      next
    }
    if (!file.exists(paths$bulk_path)) {
      message("    skip: missing bulk input ", paths$bulk_path)
      next
    }

    if (!identical(paths$config$bulk_input[[1]], recipe$bulk_input[[1]])) {
      stop(
        "Recipe ", recipe_id, " expected bulk_input=", recipe$bulk_input[[1]],
        " but prep_config_id=", recipe$prep_config_id[[1]],
        " resolves to bulk_input=", paths$config$bulk_input[[1]]
      )
    }

    truth_files <- ctse_cell_type_files(truth_dir)
    filter_frac <- read_numeric_matrix(filter_frac_path, "filter fraction")

    bulk_raw <- read_bulk(paths)
    bulk_prep <- prepare_bulk_for_deconv(bulk_raw, paths, method = "bulk baseline")
    test_samples <- restrict_to_test_samples(
      bulk = bulk_prep$bulk_expr,
      paths = paths,
      frac = filter_frac,
      use_test_samples = TRUE
    )$selected_samples
    message("    evaluation: ", length(test_samples), " test samples before cell-type filtering")
    baseline_matrices <- list()
    baseline_meta <- list()
    baseline_bulk_paths <- list()

    baseline_matrices[[recipe_id]] <- bulk_prep$bulk_expr
    baseline_bulk_paths[[recipe_id]] <- ""
    baseline_meta[[recipe_id]] <- list(
      regress_type = "none",
      regress_frac = "",
      regress_frac_path = ""
    )

    truth_regress_path <- file.path(
      dataset_root,
      "frac_input",
      paste0(recipe$truth_regress_frac[[1]], ".txt")
    )
    if (file.exists(truth_regress_path)) {
      truth_regress_frac <- read_numeric_matrix(truth_regress_path, "truth regression fraction")
      baseline_id <- paste0(recipe_id, "__regress-", recipe$truth_regress_frac[[1]])
      baseline_matrices[[baseline_id]] <- regress_out_fraction(bulk_prep$bulk_expr, truth_regress_frac)
      baseline_bulk_paths[[baseline_id]] <- file.path(
        dataset_root,
        "regressed_bulk",
        paste0(baseline_id, ".txt")
      )
      write_numeric_matrix(baseline_matrices[[baseline_id]], baseline_bulk_paths[[baseline_id]])
      message("    wrote regressed bulk: ", baseline_bulk_paths[[baseline_id]])
      baseline_meta[[baseline_id]] <- list(
        regress_type = "truth_fraction",
        regress_frac = recipe$truth_regress_frac[[1]],
        regress_frac_path = truth_regress_path
      )
    } else {
      message("    skip truth-fraction regression: missing ", truth_regress_path)
    }

    instaprism_config_ids <- strsplit(
      recipe$instaprism_regress_config_ids[[1]],
      ",",
      fixed = TRUE
    )[[1]]

    for (instaprism_config_id in instaprism_config_ids) {
      insta_frac <- tryCatch(
        read_config_fraction(dataset_root, instaprism_config_id),
        error = function(e) e
      )
      if (inherits(insta_frac, "error")) {
        message(
          "    skip InstaPrism regression for ",
          instaprism_config_id,
          ": ",
          conditionMessage(insta_frac)
        )
        next
      }

      baseline_id <- paste0(
        recipe_id,
        "__regress-InstaPrismfrac_",
        instaprism_config_id
      )
      baseline_matrices[[baseline_id]] <- regress_out_fraction(
        bulk_prep$bulk_expr,
        insta_frac$matrix
      )
      baseline_bulk_paths[[baseline_id]] <- file.path(
        dataset_root,
        "regressed_bulk",
        paste0(baseline_id, ".txt")
      )
      write_numeric_matrix(
        baseline_matrices[[baseline_id]],
        baseline_bulk_paths[[baseline_id]]
      )
      message("    wrote regressed bulk: ", baseline_bulk_paths[[baseline_id]])
      baseline_meta[[baseline_id]] <- list(
        regress_type = "InstaPrism_fraction",
        regress_frac = paste0("InstaPrismfrac_", instaprism_config_id),
        regress_frac_path = insta_frac$path
      )
    }

    for (baseline_id in names(baseline_matrices)) {
      estimate <- if (nzchar(baseline_bulk_paths[[baseline_id]])) {
        read_numeric_matrix(
          baseline_bulk_paths[[baseline_id]],
          paste0("regressed bulk baseline ", baseline_id)
        )
      } else {
        baseline_matrices[[baseline_id]]
      }
      meta <- baseline_meta[[baseline_id]]

      message("    Baseline: ", baseline_id)

      for (metric in metrics) {
        result <- evaluate_baseline_metric(
          truth_files = truth_files,
          estimate = estimate,
          filter_frac = filter_frac,
          metric = metric,
          min_frac = min_frac,
          min_n_sample = min_n_sample,
          test_samples = test_samples
        )

        out_path <- file.path(
          out_root,
          metric,
          paste0("truth-", truth_type, "__", baseline_id, ".txt")
        )
        write_metric_matrix(result, out_path, digits = digits)

        manifest_rows[[length(manifest_rows) + 1L]] <- data.frame(
          dataset = dataset_name,
          recipe_id = recipe_id,
          baseline_id = baseline_id,
          bulk_input = recipe$bulk_input[[1]],
          prep_config_id = recipe$prep_config_id[[1]],
          truth_type = truth_type,
          metric = metric,
          bulk_scale = bulk_prep$bulk_scale,
          bulk_normalization = bulk_prep$bulk_normalization,
          bulk_prep_action = bulk_prep$action,
          regress_type = meta$regress_type,
          regress_frac = meta$regress_frac,
          regress_frac_path = meta$regress_frac_path,
          regressed_bulk_file = baseline_bulk_paths[[baseline_id]],
          filter_frac = recipe$filter_frac[[1]],
          filter_frac_path = filter_frac_path,
          sample_policy = "test_samples",
          regression_sample_policy = if (meta$regress_type == "none") {
            NA_character_
          } else {
            "bulk_fraction_overlap"
          },
          min_frac = min_frac,
          min_n_sample = min_n_sample,
          digits = digits,
          n_genes = nrow(result),
          n_cell_types = ncol(result),
          metric_file = out_path,
          stringsAsFactors = FALSE
        )

        message(
          "      wrote ", out_path,
          " (", nrow(result), " genes x ", ncol(result), " cell types)"
        )
      }
    }
  }

  manifest_path <- file.path(out_root, "manifest.tsv")
  write_manifest(manifest_rows, manifest_path)
  if (length(manifest_rows) > 0L) {
    message("  wrote ", manifest_path)
  }
}

message("Done.")

