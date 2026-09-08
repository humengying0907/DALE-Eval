# Summarize the effect of high-fraction sample filtering on CTSE gene correlation.
#
# Notebook-style workflow:
#   1. Edit the settings block if needed.
#   2. Run from repo root or benchmark_summary/sample_filter_effect_summary/.
#   3. Inspect the messages and the RDS written beside this script.
#
# Comparison:
#   - all-sample correlations from spearman_cor/
#   - filtered correlations from spearman_cor_highFrac_sample/
#   - config01 uses meancpm truth and truth_cellfrac > 0.001
#   - config02 uses sumcount_cpm truth and truth_transcriptfrac > 0.001
#   - filtered samples additionally require the config-specific InstaPrism
#     fraction mask to equal 1, corresponding to InstaPrism fraction > 0.1
#
# Gene selection:
#   - top-N genes are ranked by the independent-reference limma statistic
#   - the original top N are selected before truth eligibility is applied
#   - filtered-sample truth eligibility requires at least min_n_sample samples
#     and nonconstant truth CTSE across those filtered samples
#   - one exact gene set is shared by the paired method correlations and all
#     four bulk-baseline correlations in each output row
#
# Sample counts:
#   - n_eval_samples_all is recalculated as:
#       test samples intersect truth CTSE samples intersect truth fraction > 0.001
#   - n_eval_samples_filtered additionally intersects the high-fraction mask
#   - these counts match evaluation under the benchmark contract that every
#     method exports all test samples
#
# NA summaries:
#   - unpenalized averages omit non-finite method correlations
#   - penalty averages replace non-finite method correlations with na_penalty
#   - when filtered truth has fewer than min_n_sample samples, rows are retained
#     with n_genes = 0 and correlation summaries equal to NA
#
# Output:
#   benchmark_summary/sample_filter_effect_summary/
#     all_vs_highFrac_samples_summary_list.RDS
#
# Each dataset element contains:
#   method, cell_type, group, n_genes,
#   n_eval_samples_all, n_eval_samples_filtered,
#   avg_cor_all_samples, avg_cor_all_samples_with_NA_penalty,
#   avg_cor_filtered_samples, avg_cor_filtered_samples_with_NA_penalty,
#   avg_cor_bulk_all_samples, avg_cor_bulk_filtered_samples,
#   avg_cor_bulk_truthFrac_regressed_all_samples,
#   avg_cor_bulk_truthFrac_regressed_filtered_samples,
#   config

## ---------------------------------------------------------------------------
## Settings
## ---------------------------------------------------------------------------

dataset_include <- NULL
config_include <- c("config01", "config02")
method_include <- NULL

top_n_values <- c(10L, 30L, 100L, 300L, 1000L, 3000L, 10000L)
min_frac <- 0.001
min_n_sample <- 10L
na_penalty <- 0
digits <- 6L

## ---------------------------------------------------------------------------
## Setup
## ---------------------------------------------------------------------------

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- if (length(script_arg) > 0L) {
  sub("^--file=", "", script_arg[[1]])
} else {
  ""
}
script_dir <- if (nzchar(script_path)) {
  dirname(normalizePath(script_path, mustWork = FALSE))
} else {
  normalizePath(".", mustWork = FALSE)
}

repo_root_candidates <- unique(normalizePath(
  c(".", "..", "../..", file.path(script_dir, "..", "..")),
  mustWork = FALSE
))
repo_root_matches <- repo_root_candidates[
  file.exists(file.path(repo_root_candidates, "Benchmarking_obj")) &
    file.exists(file.path(repo_root_candidates, "DALE_Eval", "configs", "deconv_configs.txt"))
]
if (length(repo_root_matches) == 0L) {
  stop("Run from repo root or benchmark_summary/sample_filter_effect_summary/")
}
repo_root <- normalizePath(repo_root_matches[[1]])

validate_values <- function(values, allowed, label) {
  bad <- setdiff(values, allowed)
  if (length(bad) > 0L) {
    stop(label, " contains unsupported values: ", paste(bad, collapse = ", "))
  }
}

validate_values(config_include, c("config01", "config02"), "config_include")

benchmark_root <- file.path(repo_root, "Benchmarking_obj")
assignment_path <- file.path(
  repo_root,
  "DALE_Eval",
  "configs",
  "benchmark_ref_assignment.txt"
)
mapping_path <- file.path(repo_root, "Indep_scReference", "cell_type_mapping.txt")
deconv_config_path <- file.path(
  repo_root,
  "DALE_Eval",
  "configs",
  "deconv_configs.txt"
)

assignments <- read.delim(
  assignment_path,
  sep = "\t",
  stringsAsFactors = FALSE,
  check.names = FALSE
)
cell_type_mapping <- read.delim(
  mapping_path,
  sep = "\t",
  stringsAsFactors = FALSE,
  check.names = FALSE
)
deconv_configs <- read.delim(
  deconv_config_path,
  sep = "\t",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

config_specs <- data.frame(
  config = c("config01", "config02"),
  truth_type = c("meancpm", "sumcount_cpm"),
  filter_frac = c("truth_cellfrac", "truth_transcriptfrac"),
  sample_mask_file = c(
    "InstaPrismfrac_config01_gt0.1.txt",
    "InstaPrismfrac_config02_gt0.1.txt"
  ),
  bulk_recipe = c("bulk-wcpm", "bulk-sumcount__norm-cpm"),
  truth_regress_frac = c("truth_cellfrac", "truth_transcriptfrac"),
  stringsAsFactors = FALSE
)
config_specs <- config_specs[
  config_specs$config %in% config_include,
  ,
  drop = FALSE
]

group_names <- c(
  paste0("top_", top_n_values),
  "all_genes"
)

datasets <- list.dirs(benchmark_root, full.names = FALSE, recursive = FALSE)
datasets <- datasets[file.exists(file.path(benchmark_root, datasets, "deconv_performance"))]
if (!is.null(dataset_include)) {
  datasets <- intersect(datasets, dataset_include)
}
datasets <- sort(datasets)

if (length(datasets) == 0L) {
  stop("No datasets selected")
}

## ---------------------------------------------------------------------------
## Helpers
## ---------------------------------------------------------------------------

config_slug <- function(config_id) {
  config_row <- deconv_configs[
    deconv_configs$config_id == config_id,
    ,
    drop = FALSE
  ]
  if (nrow(config_row) != 1L) {
    stop("Expected exactly one deconvolution config row for ", config_id)
  }

  slug <- paste0(
    config_id,
    "_bulk-", config_row$bulk_input[[1]],
    "__frac-", config_row$frac_input[[1]],
    "__ref-", config_row$refType[[1]]
  )
  if (tolower(config_row$bulk_normalization[[1]]) != "cpm") {
    slug <- paste0(
      slug,
      "__norm-",
      config_row$bulk_normalization[[1]]
    )
  }
  slug
}

read_numeric_matrix <- function(path, label, na_to_zero = FALSE) {
  if (!file.exists(path)) {
    stop(label, " not found: ", path)
  }
  x <- read.delim(
    path,
    sep = "\t",
    check.names = FALSE,
    row.names = 1
  )
  x <- as.matrix(x)
  storage.mode(x) <- "numeric"
  if (na_to_zero) {
    x[is.na(x)] <- 0
  }
  x
}

read_truth_files <- function(truth_dir) {
  paths <- list.files(
    truth_dir,
    pattern = "[.]txt([.]gz)?$",
    full.names = TRUE
  )
  stats::setNames(
    paths,
    sub("[.]txt([.]gz)?$", "", basename(paths))
  )
}

read_named_cor_files <- function(cor_dir) {
  if (!dir.exists(cor_dir)) {
    return(character())
  }
  paths <- list.files(
    cor_dir,
    pattern = "[.]txt$",
    full.names = TRUE
  )
  stats::setNames(
    paths,
    sub("[.]txt$", "", basename(paths))
  )
}

mean_correlation <- function(values) {
  values <- as.numeric(values)
  finite_values <- values[is.finite(values)]
  if (length(finite_values) == 0L) {
    return(NA_real_)
  }
  round(mean(finite_values), digits)
}

mean_correlation_with_penalty <- function(values) {
  values <- as.numeric(values)
  if (length(values) == 0L) {
    return(NA_real_)
  }
  values[!is.finite(values)] <- na_penalty
  round(mean(values), digits)
}

truth_variable_genes <- function(truth, samples) {
  if (length(samples) < min_n_sample) {
    return(character())
  }

  truth_eval <- truth[, samples, drop = FALSE]
  variable <- apply(
    truth_eval,
    1L,
    function(values) {
      values <- values[is.finite(values)]
      length(values) >= min_n_sample && length(unique(values)) > 1L
    }
  )
  rownames(truth_eval)[variable]
}

empty_summary_table <- function() {
  data.frame(
    method = character(),
    cell_type = character(),
    group = character(),
    n_genes = integer(),
    n_eval_samples_all = integer(),
    n_eval_samples_filtered = integer(),
    avg_cor_all_samples = numeric(),
    avg_cor_all_samples_with_NA_penalty = numeric(),
    avg_cor_filtered_samples = numeric(),
    avg_cor_filtered_samples_with_NA_penalty = numeric(),
    avg_cor_bulk_all_samples = numeric(),
    avg_cor_bulk_filtered_samples = numeric(),
    avg_cor_bulk_truthFrac_regressed_all_samples = numeric(),
    avg_cor_bulk_truthFrac_regressed_filtered_samples = numeric(),
    config = character(),
    stringsAsFactors = FALSE
  )
}

## ---------------------------------------------------------------------------
## Build Paired All-Sample Versus High-Fraction Summaries
## ---------------------------------------------------------------------------

summary_list <- stats::setNames(vector("list", length(datasets)), datasets)

for (dataset_name in datasets) {
  message("Dataset: ", dataset_name)

  dataset_root <- file.path(benchmark_root, dataset_name)
  split_path <- file.path(dataset_root, "self_reference", "sample_split.txt")
  assignment_row <- assignments[
    assignments$dataset == dataset_name,
    ,
    drop = FALSE
  ]

  if (!file.exists(split_path)) {
    message("  skip dataset: missing sample split ", split_path)
    summary_list[[dataset_name]] <- empty_summary_table()
    next
  }
  if (nrow(assignment_row) != 1L) {
    message("  skip dataset: expected one independent-reference assignment")
    summary_list[[dataset_name]] <- empty_summary_table()
    next
  }

  split_table <- read.delim(
    split_path,
    sep = "\t",
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  if (!all(c("group", "sampleIDs") %in% colnames(split_table))) {
    stop("Unexpected sample split columns: ", split_path)
  }
  test_samples <- unique(as.character(
    split_table$sampleIDs[split_table$group == "test"]
  ))

  indep_ref <- assignment_row$indep_ref[[1]]
  limma_path <- file.path(
    repo_root,
    "Indep_scReference",
    indep_ref,
    "limma_top_genes.csv"
  )
  if (!file.exists(limma_path)) {
    message("  skip dataset: missing independent-reference limma file ", limma_path)
    summary_list[[dataset_name]] <- empty_summary_table()
    next
  }

  limma_stats <- read.csv(
    limma_path,
    row.names = 1,
    check.names = FALSE
  )
  limma_stats <- as.matrix(limma_stats)
  storage.mode(limma_stats) <- "numeric"

  dataset_mapping <- cell_type_mapping[
    cell_type_mapping$dataset == dataset_name &
      cell_type_mapping$indep_ref == indep_ref,
    ,
    drop = FALSE
  ]

  dataset_rows <- list()

  for (config_idx in seq_len(nrow(config_specs))) {
    spec <- config_specs[config_idx, , drop = FALSE]
    config_id <- spec$config[[1]]
    truth_type <- spec$truth_type[[1]]
    recipe_id <- spec$bulk_recipe[[1]]
    regress_frac_name <- spec$truth_regress_frac[[1]]

    message("  Config: ", config_id)

    performance_dir <- file.path(
      dataset_root,
      "deconv_performance",
      paste0(config_slug(config_id), "__truth-", truth_type)
    )
    all_cor_dir <- file.path(performance_dir, "spearman_cor")
    filtered_cor_dir <- file.path(
      performance_dir,
      "spearman_cor_highFrac_sample"
    )
    truth_dir <- file.path(dataset_root, "ctse_truth", truth_type)
    filter_frac_path <- file.path(
      dataset_root,
      "frac_input",
      paste0(spec$filter_frac[[1]], ".txt")
    )
    sample_mask_path <- file.path(
      dataset_root,
      "evaluation_metadata",
      "sample_mask",
      spec$sample_mask_file[[1]]
    )

    baseline_all_dir <- file.path(
      dataset_root,
      "deconv_performance",
      "bulk_baseline",
      "spearman_cor"
    )
    baseline_filtered_dir <- file.path(
      dataset_root,
      "deconv_performance",
      "bulk_baseline",
      "spearman_cor_highFrac_sample"
    )
    baseline_raw_name <- paste0(
      "truth-", truth_type,
      "__", recipe_id,
      ".txt"
    )
    baseline_regressed_name <- paste0(
      "truth-", truth_type,
      "__", recipe_id,
      "__regress-", regress_frac_name,
      ".txt"
    )
    baseline_paths <- c(
      bulk_all = file.path(baseline_all_dir, baseline_raw_name),
      bulk_filtered = file.path(baseline_filtered_dir, baseline_raw_name),
      regressed_all = file.path(baseline_all_dir, baseline_regressed_name),
      regressed_filtered = file.path(
        baseline_filtered_dir,
        baseline_regressed_name
      )
    )

    required_paths <- c(
      "all-sample correlation directory" = all_cor_dir,
      "filtered correlation directory" = filtered_cor_dir,
      "truth directory" = truth_dir,
      "fraction filter" = filter_frac_path,
      "sample mask" = sample_mask_path,
      baseline_paths
    )
    missing_required <- names(required_paths)[!file.exists(required_paths)]
    if (length(missing_required) > 0L) {
      message(
        "    skip config: missing ",
        paste(missing_required, collapse = ", ")
      )
      next
    }

    all_method_files <- read_named_cor_files(all_cor_dir)
    filtered_method_files <- read_named_cor_files(filtered_cor_dir)
    methods <- sort(intersect(
      names(all_method_files),
      names(filtered_method_files)
    ))
    if (!is.null(method_include)) {
      methods <- intersect(methods, method_include)
    }
    if (length(methods) == 0L) {
      message("    skip config: no paired method correlation files")
      next
    }

    truth_files <- read_truth_files(truth_dir)
    filter_frac <- read_numeric_matrix(
      filter_frac_path,
      "truth fraction filter"
    )
    sample_mask <- read_numeric_matrix(
      sample_mask_path,
      "high-fraction sample mask"
    )
    baseline_bulk_all <- read_numeric_matrix(
      baseline_paths[["bulk_all"]],
      "all-sample bulk baseline"
    )
    baseline_bulk_filtered <- read_numeric_matrix(
      baseline_paths[["bulk_filtered"]],
      "filtered bulk baseline"
    )
    baseline_regressed_all <- read_numeric_matrix(
      baseline_paths[["regressed_all"]],
      "all-sample truth-fraction-regressed baseline"
    )
    baseline_regressed_filtered <- read_numeric_matrix(
      baseline_paths[["regressed_filtered"]],
      "filtered truth-fraction-regressed baseline"
    )

    mapped_cell_types <- intersect(
      dataset_mapping$target_cell_type,
      names(truth_files)
    )
    mapped_cell_types <- intersect(mapped_cell_types, colnames(filter_frac))
    mapped_cell_types <- intersect(mapped_cell_types, colnames(sample_mask))
    mapped_cell_types <- Reduce(
      intersect,
      list(
        mapped_cell_types,
        colnames(baseline_bulk_all),
        colnames(baseline_bulk_filtered),
        colnames(baseline_regressed_all),
        colnames(baseline_regressed_filtered)
      )
    )
    mapped_cell_types <- sort(unique(mapped_cell_types))

    for (cell_type in mapped_cell_types) {
      mapping_row <- dataset_mapping[
        dataset_mapping$target_cell_type == cell_type,
        ,
        drop = FALSE
      ]
      if (nrow(mapping_row) != 1L) {
        message(
          "    skip cell type ", cell_type,
          ": expected one independent-reference mapping"
        )
        next
      }

      indep_cell_type <- mapping_row$indep_ref_cell_type[[1]]
      if (!indep_cell_type %in% colnames(limma_stats)) {
        message(
          "    skip cell type ", cell_type,
          ": limma column not found for ", indep_cell_type
        )
        next
      }

      truth <- read_numeric_matrix(
        truth_files[[cell_type]],
        paste0("truth CTSE for ", cell_type),
        na_to_zero = TRUE
      )
      fraction_values <- filter_frac[, cell_type]
      fraction_pass <- rownames(filter_frac)[
        is.finite(fraction_values) & fraction_values > min_frac
      ]
      all_eval_samples <- Reduce(
        intersect,
        list(test_samples, colnames(truth), fraction_pass)
      )

      mask_values <- sample_mask[, cell_type]
      mask_pass <- rownames(sample_mask)[
        is.finite(mask_values) & mask_values == 1
      ]
      filtered_eval_samples <- intersect(all_eval_samples, mask_pass)

      eligible_genes <- truth_variable_genes(
        truth,
        filtered_eval_samples
      )

      limma_values <- limma_stats[, indep_cell_type]
      ranked_genes <- rownames(limma_stats)[order(
        limma_values,
        decreasing = TRUE,
        na.last = NA
      )]
      ranked_genes <- unique(ranked_genes)

      selected_gene_groups <- lapply(
        top_n_values,
        function(top_n) head(ranked_genes, top_n)
      )
      names(selected_gene_groups) <- paste0("top_", top_n_values)
      selected_gene_groups[["all_genes"]] <- ranked_genes

      for (method_name in methods) {
        method_all <- read_numeric_matrix(
          all_method_files[[method_name]],
          paste0("all-sample correlation for ", method_name)
        )
        method_filtered <- read_numeric_matrix(
          filtered_method_files[[method_name]],
          paste0("filtered correlation for ", method_name)
        )

        required_cell_type_matrices <- list(
          method_all,
          method_filtered
        )
        if (!all(vapply(
          required_cell_type_matrices,
          function(x) cell_type %in% colnames(x),
          logical(1)
        ))) {
          next
        }

        available_genes <- Reduce(
          intersect,
          list(
            eligible_genes,
            rownames(method_all),
            rownames(method_filtered),
            rownames(baseline_bulk_all),
            rownames(baseline_bulk_filtered),
            rownames(baseline_regressed_all),
            rownames(baseline_regressed_filtered)
          )
        )

        for (group_name in group_names) {
          genes <- intersect(
            selected_gene_groups[[group_name]],
            available_genes
          )
          n_genes <- length(genes)

          if (n_genes == 0L) {
            method_all_values <- numeric()
            method_filtered_values <- numeric()
            bulk_all_values <- numeric()
            bulk_filtered_values <- numeric()
            regressed_all_values <- numeric()
            regressed_filtered_values <- numeric()
          } else {
            method_all_values <- method_all[genes, cell_type]
            method_filtered_values <- method_filtered[genes, cell_type]
            bulk_all_values <- baseline_bulk_all[genes, cell_type]
            bulk_filtered_values <- baseline_bulk_filtered[genes, cell_type]
            regressed_all_values <- baseline_regressed_all[genes, cell_type]
            regressed_filtered_values <- baseline_regressed_filtered[
              genes,
              cell_type
            ]
          }

          dataset_rows[[length(dataset_rows) + 1L]] <- data.frame(
            method = method_name,
            cell_type = cell_type,
            group = group_name,
            n_genes = n_genes,
            n_eval_samples_all = length(all_eval_samples),
            n_eval_samples_filtered = length(filtered_eval_samples),
            avg_cor_all_samples = mean_correlation(method_all_values),
            avg_cor_all_samples_with_NA_penalty =
              mean_correlation_with_penalty(method_all_values),
            avg_cor_filtered_samples =
              mean_correlation(method_filtered_values),
            avg_cor_filtered_samples_with_NA_penalty =
              mean_correlation_with_penalty(method_filtered_values),
            avg_cor_bulk_all_samples =
              mean_correlation(bulk_all_values),
            avg_cor_bulk_filtered_samples =
              mean_correlation(bulk_filtered_values),
            avg_cor_bulk_truthFrac_regressed_all_samples =
              mean_correlation(regressed_all_values),
            avg_cor_bulk_truthFrac_regressed_filtered_samples =
              mean_correlation(regressed_filtered_values),
            config = config_id,
            stringsAsFactors = FALSE
          )
        }
      }
    }
  }

  summary_list[[dataset_name]] <- if (length(dataset_rows) > 0L) {
    do.call(rbind, dataset_rows)
  } else {
    empty_summary_table()
  }

  message(
    "  Summary rows: ",
    nrow(summary_list[[dataset_name]])
  )
}

output_path <- file.path(
  script_dir,
  "all_vs_highFrac_samples_summary_list.RDS"
)
saveRDS(summary_list, output_path)
message("Saved: ", output_path)
