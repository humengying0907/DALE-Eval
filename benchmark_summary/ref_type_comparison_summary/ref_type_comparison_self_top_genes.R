# Compare independent-reference and self-reference configurations using
# exactly the same self-reference-defined top-gene sets. The script builds
# config01 versus config07 and config02 versus config08.
#
# Notebook-style workflow:
#   1. Regenerate benchmark_dataset_info.txt and acceptable-NA masks after any
#      evaluation-policy change.
#   2. Edit dataset_include/method_include below if needed.
#   3. Source or run this script from repo root or
#      benchmark_summary/ref_type_comparison_summary/.
#
# The output objects are:
#   summary_list_config01_vs_config07[[dataset]]
#   summary_list_config02_vs_config08[[dataset]]
# `summary_list` remains an alias of the config01/config07 result.
#
# Each row represents one paired method/cell_type/top-N comparison. A top-N
# group is selected from Benchmarking_obj/<dataset>/self_reference/
# limma_top_genes.csv before taking any overlap. For CTSE methods, the selected
# genes are restricted to genes present in both config correlation matrices and
# the truth acceptable-NA mask. For bulk_InstaPrismFrac_regressed, they are
# restricted to genes present in both matching config-specific bulk-regressed
# matrices and the same mask. Genes marked as truth-side non-evaluable are
# removed, so n_genes is one shared denominator for the two values in a row.
#
# n_eval_samples is taken from the truth type associated with each pair:
# meancpm/truth_cellfrac for config01/config07 and
# sumcount_cpm/truth_transcriptfrac for config02/config08. It depends on
# dataset and cell type, but not on method, reference type, or top-N group.
#
# The avg_cor_config* columns ignore remaining correlation NAs. Their
# NA-penalized counterparts replace those NAs with zero over the shared
# n_genes denominator. The n_constant_genes_config* columns record those NAs.

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
  stop(
    "Run from repo root or ",
    "benchmark_summary/ref_type_comparison_summary/"
  )
}
repo_root <- normalizePath(repo_root_matches[[1]])

dataset_include <- NULL
method_include <- NULL

top_n_values <- c(10L, 30L, 100L, 300L, 1000L, 3000L, 10000L)

## ---------------------------------------------------------------------------
## Setup
## ---------------------------------------------------------------------------

source(file.path(repo_root, "DALE_Eval", "modules", "config_helpers.R"))
source(file.path(
  repo_root,
  "DALE_Eval",
  "modules",
  "performance_summary_helpers.R"
))

dataset_info <- read_benchmark_dataset_info(repo_root)
required_dataset_info_cols <- c(
  "dataset",
  "cell_type_raw",
  "cell_type",
  "n_eval_samples_meancpm",
  "n_eval_samples_sumcount_cpm"
)
missing_dataset_info_cols <- setdiff(
  required_dataset_info_cols,
  colnames(dataset_info)
)
if (length(missing_dataset_info_cols) > 0L) {
  stop(
    "benchmark_dataset_info.txt missing columns: ",
    paste(missing_dataset_info_cols, collapse = ", "),
    ". Regenerate it with benchmark_summary/preparation/",
    "build_benchmark_dataset_info.R."
  )
}

build_reference_type_comparison <- function(
  config_indep_id,
  config_self_id,
  truth_type,
  n_eval_samples_col
) {

  config_indep <- read_deconv_config(config_indep_id, repo_root = repo_root)
  config_self <- read_deconv_config(config_self_id, repo_root = repo_root)

  if (config_indep$refType[[1]] != "indep") {
    stop(config_indep_id, " must use refType=indep")
  }
  if (config_self$refType[[1]] != "self") {
    stop(config_self_id, " must use refType=self")
  }

  matched_config_fields <- c(
    "bulk_input",
    "bulk_scale",
    "bulk_normalization",
    "frac_input"
  )
  for (field in matched_config_fields) {
    if (!identical(config_indep[[field]][[1]], config_self[[field]][[1]])) {
      stop(
        config_indep_id, " and ", config_self_id,
        " differ in ", field, "; this is not a refType-only comparison"
      )
    }
  }

  performance_dir <- function(dataset, config_row) {
    file.path(
      repo_root,
      "Benchmarking_obj",
      dataset,
      "deconv_performance",
      paste0(config_slug(config_row), "__truth-", truth_type),
      "spearman_cor"
    )
  }

  method_file_map <- function(directory) {
    files <- sort(list.files(directory, pattern = "\\.txt$", full.names = TRUE))
    setNames(files, sub("\\.txt$", "", basename(files)))
  }

  mean_or_na <- function(values) {
    if (length(values) == 0L || all(is.na(values))) {
      return(NA_real_)
    }
    mean(values, na.rm = TRUE)
  }

  mean_with_na_penalty <- function(values) {
    if (length(values) == 0L) {
      return(NA_real_)
    }
    values[is.na(values)] <- 0
    mean(values)
  }

  baseline_path_for <- function(dataset, config_row) {
    specs <- baseline_specs_for_config(config_row, truth_type)
    matched_file <- specs$baseline_file[
      specs$value_col == "avg_cor_bulk_InstaPrismFrac_regressed"
    ]
    if (length(matched_file) != 1L || is.na(matched_file)) {
      stop("No InstaPrism-fraction-regressed baseline specification")
    }

    file.path(
      repo_root,
      "Benchmarking_obj",
      dataset,
      "deconv_performance",
      "bulk_baseline",
      "spearman_cor",
      matched_file
    )
  }

  benchmark_root <- file.path(repo_root, "Benchmarking_obj")
  datasets <- sort(list.dirs(
    benchmark_root,
    full.names = FALSE,
    recursive = FALSE
  ))
  datasets <- datasets[
    dir.exists(file.path(
      benchmark_root,
      datasets,
      "self_reference"
    ))
  ]

  if (!is.null(dataset_include)) {
    datasets <- intersect(datasets, dataset_include)
  }
  if (length(datasets) == 0L) {
    stop("No datasets selected")
  }

  summary_list <- list()

  ## ---------------------------------------------------------------------------
  ## Build Paired Reference-Type Comparisons
  ## ---------------------------------------------------------------------------

  for (dataset_name in datasets) {
    message("Dataset: ", dataset_name)

    cor_dir_indep <- performance_dir(dataset_name, config_indep)
    cor_dir_self <- performance_dir(dataset_name, config_self)
    if (!dir.exists(cor_dir_indep) || !dir.exists(cor_dir_self)) {
      message(
        "  skip: ", config_indep_id, "/", config_self_id,
        " Spearman directories are not both present"
      )
      next
    }

    files_indep <- method_file_map(cor_dir_indep)
    files_self <- method_file_map(cor_dir_self)
    methods <- sort(intersect(names(files_indep), names(files_self)))
    if (!is.null(method_include)) {
      methods <- intersect(methods, method_include)
    }
    if (length(methods) == 0L) {
      message("  skip: no methods are present in both configurations")
      next
    }

    mask <- read_cor_na_acceptable_mask(
      dataset = dataset_name,
      truth_type = truth_type,
      repo_root = repo_root
    )
    self_limma <- read_limma_top_genes_csv(file.path(
      benchmark_root,
      dataset_name,
      "self_reference"
    ))
    bulk_instaprism_regressed_cor_indep <- read_summary_numeric_matrix(
      baseline_path_for(dataset_name, config_indep),
      required = TRUE
    )
    bulk_instaprism_regressed_cor_self <- read_summary_numeric_matrix(
      baseline_path_for(dataset_name, config_self),
      required = TRUE
    )

    methods <- c(methods, "bulk_InstaPrismFrac_regressed")
    dataset_info_one <- dataset_info[
      dataset_info$dataset == dataset_name,
      ,
      drop = FALSE
    ]
    dataset_rows <- list()

    for (method in methods) {
      if (method == "bulk_InstaPrismFrac_regressed") {
        cor_indep <- bulk_instaprism_regressed_cor_indep
        cor_self <- bulk_instaprism_regressed_cor_self
      } else {
        cor_indep <- read_summary_numeric_matrix(
          files_indep[[method]],
          required = TRUE
        )
        cor_self <- read_summary_numeric_matrix(
          files_self[[method]],
          required = TRUE
        )
      }

      common_cell_types <- Reduce(
        intersect,
        list(
          colnames(cor_indep),
          colnames(cor_self),
          colnames(mask),
          colnames(self_limma)
        )
      )

      for (cell_type_raw in common_cell_types) {
        ranked_self_genes <- rank_limma_genes_for_cell_type(
          self_limma,
          cell_type_raw
        )
        if (length(ranked_self_genes) == 0L) {
          next
        }

        shared_available_genes <- Reduce(
          intersect,
          list(
            rownames(cor_indep),
            rownames(cor_self),
            rownames(mask)
          )
        )

        info_index <- match(cell_type_raw, dataset_info_one$cell_type_raw)
        if (is.na(info_index)) {
          stop(
            dataset_name,
            ": cell type missing from benchmark_dataset_info.txt: ",
            cell_type_raw
          )
        }
        cell_type <- dataset_info_one$cell_type[[info_index]]
        n_eval_samples <- dataset_info_one[[n_eval_samples_col]][[info_index]]

        for (top_n in top_n_values) {
          selected_self_genes <- ranked_self_genes[
            seq_len(min(top_n, length(ranked_self_genes)))
          ]
          genes <- selected_self_genes[
            selected_self_genes %in% shared_available_genes
          ]
          genes <- genes[mask[genes, cell_type_raw] != 1]

          values_indep <- cor_indep[genes, cell_type_raw]
          values_self <- cor_self[genes, cell_type_raw]

          summary_row <- data.frame(
            method = method,
            cell_type = cell_type,
            group = paste0("top_", top_n),
            n_genes = length(genes),
            n_eval_samples = n_eval_samples,
            stringsAsFactors = FALSE
          )
          summary_row[[paste0(
            "n_constant_genes_", config_indep_id
          )]] <- sum(is.na(values_indep))
          summary_row[[paste0(
            "n_constant_genes_", config_self_id
          )]] <- sum(is.na(values_self))
          summary_row[[paste0("avg_cor_", config_indep_id)]] <- mean_or_na(
            values_indep
          )
          summary_row[[paste0("avg_cor_", config_self_id)]] <- mean_or_na(
            values_self
          )
          summary_row[[paste0(
            "avg_cor_with_NA_penalty_", config_indep_id
          )]] <- mean_with_na_penalty(values_indep)
          summary_row[[paste0(
            "avg_cor_with_NA_penalty_", config_self_id
          )]] <- mean_with_na_penalty(values_self)

          dataset_rows[[length(dataset_rows) + 1L]] <- summary_row
        }
      }
    }

    if (length(dataset_rows) == 0L) {
      summary_list[[dataset_name]] <- data.frame()
    } else {
      summary_list[[dataset_name]] <- do.call(rbind, dataset_rows)
      rownames(summary_list[[dataset_name]]) <- NULL
    }
  }

  message(
    "Built ", config_indep_id, "/", config_self_id,
    " reference-type comparison for ", length(summary_list), " dataset(s)."
  )
  summary_list
}

## ---------------------------------------------------------------------------
## Build And Save Both Reference-Type Comparisons
## ---------------------------------------------------------------------------

summary_list_config01_vs_config07 <- build_reference_type_comparison(
  config_indep_id = "config01",
  config_self_id = "config07",
  truth_type = "meancpm",
  n_eval_samples_col = "n_eval_samples_meancpm"
)

summary_output_path_config01_vs_config07 <- file.path(
  repo_root,
  "benchmark_summary",
  "ref_type_comparison_summary",
  "config01_vs_config07_self_top_genes_summary_list.RDS"
)
saveRDS(
  summary_list_config01_vs_config07,
  file = summary_output_path_config01_vs_config07
)
message("Saved summary list: ", summary_output_path_config01_vs_config07)

summary_list_config02_vs_config08 <- build_reference_type_comparison(
  config_indep_id = "config02",
  config_self_id = "config08",
  truth_type = "sumcount_cpm",
  n_eval_samples_col = "n_eval_samples_sumcount_cpm"
)

summary_output_path_config02_vs_config08 <- file.path(
  repo_root,
  "benchmark_summary",
  "ref_type_comparison_summary",
  "config02_vs_config08_self_top_genes_summary_list.RDS"
)
saveRDS(
  summary_list_config02_vs_config08,
  file = summary_output_path_config02_vs_config08
)
message("Saved summary list: ", summary_output_path_config02_vs_config08)

summary_list <- summary_list_config01_vs_config07
message("Examples: summary_list_config01_vs_config07[[1]] and ",
        "summary_list_config02_vs_config08[[1]]")
