# Summarize method Spearman correlations by limma top-gene groups.
#
# Notebook-style workflow:
#   1. Edit dataset_include/config_include in Settings.
#   2. Source or run this script from repo root or benchmark_summary/ctse_spearman_cor_summary/.
#   3. Inspect summary_list; the final saveRDS line can be commented out for inspection-only runs.
#
# The output object is:
#   summary_list[[dataset]]
#
# Each element is a long table with one row per method/cell_type/group/config/truth.
# Matching bulk baseline averages are repeated as annotation columns.

## ---------------------------------------------------------------------------
## Settings
## ---------------------------------------------------------------------------

# config01, config02 → independent-reference limma rankings
# config07, config08, config09 → self-reference limma rankings
# n_genes: top n genes overlap with method output, after removing genes marked as non-evaluable by the acceptable-NA mask
# n_constant_genes: among n_genes, genes with Spearman_cor = NA
# n_eval_samples: testing samples present in the matching CTSE truth with truth fraction > 0.001;
#                 it depends on dataset/cell type/truth type, not method or top-n group
#                 (meancpm uses truth_cellfrac; sumcount_cpm uses truth_transcriptfrac)
# An acceptable NA is a gene correlation that is undefined for a truth-side reason and therefore should not penalize the method
# For a gene × cell type, the mask marks acceptable NA = 1 when:
# - Fewer than 10 test samples have the corresponding fraction > 0.001, or
# - Truth CTSE is constant across those eligible test samples.

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
    "benchmark_summary/ctse_spearman_cor_summary/"
  )
}
repo_root <- normalizePath(repo_root_matches[[1]])

# Set either value to NULL to summarize all available datasets/configs.
# Keep the defaults small for quick notebook testing.
#dataset_include <- c("BRCA_Bassez2021")
#config_include <- c("config01",'config02')

dataset_include = NULL
config_include = NULL

# Optional filters. Set to NULL to include all available methods/truth types.
method_include <- NULL
truth_type_include <- NULL

top_n_values <- c(10L, 30L, 100L, 300L, 1000L, 3000L, 10000L)

## ---------------------------------------------------------------------------
## Setup
## ---------------------------------------------------------------------------

source(file.path(repo_root, "DALE_Eval", "modules", "config_helpers.R"))
source(file.path(repo_root, "DALE_Eval", "modules", "performance_summary_helpers.R"))

dataset_info <- read_benchmark_dataset_info(repo_root)

eval_sample_col_by_truth <- c(
  meancpm = "n_eval_samples_meancpm",
  sumcount_cpm = "n_eval_samples_sumcount_cpm"
)
required_eval_sample_cols <- unname(eval_sample_col_by_truth)
missing_eval_sample_cols <- setdiff(
  required_eval_sample_cols,
  colnames(dataset_info)
)
if (length(missing_eval_sample_cols) > 0L) {
  stop(
    "benchmark_dataset_info.txt missing evaluation-sample columns: ",
    paste(missing_eval_sample_cols, collapse = ", "),
    ". Regenerate it with benchmark_summary/preparation/",
    "build_benchmark_dataset_info.R."
  )
}

deconv_config_path <- file.path(repo_root, "DALE_Eval", "configs", "deconv_configs.txt")
deconv_configs <- read.delim(deconv_config_path, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)

required_config_cols <- c("config_id", "bulk_input", "bulk_scale", "bulk_normalization", "frac_input", "refType")
missing_config_cols <- setdiff(required_config_cols, colnames(deconv_configs))
if (length(missing_config_cols) > 0L) {
  stop("deconv_configs.txt missing columns: ", paste(missing_config_cols, collapse = ", "))
}

benchmark_root <- file.path(repo_root, "Benchmarking_obj")
datasets <- sort(list.dirs(benchmark_root, full.names = FALSE, recursive = FALSE))
datasets <- datasets[file.exists(file.path(benchmark_root, datasets, "deconv_performance"))]

if (!is.null(dataset_include)) {
  datasets <- intersect(datasets, dataset_include)
}

if (length(datasets) == 0L) {
  stop("No datasets selected")
}

summary_list <- list()

## ---------------------------------------------------------------------------
## Build Dataset-Level Summary Tables
## ---------------------------------------------------------------------------

for (dataset_name in datasets) {
  message("Dataset: ", dataset_name)

  config_rows <- deconv_configs
  if (!is.null(config_include)) {
    config_rows <- config_rows[config_rows$config_id %in% config_include, , drop = FALSE]
  }

  dataset_tables <- list()

  for (i in seq_len(nrow(config_rows))) {
    config_row <- config_rows[i, , drop = FALSE]
    config_id <- config_row$config_id[[1]]
    slug <- config_slug(config_row)

    performance_dirs <- find_performance_dirs_for_config(
      dataset = dataset_name,
      config_row = config_row,
      truth_type_include = truth_type_include,
      repo_root = repo_root
    )

    if (nrow(performance_dirs) == 0L) {
      next
    }

    for (j in seq_len(nrow(performance_dirs))) {
      performance_dir <- performance_dirs$performance_dir[[j]]
      truth_type <- performance_dirs$truth_type[[j]]
      cor_dir <- file.path(performance_dir, "spearman_cor")

      if (!dir.exists(cor_dir)) {
        next
      }

      cor_files <- sort(list.files(cor_dir, pattern = "\\.txt$", full.names = TRUE))
      if (length(cor_files) == 0L) {
        next
      }

      methods <- sub("\\.txt$", "", basename(cor_files))
      if (!is.null(method_include)) {
        keep <- methods %in% method_include
        cor_files <- cor_files[keep]
        methods <- methods[keep]
      }

      if (length(cor_files) == 0L) {
        next
      }

      mask <- read_cor_na_acceptable_mask(
        dataset = dataset_name,
        truth_type = truth_type,
        repo_root = repo_root
      )

      all_cell_types <- sort(unique(unlist(lapply(cor_files, function(path) {
        colnames(read_summary_numeric_matrix(path, required = TRUE))
      }), use.names = FALSE)))

      top_groups_by_cell_type <- build_top_gene_groups_by_cell_type(
        dataset = dataset_name,
        config_row = config_row,
        cell_types_raw = all_cell_types,
        top_n_values = top_n_values,
        repo_root = repo_root
      )

      baseline_summary <- summarize_baselines_for_config(
        dataset = dataset_name,
        config_row = config_row,
        truth_type = truth_type,
        mask = mask,
        dataset_info = dataset_info,
        top_groups_by_cell_type = top_groups_by_cell_type,
        top_n_values = top_n_values,
        repo_root = repo_root
      )

      message("  ", config_id, " / ", truth_type, " / methods=", length(methods))

      for (k in seq_along(cor_files)) {
        table_one <- summarize_spearman_file(
          dataset = dataset_name,
          method = methods[[k]],
          cor_path = cor_files[[k]],
          mask = mask,
          config_id = config_id,
          config_row = config_row,
          truth_type = truth_type,
          dataset_info = dataset_info,
          top_groups_by_cell_type = top_groups_by_cell_type,
          top_n_values = top_n_values
        )

        table_one <- merge(
          table_one,
          baseline_summary,
          by = c("cell_type", "group"),
          all.x = TRUE,
          sort = FALSE
        )

        eval_sample_col <- unname(eval_sample_col_by_truth[truth_type])
        if (is.na(eval_sample_col)) {
          stop(
            "No evaluation-sample column configured for truth_type: ",
            truth_type
          )
        }
        dataset_eval_info <- dataset_info[
          dataset_info$dataset == dataset_name,
          ,
          drop = FALSE
        ]
        table_one$n_eval_samples <- dataset_eval_info[[eval_sample_col]][
          match(table_one$cell_type, dataset_eval_info$cell_type)
        ]

        table_one <- table_one[, c(
          "method",
          "cell_type",
          "group",
          "avg_cor",
          "avg_cor_with_NA_penalty",
          "n_genes",
          "n_constant_genes",
          "n_eval_samples",
          "avg_cor_bulk",
          "avg_cor_bulk_truthFrac_regressed",
          "avg_cor_bulk_InstaPrismFrac_regressed",
          "refType",
          "bulk_input",
          "bulk_normalization",
          "frac_input",
          "config_id",
          "truth_type"
        )]

        dataset_tables[[length(dataset_tables) + 1L]] <- table_one
      }
    }
  }

  if (length(dataset_tables) == 0L) {
    summary_list[[dataset_name]] <- data.frame()
  } else {
    summary_list[[dataset_name]] <- do.call(rbind, dataset_tables)
    rownames(summary_list[[dataset_name]]) <- NULL
  }
}

message("Built summary_list for ", length(summary_list), " dataset(s).")
message("Example: summary_list[[1]]")

summary_output_path <- file.path(
  repo_root,
  "benchmark_summary",
  "ctse_spearman_cor_summary",
  "spearman_by_top_genes_summary_list.RDS"
)
saveRDS(summary_list, file = summary_output_path)
message("Saved summary list: ", summary_output_path)
