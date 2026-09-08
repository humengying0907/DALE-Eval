# Summarize post-hoc versus independent-reference gene prioritization.
#
# Notebook-style workflow:
#   1. Edit the settings block below.
#   2. Run from the repository root or this summary directory.
#   3. Inspect summary_list, summary_df, and availability_summary.
#   4. The final section writes posthoc_vs_indep_ref_summary_list.RDS.
#
# Comparison contract:
#   - all selected method-specific post-hoc rankings and adapted GeneSigTest
#     are compared with the mapped independent-reference limma ranking;
#   - within each comparison, the alternative and independent-reference
#     rankings use the same truth-independent candidate universe;
#   - top N selection occurs before truth-derived acceptable-NA filtering;
#   - selected genes are evaluated with the existing method Spearman matrix;
#   - larger post-hoc and limma values rank first, while smaller GeneSigTest
#     FDR values rank first;
#   - ordinary and zero-penalized average correlations are both retained.

## ---------------------------------------------------------------------------
## Settings
## ---------------------------------------------------------------------------

# NULL includes every dataset with a gene_prioritization directory.
dataset_include = NULL

# Defaults to config01, which currently has the broadest method/dataset
# coverage. Add "config02" here to summarize both input-generation settings.
config_include = c("config01")

# NULL includes every method with a Spearman matrix. Missing post-hoc sources
# are skipped individually.
method_include = NULL

# Selected method-specific post-hoc scores. The default retains all currently
# generated post-hoc rankings. GeneSigTest is config-level and
# method-independent, and is included separately when requested below.
posthoc_scores = c(
  "meanlog_margin",
  "meanlog_mean_contrast",
  "logmean_margin",
  "logmean_contrast",
  "partial_r2"
)
include_genesigtest = TRUE

top_n_values = c(
  10L,
  30L,
  100L,
  300L,
  1000L,
  3000L
)

output_file = "posthoc_vs_indep_ref_summary_list.RDS"


## ---------------------------------------------------------------------------
## Repository setup
## ---------------------------------------------------------------------------

repo_root_candidates = unique(normalizePath(
  c(getwd(), file.path(getwd(), ".."), file.path(getwd(), "../..")),
  mustWork = FALSE
))

repo_root = repo_root_candidates[file.exists(file.path(
  repo_root_candidates,
  "DALE_Eval", "configs", "deconv_configs.txt"
))][1]

if (is.na(repo_root)) {
  stop("Could not locate the ctse_benchmark repository root")
}

source(file.path(
  repo_root,
  "DALE_Eval",
  "modules",
  "config_helpers.R"
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
source(file.path(
  repo_root,
  "DALE_Eval",
  "modules",
  "performance_summary_helpers.R"
))

supported_posthoc_scores = c(
  "meanlog_margin",
  "meanlog_mean_contrast",
  "logmean_margin",
  "logmean_contrast",
  "partial_r2"
)

if (
  length(posthoc_scores) == 0L ||
    anyNA(posthoc_scores) ||
    any(!posthoc_scores %in% supported_posthoc_scores)
) {
  stop(
    "posthoc_scores must contain one or more of: ",
    paste(supported_posthoc_scores, collapse = ", ")
  )
}

posthoc_scores = unique(posthoc_scores)

if (
  length(include_genesigtest) != 1L ||
    !is.logical(include_genesigtest) ||
    is.na(include_genesigtest)
) {
  stop("include_genesigtest must be TRUE or FALSE")
}

if (
  length(top_n_values) == 0L ||
    anyNA(top_n_values) ||
    any(top_n_values < 1L)
) {
  stop("top_n_values must contain positive integers")
}

top_n_values = unique(as.integer(top_n_values))

posthoc_score_labels = c(
  "meanlog_margin" = "Post-hoc mean-log margin",
  "meanlog_mean_contrast" = "Post-hoc mean-log mean contrast",
  "logmean_margin" = "Post-hoc log-mean margin",
  "logmean_contrast" = "Post-hoc log-mean contrast",
  "partial_r2" = "Post-hoc partial R2"
)

ranking_score_labels = c(
  "indep_ref_limma" = "Independent reference",
  posthoc_score_labels,
  "genesigtest_adapted_fdr" = "Adapted GeneSigTest FDR"
)

ranking_source_by_score = c(
  "indep_ref_limma" = "indep_ref",
  "meanlog_margin" = "posthoc",
  "meanlog_mean_contrast" = "posthoc",
  "logmean_margin" = "posthoc",
  "logmean_contrast" = "posthoc",
  "partial_r2" = "posthoc",
  "genesigtest_adapted_fdr" = "genesigtest"
)

ranking_decreasing_by_score = c(
  "indep_ref_limma" = TRUE,
  "meanlog_margin" = TRUE,
  "meanlog_mean_contrast" = TRUE,
  "logmean_margin" = TRUE,
  "logmean_contrast" = TRUE,
  "partial_r2" = TRUE,
  "genesigtest_adapted_fdr" = FALSE
)

deconv_configs = read.delim(
  file.path(
    repo_root,
    "DALE_Eval",
    "configs",
    "deconv_configs.txt"
  ),
  sep = "\t",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

eval_configs = read.delim(
  file.path(
    repo_root,
    "DALE_Eval",
    "configs",
    "eval_configs.txt"
  ),
  sep = "\t",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

dataset_info = read_benchmark_dataset_info(repo_root)

eval_sample_col_by_truth = c(
  "meancpm" = "n_eval_samples_meancpm",
  "sumcount_cpm" = "n_eval_samples_sumcount_cpm"
)

benchmark_root = file.path(repo_root, "Benchmarking_obj")
datasets = sort(list.dirs(
  benchmark_root,
  recursive = FALSE,
  full.names = FALSE
))
datasets = datasets[dir.exists(file.path(
  benchmark_root,
  datasets,
  "gene_prioritization"
))]

if (!is.null(dataset_include)) {
  datasets = datasets[datasets %in% dataset_include]
}

if (length(datasets) == 0L) {
  stop("No datasets selected")
}


## ---------------------------------------------------------------------------
## Summary helpers
## ---------------------------------------------------------------------------

empty_summary_table = function() {
  data.frame(
    dataset = character(),
    method = character(),
    cell_type_raw = character(),
    cell_type = character(),
    group = character(),
    top_n = integer(),
    ranking_source = character(),
    ranking_score = character(),
    comparison_score = character(),
    ranking_label = character(),
    avg_cor = numeric(),
    avg_cor_with_NA_penalty = numeric(),
    n_candidate_genes = integer(),
    n_selected_genes = integer(),
    n_cor_available_genes = integer(),
    n_evaluated_genes = integer(),
    n_constant_genes = integer(),
    n_eval_samples = integer(),
    tissue = character(),
    ct_abundance = numeric(),
    config_id = character(),
    config_slug = character(),
    truth_type = character(),
    indep_ref = character(),
    stringsAsFactors = FALSE
  )
}

rank_score_vector = function(
  score_values,
  candidate_genes,
  decreasing = TRUE
) {
  score_values = score_values[candidate_genes]
  score_values = score_values[is.finite(score_values)]

  if (length(score_values) == 0L) {
    return(character())
  }

  gene_names = names(score_values)
  score_order = if (decreasing) {
    order(-score_values, gene_names)
  } else {
    order(score_values, gene_names)
  }

  gene_names[score_order]
}

summarize_one_ranking = function(
  ranked_genes,
  ranking_source,
  ranking_score,
  ranking_label,
  cor_values,
  comparison_score,
  acceptable_mask_values,
  n_candidate_genes,
  dataset_name,
  method_name,
  cell_type_raw,
  cell_type,
  n_eval_samples,
  tissue,
  ct_abundance,
  config_id,
  config_slug_value,
  truth_type,
  indep_ref
) {
  rows = lapply(top_n_values, function(top_n) {
    n_selected = min(top_n, length(ranked_genes))
    selected_genes = if (n_selected == 0L) {
      character()
    } else {
      ranked_genes[seq_len(n_selected)]
    }

    cor_available_genes = selected_genes[
      selected_genes %in% names(cor_values)
    ]

    acceptable = acceptable_mask_values[cor_available_genes]
    acceptable[is.na(acceptable)] = FALSE
    evaluated_genes = cor_available_genes[!acceptable]
    values = cor_values[evaluated_genes]

    finite_values = values[is.finite(values)]
    avg_cor = if (length(finite_values) == 0L) {
      NA_real_
    } else {
      mean(finite_values)
    }

    values_with_penalty = values
    values_with_penalty[!is.finite(values_with_penalty)] = 0
    avg_cor_with_NA_penalty = if (
      length(values_with_penalty) == 0L
    ) {
      NA_real_
    } else {
      mean(values_with_penalty)
    }

    data.frame(
      dataset = dataset_name,
      method = method_name,
      cell_type_raw = cell_type_raw,
      cell_type = cell_type,
      group = paste0("top_", top_n),
      top_n = top_n,
      ranking_source = ranking_source,
      ranking_score = ranking_score,
      comparison_score = comparison_score,
      ranking_label = ranking_label,
      avg_cor = avg_cor,
      avg_cor_with_NA_penalty = avg_cor_with_NA_penalty,
      n_candidate_genes = n_candidate_genes,
      n_selected_genes = length(selected_genes),
      n_cor_available_genes = length(cor_available_genes),
      n_evaluated_genes = length(evaluated_genes),
      n_constant_genes = sum(!is.finite(values)),
      n_eval_samples = n_eval_samples,
      tissue = tissue,
      ct_abundance = ct_abundance,
      config_id = config_id,
      config_slug = config_slug_value,
      truth_type = truth_type,
      indep_ref = indep_ref,
      stringsAsFactors = FALSE
    )
  })

  do.call(rbind, rows)
}

cell_type_metadata = function(
  dataset_name,
  cell_type_raw,
  cell_type,
  truth_type
) {
  x = dataset_info[
    dataset_info$dataset == dataset_name &
      (
        dataset_info$cell_type_raw == cell_type_raw |
          dataset_info$cell_type == cell_type
      ),
    ,
    drop = FALSE
  ]

  if (nrow(x) == 0L) {
    return(list(
      n_eval_samples = NA_integer_,
      tissue = NA_character_,
      ct_abundance = NA_real_
    ))
  }

  eval_sample_col = unname(eval_sample_col_by_truth[truth_type])
  n_eval_samples = if (
    is.na(eval_sample_col) ||
      !eval_sample_col %in% colnames(x)
  ) {
    NA_integer_
  } else {
    as.integer(x[[eval_sample_col]][1])
  }

  list(
    n_eval_samples = n_eval_samples,
    tissue = as.character(x$tissue[1]),
    ct_abundance = as.numeric(x$ct_abundance[1])
  )
}


## ---------------------------------------------------------------------------
## Build dataset-level summary tables
## ---------------------------------------------------------------------------

summary_list = setNames(vector("list", length(datasets)), datasets)
availability_rows = list()

for (dataset_name in datasets) {
  message("Dataset: ", dataset_name)
  dataset_rows = list()

  for (config_id in config_include) {
    config_row = deconv_configs[
      deconv_configs$config_id == config_id,
      ,
      drop = FALSE
    ]

    if (nrow(config_row) != 1L) {
      stop(
        "Expected exactly one deconvolution config for ",
        config_id
      )
    }

    if (config_row$refType[[1]] != "indep") {
      stop(
        "This summary requires independent-reference configs: ",
        config_id
      )
    }

    eval_row = eval_configs[
      eval_configs$config_id == config_id,
      ,
      drop = FALSE
    ]
    if (nrow(eval_row) != 1L) {
      stop(
        "Expected exactly one evaluation config for ",
        config_id
      )
    }

    truth_type = eval_row$expected_truth_type[[1]]
    paths = resolve_deconv_paths(
      dataset_name,
      config_id,
      repo_root = repo_root
    )
    config_slug_value = config_slug(config_row)

    performance_dir = file.path(
      paths$obj_dir,
      "deconv_performance",
      paste0(
        config_slug_value,
        "__truth-",
        truth_type
      )
    )
    cor_dir = file.path(performance_dir, "spearman_cor")

    if (!dir.exists(cor_dir)) {
      message("  ", config_id, ": no Spearman directory")
      next
    }

    cor_files = sort(list.files(
      cor_dir,
      pattern = "\\.txt$",
      full.names = TRUE
    ))
    methods = sub("\\.txt$", "", basename(cor_files))

    if (!is.null(method_include)) {
      keep = methods %in% method_include
      cor_files = cor_files[keep]
      methods = methods[keep]
    }

    if (length(methods) == 0L) {
      next
    }

    mask = read_cor_na_acceptable_mask(
      dataset_name,
      truth_type,
      repo_root = repo_root
    )
    bulk_genes = read_bulk_gene_names(
      dataset_name,
      config_row$bulk_input[[1]],
      repo_root = repo_root
    )

    ref_dir = resolve_indep_ref_dir(
      dataset_name,
      repo_root = repo_root
    )
    indep_ref = basename(ref_dir)
    limma_stats = read_limma_top_genes_csv(ref_dir)
    indep_mapping = read_summary_indep_mapping(
      dataset_name,
      indep_ref,
      repo_root = repo_root
    )

    genesigtest_path =
      gene_prioritization_genesigtest_output_path(paths)
    genesigtest_available =
      include_genesigtest && file.exists(genesigtest_path)
    genesigtest_mat = if (genesigtest_available) {
      read_summary_numeric_matrix(genesigtest_path, required = TRUE)
    } else {
      NULL
    }

    if (include_genesigtest && !genesigtest_available) {
      message("  ", config_id, ": GeneSigTest result unavailable")
    }

    n_methods_available = 0L

    for (method_index in seq_along(methods)) {
      method_name = methods[method_index]
      cor_path = cor_files[method_index]

      score_paths = gene_prioritization_method_output_paths(
        paths,
        method_name
      )
      posthoc_paths = unlist(
        score_paths[posthoc_scores],
        use.names = TRUE
      )
      posthoc_paths = posthoc_paths[file.exists(posthoc_paths)]

      if (length(posthoc_paths) == 0L && !genesigtest_available) {
        next
      }

      cor_mat = read_summary_numeric_matrix(
        cor_path,
        required = TRUE
      )
      posthoc_mats = lapply(
        posthoc_paths,
        read_summary_numeric_matrix,
        required = TRUE
      )
      names(posthoc_mats) = names(posthoc_paths)
      ranking_mats = posthoc_mats
      if (genesigtest_available) {
        ranking_mats[["genesigtest_adapted_fdr"]] = genesigtest_mat
      }

      cell_types_raw = colnames(cor_mat)
      if (length(cell_types_raw) == 0L) {
        next
      }

      n_methods_available = n_methods_available + 1L
      cell_types = standardize_summary_cell_types(
        dataset_name,
        cell_types_raw,
        dataset_info
      )

      for (cell_type_index in seq_along(cell_types_raw)) {
        cell_type_raw = cell_types_raw[cell_type_index]
        cell_type = cell_types[cell_type_index]

        limma_cell_type = resolve_limma_cell_type(
          cell_type_raw = cell_type_raw,
          ref_type = "indep",
          limma_stats = limma_stats,
          indep_mapping = indep_mapping
        )
        if (is.na(limma_cell_type)) {
          next
        }

        cell_type_ranking_mats = ranking_mats[
          vapply(
            ranking_mats,
            function(x) cell_type_raw %in% colnames(x),
            logical(1)
          )
        ]
        if (length(cell_type_ranking_mats) == 0L) {
          next
        }

        ranking_values = lapply(cell_type_ranking_mats, function(x) {
          values = x[, cell_type_raw]
          names(values) = rownames(x)
          values
        })

        limma_values = limma_stats[, limma_cell_type]
        names(limma_values) = rownames(limma_stats)

        limma_finite_genes = names(limma_values)[
          is.finite(limma_values)
        ]

        comparison_data = lapply(
          names(ranking_values),
          function(comparison_score) {
            comparison_values =
              ranking_values[[comparison_score]]
            comparison_finite_genes = names(comparison_values)[
              is.finite(comparison_values)
            ]
            candidate_genes = Reduce(
              intersect,
              list(
                comparison_finite_genes,
                limma_finite_genes,
                bulk_genes
              )
            )

            pair_values = list(
              "indep_ref_limma" = limma_values
            )
            pair_values[[comparison_score]] = comparison_values
            ranked_genes = lapply(
              names(pair_values),
              function(ranking_score) {
                rank_score_vector(
                  pair_values[[ranking_score]],
                  candidate_genes,
                  decreasing = unname(
                    ranking_decreasing_by_score[ranking_score]
                  )
                )
              }
            )
            names(ranked_genes) = names(pair_values)

            list(
              candidate_genes = candidate_genes,
              ranked_genes = ranked_genes
            )
          }
        )
        names(comparison_data) = names(ranking_values)

        cor_values = cor_mat[, cell_type_raw]
        names(cor_values) = rownames(cor_mat)

        acceptable_mask_values = rep(
          FALSE,
          length(cor_values)
        )
        names(acceptable_mask_values) = names(cor_values)

        if (cell_type_raw %in% colnames(mask)) {
          common_mask_genes = intersect(
            names(cor_values),
            rownames(mask)
          )
          acceptable_mask_values[common_mask_genes] =
            mask[common_mask_genes, cell_type_raw] == 1
        }

        metadata = cell_type_metadata(
          dataset_name,
          cell_type_raw,
          cell_type,
          truth_type
        )

        comparison_rows = lapply(
          names(comparison_data),
          function(comparison_score) {
            comparison = comparison_data[[comparison_score]]
            ranking_rows = lapply(
              names(comparison$ranked_genes),
              function(ranking_score) {
                summarize_one_ranking(
                  ranked_genes =
                    comparison$ranked_genes[[ranking_score]],
                  ranking_source = unname(ranking_source_by_score[
                    ranking_score
                  ]),
                  ranking_score = ranking_score,
                  comparison_score = comparison_score,
                  ranking_label = unname(ranking_score_labels[
                    ranking_score
                  ]),
                  cor_values = cor_values,
                  acceptable_mask_values = acceptable_mask_values,
                  n_candidate_genes =
                    length(comparison$candidate_genes),
                  dataset_name = dataset_name,
                  method_name = method_name,
                  cell_type_raw = cell_type_raw,
                  cell_type = cell_type,
                  n_eval_samples = metadata$n_eval_samples,
                  tissue = metadata$tissue,
                  ct_abundance = metadata$ct_abundance,
                  config_id = config_id,
                  config_slug_value = config_slug_value,
                  truth_type = truth_type,
                  indep_ref = indep_ref
                )
              }
            )
            do.call(rbind, ranking_rows)
          }
        )

        dataset_rows[[length(dataset_rows) + 1L]] =
          do.call(rbind, comparison_rows)
      }
    }

    availability_rows[[length(availability_rows) + 1L]] =
      data.frame(
        dataset = dataset_name,
        config_id = config_id,
        posthoc_scores = paste(posthoc_scores, collapse = ","),
        genesigtest_available = genesigtest_available,
        n_methods_available = n_methods_available,
        stringsAsFactors = FALSE
      )

    message(
      "  ",
      config_id,
      ": methods with available rankings = ",
      n_methods_available,
      "; GeneSigTest = ",
      if (genesigtest_available) "available" else "unavailable"
    )
  }

  summary_list[[dataset_name]] = if (
    length(dataset_rows) == 0L
  ) {
    empty_summary_table()
  } else {
    x = do.call(rbind, dataset_rows)
    rownames(x) = NULL
    x
  }
}

summary_df = if (length(summary_list) == 0L) {
  empty_summary_table()
} else {
  do.call(rbind, summary_list)
}
rownames(summary_df) = NULL

availability_summary = if (length(availability_rows) == 0L) {
  data.frame(
    dataset = character(),
    config_id = character(),
    posthoc_scores = character(),
    genesigtest_available = logical(),
    n_methods_available = integer(),
    stringsAsFactors = FALSE
  )
} else {
  do.call(rbind, availability_rows)
}
rownames(availability_summary) = NULL


## ---------------------------------------------------------------------------
## Save reusable summary
## ---------------------------------------------------------------------------

output_path = file.path(
  repo_root,
  "benchmark_summary",
  "gene_prioritization_summary",
  output_file
)

saveRDS(summary_list, output_path)

message("Saved summary list: ", output_path)
message("Rows: ", nrow(summary_df))
