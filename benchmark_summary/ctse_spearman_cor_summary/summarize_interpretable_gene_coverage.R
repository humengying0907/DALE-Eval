# Summarize interpretable CTSE gene coverage over a reference-aware universe.
#
# Notebook-style workflow:
#   1. Edit dataset_include, config_include, or method_include below.
#   2. Run from the repository root or
#      benchmark_summary/ctse_spearman_cor_summary/.
#   3. Inspect coverage_summary_list and coverage_summary.
#   4. The final section writes
#      interpretable_gene_coverage_summary_list.RDS.
#
# Definition for one dataset, truth type, and cell type:
#
#   truth/reference/variable universe =
#     truth genes
#     intersect assigned independent-reference genes
#     intersect genes evaluable in the truth-side acceptable-NA mask
#
#   interpretable method genes = genes in that universe with a finite
#     method Spearman correlation
#
#   interpretable gene coverage =
#     n_interpretable_genes / n_truth_indep_ref_variable_genes
#
# Missing method genes and present genes with non-finite correlations both
# reduce coverage.


## settings ----

# Set to NULL to include every dataset in the current Spearman summary.
dataset_include = NULL

# This reference-aware summary is intended for independent-reference configs.
# The expected truth type is obtained from eval_configs.txt for each config.
config_include = c("config01", "config02")

# Set to NULL to retain every method in the current all_genes summary rows.
method_include = NULL


## repository inputs ----

repo_root_candidates = unique(normalizePath(
  c(getwd(), file.path(getwd(), ".."), file.path(getwd(), "../..")),
  mustWork = FALSE
))

repo_root = repo_root_candidates[file.exists(file.path(
  repo_root_candidates,
  "benchmark_summary",
  "ctse_spearman_cor_summary",
  "spearman_by_top_genes_summary_list.RDS"
))][1]

if (is.na(repo_root)) {
  stop("Could not locate the ctse_benchmark repository root")
}

source(file.path(repo_root, "DALE_Eval", "modules", "config_helpers.R"))
source(file.path(
  repo_root,
  "DALE_Eval",
  "modules",
  "performance_summary_helpers.R"
))

spearman_summary_path = file.path(
  repo_root,
  "benchmark_summary",
  "ctse_spearman_cor_summary",
  "spearman_by_top_genes_summary_list.RDS"
)

coverage_output_path = file.path(
  repo_root,
  "benchmark_summary",
  "ctse_spearman_cor_summary",
  "interpretable_gene_coverage_summary_list.RDS"
)

deconv_configs_path = file.path(
  repo_root,
  "DALE_Eval",
  "configs",
  "deconv_configs.txt"
)

eval_configs_path = file.path(
  repo_root,
  "DALE_Eval",
  "configs",
  "eval_configs.txt"
)

ref_assignment_path = file.path(
  repo_root,
  "DALE_Eval",
  "configs",
  "benchmark_ref_assignment.txt"
)

spearman_summary_list = readRDS(spearman_summary_path)
dataset_info = read_benchmark_dataset_info(repo_root)

deconv_configs = read.delim(
  deconv_configs_path,
  sep = "\t",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

eval_configs = read.delim(
  eval_configs_path,
  sep = "\t",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

ref_assignment = read.delim(
  ref_assignment_path,
  sep = "\t",
  stringsAsFactors = FALSE,
  check.names = FALSE
)


## selected datasets and config contracts ----

required_summary_columns = c(
  "method",
  "cell_type",
  "group",
  "config_id",
  "truth_type"
)

config_contract = merge(
  deconv_configs,
  eval_configs[, c("config_id", "expected_truth_type")],
  by = "config_id",
  all = FALSE,
  sort = FALSE
)

config_contract = config_contract[
  config_contract$refType == "indep",
  ,
  drop = FALSE
]

if (!is.null(config_include)) {
  missing_configs = setdiff(config_include, config_contract$config_id)
  if (length(missing_configs) > 0L) {
    stop(
      "Requested configs are missing or are not independent-reference configs: ",
      paste(missing_configs, collapse = ", ")
    )
  }

  config_contract = config_contract[
    match(config_include, config_contract$config_id),
    ,
    drop = FALSE
  ]
}

datasets = names(spearman_summary_list)

if (!is.null(dataset_include)) {
  missing_datasets = setdiff(dataset_include, datasets)
  if (length(missing_datasets) > 0L) {
    stop(
      "Datasets missing from the Spearman summary: ",
      paste(missing_datasets, collapse = ", ")
    )
  }
  datasets = dataset_include
}


## small helpers ----

empty_coverage_table = function() {
  data.frame(
    dataset = character(),
    tissue = character(),
    cell_type_raw = character(),
    cell_type = character(),
    method = character(),
    config_id = character(),
    truth_type = character(),
    refType = character(),
    n_truth_indep_ref_variable_genes = integer(),
    n_method_genes_in_universe = integer(),
    n_missing_method_genes = integer(),
    n_nonfinite_method_genes = integer(),
    n_interpretable_genes = integer(),
    interpretable_gene_coverage = numeric(),
    stringsAsFactors = FALSE
  )
}

read_independent_reference_genes = function(dataset_name) {
  assignment_row = ref_assignment[
    ref_assignment$dataset == dataset_name,
    ,
    drop = FALSE
  ]

  if (nrow(assignment_row) != 1L) {
    stop(
      dataset_name,
      ": expected exactly one assigned independent reference"
    )
  }

  ref_path = file.path(
    repo_root,
    "Indep_scReference",
    assignment_row$indep_ref[[1]],
    "rowMeans_sig.csv"
  )

  if (!file.exists(ref_path)) {
    stop("Independent-reference signature not found: ", ref_path)
  }

  ref_stub = read.csv(
    ref_path,
    row.names = 1,
    check.names = FALSE
  )

  unique(rownames(ref_stub))
}


## build reference-aware coverage ----

coverage_summary_list = setNames(
  vector("list", length(datasets)),
  datasets
)

for (dataset_name in datasets) {
  message("Dataset: ", dataset_name)

  summary_one = spearman_summary_list[[dataset_name]]
  if (nrow(summary_one) == 0L) {
    coverage_summary_list[[dataset_name]] = empty_coverage_table()
    next
  }

  missing_summary_columns = setdiff(
    required_summary_columns,
    colnames(summary_one)
  )
  if (length(missing_summary_columns) > 0L) {
    stop(
      dataset_name,
      ": Spearman summary missing columns: ",
      paste(missing_summary_columns, collapse = ", ")
    )
  }

  summary_one = merge(
    summary_one,
    config_contract[, c(
      "config_id",
      "expected_truth_type",
      "refType",
      "bulk_input",
      "bulk_normalization",
      "frac_input"
    )],
    by = "config_id",
    all = FALSE,
    sort = FALSE,
    suffixes = c("", "_contract")
  )

  summary_one = summary_one[
    summary_one$group == "all_genes" &
      summary_one$truth_type == summary_one$expected_truth_type,
    ,
    drop = FALSE
  ]

  if (!is.null(method_include)) {
    summary_one = summary_one[
      summary_one$method %in% method_include,
      ,
      drop = FALSE
    ]
  }

  if (nrow(summary_one) == 0L) {
    coverage_summary_list[[dataset_name]] = empty_coverage_table()
    next
  }

  summary_keys = paste(
    summary_one$config_id,
    summary_one$truth_type,
    summary_one$method,
    summary_one$cell_type,
    sep = "\r"
  )
  if (anyDuplicated(summary_keys)) {
    stop(dataset_name, ": duplicated all_genes summary rows")
  }

  dataset_annotation = dataset_info[
    dataset_info$dataset == dataset_name,
    c("cell_type_raw", "cell_type", "tissue"),
    drop = FALSE
  ]
  if (anyDuplicated(dataset_annotation$cell_type)) {
    stop(dataset_name, ": standardized cell types are not unique")
  }

  independent_reference_genes = read_independent_reference_genes(
    dataset_name
  )

  dataset_rows = list()

  selected_config_ids = unique(summary_one$config_id)
  for (config_id in selected_config_ids) {
    config_row = config_contract[
      config_contract$config_id == config_id,
      ,
      drop = FALSE
    ]
    truth_type = config_row$expected_truth_type[[1]]

    summary_config = summary_one[
      summary_one$config_id == config_id &
        summary_one$truth_type == truth_type,
      ,
      drop = FALSE
    ]

    acceptable_mask = read_cor_na_acceptable_mask(
      dataset = dataset_name,
      truth_type = truth_type,
      repo_root = repo_root
    )

    valid_mask = is.finite(acceptable_mask) &
      acceptable_mask %in% c(0, 1)
    if (!all(valid_mask)) {
      stop(
        dataset_name,
        " / ", truth_type,
        ": acceptable-NA mask contains values other than 0 or 1"
      )
    }

    performance_dir = file.path(
      repo_root,
      "Benchmarking_obj",
      dataset_name,
      "deconv_performance",
      paste0(config_slug(config_row), "__truth-", truth_type),
      "spearman_cor"
    )

    for (method_name in unique(summary_config$method)) {
      cor_path = file.path(
        performance_dir,
        paste0(method_name, ".txt")
      )
      method_cor = read_summary_numeric_matrix(cor_path, required = TRUE)

      summary_method = summary_config[
        summary_config$method == method_name,
        ,
        drop = FALSE
      ]

      for (row_index in seq_len(nrow(summary_method))) {
        source_row = summary_method[row_index, , drop = FALSE]
        cell_type = source_row$cell_type[[1]]
        annotation_row = dataset_annotation[
          dataset_annotation$cell_type == cell_type,
          ,
          drop = FALSE
        ]

        if (nrow(annotation_row) != 1L) {
          stop(
            dataset_name,
            " / ", cell_type,
            ": expected one benchmark dataset annotation row"
          )
        }

        cell_type_raw = annotation_row$cell_type_raw[[1]]
        if (!cell_type_raw %in% colnames(acceptable_mask)) {
          stop(
            dataset_name,
            " / ", truth_type,
            " / ", cell_type_raw,
            ": cell type missing from acceptable-NA mask"
          )
        }
        if (!cell_type_raw %in% colnames(method_cor)) {
          stop(
            dataset_name,
            " / ", config_id,
            " / ", method_name,
            " / ", cell_type_raw,
            ": cell type missing from method correlation matrix"
          )
        }

        truth_variable_genes = rownames(acceptable_mask)[
          acceptable_mask[, cell_type_raw] == 0
        ]
        truth_indep_ref_variable_genes = intersect(
          truth_variable_genes,
          independent_reference_genes
        )

        method_genes_in_universe = intersect(
          truth_indep_ref_variable_genes,
          rownames(method_cor)
        )
        method_values = method_cor[
          method_genes_in_universe,
          cell_type_raw
        ]

        n_truth_indep_ref_variable_genes = length(
          truth_indep_ref_variable_genes
        )
        n_method_genes_in_universe = length(method_genes_in_universe)
        n_nonfinite_method_genes = sum(!is.finite(method_values))
        n_interpretable_genes = sum(is.finite(method_values))
        n_missing_method_genes =
          n_truth_indep_ref_variable_genes - n_method_genes_in_universe

        interpretable_gene_coverage = if (
          n_truth_indep_ref_variable_genes > 0L
        ) {
          n_interpretable_genes / n_truth_indep_ref_variable_genes
        } else {
          NA_real_
        }

        if (
          is.finite(interpretable_gene_coverage) &&
            (
              interpretable_gene_coverage < 0 ||
                interpretable_gene_coverage > 1
            )
        ) {
          stop(
            dataset_name,
            " / ", config_id,
            " / ", method_name,
            " / ", cell_type,
            ": interpretable coverage is outside [0, 1]"
          )
        }


        dataset_rows[[length(dataset_rows) + 1L]] = data.frame(
          dataset = dataset_name,
          tissue = annotation_row$tissue[[1]],
          cell_type_raw = cell_type_raw,
          cell_type = cell_type,
          method = method_name,
          config_id = config_id,
          truth_type = truth_type,
          refType = config_row$refType[[1]],
          n_truth_indep_ref_variable_genes =
            n_truth_indep_ref_variable_genes,
          n_method_genes_in_universe = n_method_genes_in_universe,
          n_missing_method_genes = n_missing_method_genes,
          n_nonfinite_method_genes = n_nonfinite_method_genes,
          n_interpretable_genes = n_interpretable_genes,
          interpretable_gene_coverage = interpretable_gene_coverage,
          stringsAsFactors = FALSE
        )
      }
    }
  }

  coverage_summary_list[[dataset_name]] = if (
    length(dataset_rows) > 0L
  ) {
    result = do.call(rbind, dataset_rows)
    rownames(result) = NULL
    result
  } else {
    empty_coverage_table()
  }

  message(
    "  Coverage rows: ",
    nrow(coverage_summary_list[[dataset_name]])
  )
}

coverage_summary = do.call(rbind, coverage_summary_list)
rownames(coverage_summary) = NULL


## export ----

saveRDS(
  coverage_summary_list,
  file = coverage_output_path
)

message("Saved: ", coverage_output_path)
