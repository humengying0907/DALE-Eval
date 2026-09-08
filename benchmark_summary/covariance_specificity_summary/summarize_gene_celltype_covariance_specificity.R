# Gene-by-cell-type covariance-specificity summary for one dataset or all
# datasets with pairwise-correlation results for the selected configuration.
#
# The operational definition uses the median absolute Pearson correlation
# between one focal cell type and its partner cell types. For each gene and
# focal cell type, the method and truth summaries use exactly the same finite
# partner-pair values.
#
# Notebook-style usage:
# 1. Edit the settings below.
# 2. Run from the repo root or this summary directory.
# 3. Inspect dataset_names, gene_universe_list,
#    gene_celltype_specificity, gene_coverage, and
#    method_celltype_summary.
# 4. The final section writes detailed and method-by-cell-type RDS results.

library(dplyr)


## settings ----

# dataset_include = "PBMC_1k1k"
dataset_include = "all"

config_include = "config01"

# NULL includes every evaluated method available for each selected dataset.
method_include = NULL

# TRUE requires every partner cell type to have finite truth and method
# correlations before a gene-by-cell-type specificity score is calculated.
# FALSE requires at least one matched finite partner pair and records the
# actual number used in n_matched_partner_pairs.
require_all_partner_pairs = TRUE


## paths and dataset selection ----

repo_root_candidates = unique(normalizePath(
  c(getwd(), file.path(getwd(), ".."), file.path(getwd(), "../..")),
  mustWork = FALSE
))

repo_root = repo_root_candidates[file.exists(file.path(
  repo_root_candidates,
  "DALE_Eval",
  "configs",
  "benchmark_ref_assignment.txt"
))][1]

if (is.na(repo_root)) {
  stop("Could not locate the ctse_benchmark repository root")
}

source(file.path(repo_root, "DALE_Eval", "modules", "config_helpers.R"))
source(file.path(
  repo_root,
  "DALE_Eval",
  "modules",
  "celltype_pairwise_cor_helpers.R"
))

config_row = read_deconv_config(config_include, repo_root = repo_root)
eval_config = read_pairwise_eval_config(
  config_include,
  repo_root = repo_root
)
expected_truth_type = eval_config$expected_truth_type[[1]]

performance_slug = paste0(
  config_slug(config_row),
  "__truth-",
  expected_truth_type
)

benchmark_root = file.path(repo_root, "Benchmarking_obj")
output_dir = file.path(
  repo_root,
  "benchmark_summary",
  "covariance_specificity_summary"
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

input_paths_for_dataset = function(dataset_name) {
  list(
    method_dir = file.path(
      benchmark_root,
      dataset_name,
      "deconv_performance",
      performance_slug,
      "celltype_pairwise_pearson_cor"
    ),
    truth_path = file.path(
      benchmark_root,
      dataset_name,
      "deconv_performance",
      paste0("truth_", expected_truth_type),
      "celltype_pairwise_pearson_cor",
      "Z_truth.txt"
    )
  )
}

dataset_has_pairwise_inputs = function(dataset_name) {
  paths = input_paths_for_dataset(dataset_name)
  dir.exists(paths$method_dir) && file.exists(paths$truth_path)
}

dataset_candidates = list.files(
  benchmark_root,
  full.names = FALSE,
  no.. = TRUE
)
dataset_candidates = sort(dataset_candidates[dir.exists(file.path(
  benchmark_root,
  dataset_candidates
))])

if (identical(dataset_include, "all")) {
  dataset_names = dataset_candidates[vapply(
    dataset_candidates,
    dataset_has_pairwise_inputs,
    logical(1)
  )]
} else {
  dataset_names = dataset_include
}

missing_dataset_inputs = dataset_names[!vapply(
  dataset_names,
  dataset_has_pairwise_inputs,
  logical(1)
)]

if (length(missing_dataset_inputs) > 0) {
  stop(
    "Missing pairwise method or truth inputs for config ",
    config_include,
    ": ",
    paste(missing_dataset_inputs, collapse = ", ")
  )
}

if (length(dataset_names) == 0) {
  stop("No datasets have complete pairwise inputs for config ", config_include)
}


## summarize one dataset ----

summarize_dataset_specificity = function(dataset_name) {
  input_paths = input_paths_for_dataset(dataset_name)

  indep_ref_dir = resolve_indep_ref_dir(
    dataset_name,
    repo_root = repo_root
  )
  indep_ref_signature_path = file.path(
    indep_ref_dir,
    "rowMeans_sig.csv"
  )

  truth_cor = read_pairwise_ctse_matrix(
    input_paths$truth_path,
    label = paste0(
      dataset_name,
      " truth ",
      expected_truth_type,
      " pairwise correlation"
    )
  )

  indep_ref_signature = read.csv(
    indep_ref_signature_path,
    row.names = 1,
    check.names = FALSE
  )

  # Match Fig. 5c: truth genes intersected with genes represented in the
  # assigned independent-reference signature.
  gene_universe = intersect(
    rownames(truth_cor),
    rownames(indep_ref_signature)
  )

  truth_pair_table = parse_celltype_pair_names(colnames(truth_cor))
  truth_cell_types = unique(c(
    truth_pair_table$cell_type_a,
    truth_pair_table$cell_type_b
  ))

  method_paths = list.files(
    input_paths$method_dir,
    pattern = "\\.txt$",
    full.names = TRUE
  )
  names(method_paths) = sub("\\.txt$", "", basename(method_paths))

  if (length(method_paths) == 0) {
    stop(
      "No evaluated method pairwise-correlation files found for ",
      dataset_name
    )
  }

  if (!is.null(method_include)) {
    missing_methods = setdiff(method_include, names(method_paths))
    if (length(missing_methods) > 0) {
      stop(
        "Selected methods do not have pairwise-correlation results for ",
        dataset_name,
        ": ",
        paste(missing_methods, collapse = ", ")
      )
    }
    method_paths = method_paths[method_include]
  }

  summarize_method_specificity = function(method_name) {
    method_cor = read_pairwise_ctse_matrix(
      method_paths[[method_name]],
      label = paste(dataset_name, method_name, "pairwise correlation")
    )

    pair_keep = intersect(colnames(truth_cor), colnames(method_cor))
    if (length(pair_keep) == 0) {
      stop(
        dataset_name,
        " ",
        method_name,
        " has no cell-type pairs shared with truth"
      )
    }

    method_pair_table = parse_celltype_pair_names(pair_keep)
    method_cell_types = unique(c(
      method_pair_table$cell_type_a,
      method_pair_table$cell_type_b
    ))
    missing_truth_cell_types = setdiff(
      truth_cell_types,
      method_cell_types
    )

    gene_keep = intersect(gene_universe, rownames(method_cor))

    cell_type_results = lapply(
      method_cell_types,
      function(focal_cell_type) {
        partner_pairs = celltype_pair_columns(
          pair_keep,
          focal_cell_type
        )

        truth_x = truth_cor[
          gene_keep,
          partner_pairs,
          drop = FALSE
        ]
        method_x = method_cor[
          gene_keep,
          partner_pairs,
          drop = FALSE
        ]

        matched_finite = is.finite(truth_x) & is.finite(method_x)
        n_matched_partner_pairs = rowSums(matched_finite)

        n_required_pairs = if (require_all_partner_pairs) {
          length(partner_pairs)
        } else {
          1L
        }
        score_keep = n_matched_partner_pairs >= n_required_pairs
        score_index = which(score_keep)

        truth_x[!matched_finite] = NA_real_
        method_x[!matched_finite] = NA_real_

        if (length(score_index) > 0) {
          truth_median_abs_cor = apply(
            abs(truth_x[score_index, , drop = FALSE]),
            1,
            median,
            na.rm = TRUE
          )
          method_median_abs_cor = apply(
            abs(method_x[score_index, , drop = FALSE]),
            1,
            median,
            na.rm = TRUE
          )
        } else {
          truth_median_abs_cor = numeric()
          method_median_abs_cor = numeric()
        }

        truth_raw_specificity = 1 - truth_median_abs_cor
        method_raw_specificity = 1 - method_median_abs_cor

        # A truth-matching method scores 1. Methods with lower absolute
        # correlation than truth are also capped at 1. When truth has absolute
        # correlation 1, no method can have a larger value, so its relative
        # score is set to 1.
        relative_covariance_specificity = ifelse(
          truth_raw_specificity > 0,
          pmin(
            1,
            pmax(0, method_raw_specificity / truth_raw_specificity)
          ),
          1
        )

        gene_df = data.frame(
          dataset = rep(dataset_name, length(score_index)),
          config_id = rep(config_include, length(score_index)),
          truth_type = rep(expected_truth_type, length(score_index)),
          method = rep(method_name, length(score_index)),
          gene = gene_keep[score_index],
          focal_cell_type = rep(
            focal_cell_type,
            length(score_index)
          ),
          n_partner_cell_types = rep(
            length(partner_pairs),
            length(score_index)
          ),
          n_matched_partner_pairs = n_matched_partner_pairs[score_index],
          truth_median_abs_cor = truth_median_abs_cor,
          method_median_abs_cor = method_median_abs_cor,
          truth_raw_specificity = truth_raw_specificity,
          method_raw_specificity = method_raw_specificity,
          delta_abs_cor = method_median_abs_cor - truth_median_abs_cor,
          excess_abs_cor = pmax(
            0,
            method_median_abs_cor - truth_median_abs_cor
          ),
          relative_covariance_specificity =
            relative_covariance_specificity,
          stringsAsFactors = FALSE
        )

        coverage_df = data.frame(
          dataset = dataset_name,
          config_id = config_include,
          truth_type = expected_truth_type,
          method = method_name,
          focal_cell_type = focal_cell_type,
          n_truth_cell_types = length(truth_cell_types),
          n_method_cell_types = length(method_cell_types),
          missing_truth_cell_types = paste(
            missing_truth_cell_types,
            collapse = ";"
          ),
          n_truth_celltype_pairs = ncol(truth_cor),
          n_method_celltype_pairs = length(pair_keep),
          n_truth_indep_ref_genes = length(gene_universe),
          n_method_genes_in_universe = length(gene_keep),
          n_partner_cell_types = length(partner_pairs),
          n_genes_at_least_one_pair = sum(
            n_matched_partner_pairs >= 1L
          ),
          n_genes_all_partner_pairs = sum(
            n_matched_partner_pairs == length(partner_pairs)
          ),
          n_genes_scored = sum(score_keep),
          gene_score_coverage = sum(score_keep) / length(gene_universe),
          stringsAsFactors = FALSE
        )

        list(gene = gene_df, coverage = coverage_df)
      }
    )

    list(
      gene = bind_rows(lapply(cell_type_results, `[[`, "gene")),
      coverage = bind_rows(lapply(cell_type_results, `[[`, "coverage"))
    )
  }

  method_results = lapply(
    names(method_paths),
    summarize_method_specificity
  )
  names(method_results) = names(method_paths)

  list(
    gene_universe = gene_universe,
    gene = bind_rows(lapply(method_results, `[[`, "gene")),
    coverage = bind_rows(lapply(method_results, `[[`, "coverage"))
  )
}


## selected-dataset results ----

dataset_results = lapply(
  dataset_names,
  summarize_dataset_specificity
)
names(dataset_results) = dataset_names

gene_universe_list = lapply(
  dataset_results,
  `[[`,
  "gene_universe"
)

gene_celltype_specificity = bind_rows(lapply(
  dataset_results,
  `[[`,
  "gene"
))

gene_coverage = bind_rows(lapply(
  dataset_results,
  `[[`,
  "coverage"
))


## method-by-cell-type summary ----

gene_score_summary = gene_celltype_specificity %>%
  group_by(
    dataset,
    config_id,
    truth_type,
    method,
    focal_cell_type
  ) %>%
  summarise(
    median_n_matched_partner_pairs = median(n_matched_partner_pairs),
    median_truth_abs_cor = median(truth_median_abs_cor),
    median_method_abs_cor = median(method_median_abs_cor),
    median_delta_abs_cor = median(delta_abs_cor),
    median_excess_abs_cor = median(excess_abs_cor),
    fraction_no_worse_than_truth = mean(delta_abs_cor <= 0),
    median_relative_covariance_specificity = median(
      relative_covariance_specificity
    ),
    mean_relative_covariance_specificity = mean(
      relative_covariance_specificity
    ),
    .groups = "drop"
  )

method_celltype_summary = gene_coverage %>%
  left_join(
    gene_score_summary,
    by = c(
      "dataset",
      "config_id",
      "truth_type",
      "method",
      "focal_cell_type"
    )
  ) %>%
  mutate(n_genes = n_genes_scored) %>%
  relocate(n_genes, .after = focal_cell_type) %>%
  arrange(dataset, focal_cell_type, method)

rm(
  repo_root_candidates, config_row, eval_config, performance_slug,
  benchmark_root, input_paths_for_dataset, dataset_has_pairwise_inputs,
  dataset_candidates, missing_dataset_inputs, dataset_results,
  gene_score_summary
)


## export ----

saveRDS(
  gene_celltype_specificity,
  file = file.path(
    output_dir,
    paste0(config_include, "__gene_celltype_specificity.RDS")
  )
)

saveRDS(
  method_celltype_summary,
  file = file.path(
    output_dir,
    paste0(
      config_include,
      "__gene_celltype_specificity_method_summary.RDS"
    )
  )
)
