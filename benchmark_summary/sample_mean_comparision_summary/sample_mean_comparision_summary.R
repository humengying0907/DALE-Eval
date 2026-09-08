# Config01 gene-level and intra-cell-type sample-mean comparison summaries.
#
# Notebook-style usage:
# 1. Edit dataset_include, method_include, or config_include if needed.
# 2. Run from the repository root or this summary directory.
# 3. Inspect dataset_result_list,
#    truth_weighted_ordering_score_summary_list, ccc_summary_list,
#    intra_cell_type_cor_summary, and method_status_summary.
# 4. The export section writes two RDS lists and one tab-separated text file.


## settings ----

dataset_include = c(
  "BRCA_Bassez2021",
  "CRC_Pelka2021",
  "LUAD_Kim2020",
  "PBMC_1k1k",
  "PBMC_Perez2022",
  "ROSMAP_AD430_Mathys2023",
  "ROSMAP_AD92_Xiong2023"
)

method_include = c(
  "InstaPrism",
  "EPICunmix",
  "ENIGMAtrace",
  "Unico",
  "TCA",
  "CIBERSORTx",
  "bMIND",
  "ENIGMAL2",
  "BLUE",
  "scTAPE"
)

config_include = "config01"

truth_weighted_ordering_score_filename =
  "truth_weighted_ordering_score_summary_list.RDS"
ccc_filename = "ccc_summary_list.RDS"
intra_cell_type_cor_filename = "intra_cell_type_cor_summary.txt"


## repository inputs ----

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

summary_dir = file.path(
  repo_root,
  "benchmark_summary",
  "sample_mean_comparision_summary"
)

source(file.path(repo_root, "DALE_Eval", "modules", "config_helpers.R"))

config_row = read_deconv_config(config_include, repo_root = repo_root)

eval_configs = read.delim(
  file.path(repo_root, "DALE_Eval", "configs", "eval_configs.txt"),
  sep = "\t",
  stringsAsFactors = FALSE,
  check.names = FALSE
)
eval_config = eval_configs[
  eval_configs$config_id == config_include,
  ,
  drop = FALSE
]
if (nrow(eval_config) != 1L) {
  stop("Expected exactly one evaluation config for ", config_include)
}

expected_truth_type = eval_config$expected_truth_type[[1]]
performance_slug = paste0(
  config_slug(config_row),
  "__truth-",
  expected_truth_type
)

dataset_info_all = read.delim(
  file.path(
    repo_root,
    "DALE_Eval",
    "configs",
    "benchmark_dataset_info.txt"
  ),
  sep = "\t",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

cell_type_mapping_all = read.delim(
  file.path(
    repo_root,
    "Indep_scReference",
    "cell_type_mapping.txt"
  ),
  sep = "\t",
  stringsAsFactors = FALSE,
  check.names = FALSE
)


## helpers ----

read_sample_mean_matrix = function(path) {
  x = as.matrix(read.delim(
    path,
    row.names = 1,
    stringsAsFactors = FALSE,
    check.names = FALSE
  ))
  suppressWarnings(storage.mode(x) <- "double")
  x
}

# For each gene, compare every matched cell-type pair. Truth differences
# determine the weights; an inferred tie contributes zero, a preserved order
# contributes positively, and a reversed order contributes negatively.
calculate_truth_weighted_ordering_score = function(
    truth_matrix,
    inferred_matrix
) {
  output = setNames(rep(NA_real_, nrow(truth_matrix)), rownames(truth_matrix))
  complete = rowSums(
    is.finite(truth_matrix) & is.finite(inferred_matrix)
  ) == ncol(truth_matrix)

  if (!any(complete)) {
    return(output)
  }

  truth_complete = truth_matrix[complete, , drop = FALSE]
  inferred_complete = inferred_matrix[complete, , drop = FALSE]
  pair_index = utils::combn(seq_len(ncol(truth_complete)), 2L)

  truth_difference =
    truth_complete[, pair_index[1L, ], drop = FALSE] -
    truth_complete[, pair_index[2L, ], drop = FALSE]
  inferred_difference =
    inferred_complete[, pair_index[1L, ], drop = FALSE] -
    inferred_complete[, pair_index[2L, ], drop = FALSE]

  truth_weight = abs(truth_difference)
  denominator = rowSums(truth_weight)
  numerator = rowSums(
    truth_weight *
      sign(truth_difference) *
      sign(inferred_difference)
  )

  score = rep(NA_real_, nrow(truth_complete))
  truth_has_ordering = is.finite(denominator) & denominator > 0
  score[truth_has_ordering] =
    numerator[truth_has_ordering] / denominator[truth_has_ordering]

  output[rownames(truth_complete)] = score
  output
}

# Lin's concordance correlation coefficient across matched cell types, using
# population covariance and variances for the finite paired profile.
calculate_ccc = function(truth_matrix, inferred_matrix) {
  output = setNames(rep(NA_real_, nrow(truth_matrix)), rownames(truth_matrix))
  complete = rowSums(
    is.finite(truth_matrix) & is.finite(inferred_matrix)
  ) == ncol(truth_matrix)

  if (!any(complete)) {
    return(output)
  }

  truth_complete = truth_matrix[complete, , drop = FALSE]
  inferred_complete = inferred_matrix[complete, , drop = FALSE]
  truth_mean = rowMeans(truth_complete)
  inferred_mean = rowMeans(inferred_complete)
  truth_centered = truth_complete - truth_mean
  inferred_centered = inferred_complete - inferred_mean

  truth_variance = rowMeans(truth_centered^2)
  inferred_variance = rowMeans(inferred_centered^2)
  covariance = rowMeans(truth_centered * inferred_centered)
  denominator =
    truth_variance +
    inferred_variance +
    (truth_mean - inferred_mean)^2

  ccc = rep(NA_real_, nrow(truth_complete))
  defined = is.finite(denominator) & denominator > 0
  ccc[defined] = 2 * covariance[defined] / denominator[defined]

  negligible_boundary_error =
    is.finite(ccc) & abs(ccc) > 1 & abs(ccc) <= 1 + 1e-12
  ccc[negligible_boundary_error] = pmax(
    -1,
    pmin(1, ccc[negligible_boundary_error])
  )
  ccc[is.finite(ccc) & abs(ccc) > 1] = NA_real_

  output[rownames(truth_complete)] = ccc
  output
}

safe_spearman = function(truth_values, inferred_values) {
  keep = is.finite(truth_values) & is.finite(inferred_values)
  if (
    sum(keep) < 2L ||
      length(unique(truth_values[keep])) < 2L ||
      length(unique(inferred_values[keep])) < 2L
  ) {
    return(NA_real_)
  }

  suppressWarnings(stats::cor(
    truth_values[keep],
    inferred_values[keep],
    method = "spearman"
  ))
}


## prepare one dataset ----

build_dataset_summary = function(dataset_name) {
  method_sample_mean_dir = file.path(
    repo_root,
    "Benchmarking_obj",
    dataset_name,
    "deconv_performance",
    performance_slug,
    "sample_mean"
  )
  truth_sample_mean_path = file.path(
    repo_root,
    "Benchmarking_obj",
    dataset_name,
    "deconv_performance",
    paste0("truth_", expected_truth_type),
    "sample_mean",
    "Z_truth.txt"
  )
  truth_sample_mean = read_sample_mean_matrix(truth_sample_mean_path)

  indep_ref_dir = resolve_indep_ref_dir(
    dataset_name,
    repo_root = repo_root
  )
  indep_ref_name = basename(indep_ref_dir)
  indep_ref_signature = read.csv(
    file.path(indep_ref_dir, "rowMeans_sig.csv"),
    row.names = 1,
    check.names = FALSE
  )

  dataset_mapping = cell_type_mapping_all[
    cell_type_mapping_all$dataset == dataset_name &
      cell_type_mapping_all$indep_ref == indep_ref_name,
    ,
    drop = FALSE
  ]
  mapped_cell_types = unique(dataset_mapping$target_cell_type[
    dataset_mapping$indep_ref_cell_type %in% colnames(indep_ref_signature)
  ])
  mapped_cell_types = mapped_cell_types[
    !is.na(mapped_cell_types) & nzchar(mapped_cell_types)
  ]
  matched_cell_types = intersect(
    colnames(truth_sample_mean),
    mapped_cell_types
  )

  abundance_rows = dataset_info_all[
    dataset_info_all$dataset == dataset_name &
      dataset_info_all$cell_type_raw %in% matched_cell_types,
    c("cell_type_raw", "ct_abundance")
  ]
  cell_type_meta = aggregate(
    ct_abundance ~ cell_type_raw,
    data = abundance_rows,
    FUN = sum,
    na.rm = TRUE
  )
  cell_type_meta = cell_type_meta[order(
    -cell_type_meta$ct_abundance,
    cell_type_meta$cell_type_raw
  ), ]
  missing_abundance = setdiff(
    matched_cell_types,
    cell_type_meta$cell_type_raw
  )
  if (length(missing_abundance) > 0L) {
    stop(
      dataset_name,
      " is missing truth-abundance annotations for: ",
      paste(missing_abundance, collapse = ", ")
    )
  }
  matched_cell_types = cell_type_meta$cell_type_raw

  if (length(matched_cell_types) < 2L) {
    stop(dataset_name, " has fewer than two matched cell types")
  }

  truth_indep_ref_genes = intersect(
    rownames(truth_sample_mean),
    rownames(indep_ref_signature)
  )

  method_paths = file.path(
    method_sample_mean_dir,
    paste0(method_include, ".txt")
  )
  names(method_paths) = method_include
  method_sample_mean_list = lapply(method_paths, function(path) {
    if (file.exists(path)) {
      read_sample_mean_matrix(path)
    } else {
      NULL
    }
  })

  ordering_score_matrix = matrix(
    NA_real_,
    nrow = length(truth_indep_ref_genes),
    ncol = length(method_include),
    dimnames = list(truth_indep_ref_genes, method_include)
  )
  ccc_matrix = ordering_score_matrix

  intra_rows = vector(
    "list",
    length(method_include) * length(matched_cell_types)
  )
  status_rows = vector("list", length(method_include))
  intra_index = 0L

  for (method_name in method_include) {
    method_sample_mean = method_sample_mean_list[[method_name]]
    method_file_available = !is.null(method_sample_mean)
    method_cell_types = if (method_file_available) {
      intersect(matched_cell_types, colnames(method_sample_mean))
    } else {
      character(0)
    }
    all_cell_types_available = identical(
      sort(method_cell_types),
      sort(matched_cell_types)
    )
    method_genes = if (method_file_available) {
      intersect(truth_indep_ref_genes, rownames(method_sample_mean))
    } else {
      character(0)
    }

    if (all_cell_types_available && length(method_genes) > 0L) {
      truth_profiles = truth_sample_mean[
        method_genes,
        matched_cell_types,
        drop = FALSE
      ]
      inferred_profiles = method_sample_mean[
        method_genes,
        matched_cell_types,
        drop = FALSE
      ]

      ordering_score_matrix[method_genes, method_name] =
        calculate_truth_weighted_ordering_score(
          truth_profiles,
          inferred_profiles
        )
      ccc_matrix[method_genes, method_name] = calculate_ccc(
        truth_profiles,
        inferred_profiles
      )
    }

    for (cell_type in matched_cell_types) {
      intra_index = intra_index + 1L
      n_genes = 0L
      intra_cell_type_cor = NA_real_

      if (
        method_file_available &&
          cell_type %in% colnames(method_sample_mean) &&
          length(method_genes) > 0L
      ) {
        truth_values = truth_sample_mean[method_genes, cell_type]
        inferred_values = method_sample_mean[method_genes, cell_type]
        complete = is.finite(truth_values) & is.finite(inferred_values)
        n_genes = sum(complete)
        intra_cell_type_cor = safe_spearman(
          truth_values,
          inferred_values
        )
      }

      intra_rows[[intra_index]] = data.frame(
        method = method_name,
        cell_type = cell_type,
        n_genes = n_genes,
        intra_cell_type_cor = intra_cell_type_cor,
        dataset = dataset_name,
        stringsAsFactors = FALSE
      )
    }

    status_rows[[match(method_name, method_include)]] = data.frame(
      dataset = dataset_name,
      method = method_name,
      method_file_available = method_file_available,
      all_matched_cell_types_available = all_cell_types_available,
      n_matched_cell_types = length(matched_cell_types),
      n_method_cell_types_available = length(method_cell_types),
      n_truth_indep_ref_genes = length(truth_indep_ref_genes),
      n_method_genes_available = length(method_genes),
      stringsAsFactors = FALSE
    )
  }

  list(
    truth_weighted_ordering_score = as.data.frame(
      ordering_score_matrix,
      check.names = FALSE
    ),
    ccc = as.data.frame(ccc_matrix, check.names = FALSE),
    intra_cell_type_cor = do.call(rbind, intra_rows),
    method_status = do.call(rbind, status_rows),
    matched_cell_types = matched_cell_types,
    truth_indep_ref_genes = truth_indep_ref_genes
  )
}


## combine dataset summaries ----

dataset_result_list = lapply(dataset_include, build_dataset_summary)
names(dataset_result_list) = dataset_include

truth_weighted_ordering_score_summary_list = lapply(
  dataset_result_list,
  `[[`,
  "truth_weighted_ordering_score"
)
ccc_summary_list = lapply(dataset_result_list, `[[`, "ccc")

intra_cell_type_cor_summary = do.call(rbind, lapply(
  dataset_result_list,
  `[[`,
  "intra_cell_type_cor"
))
rownames(intra_cell_type_cor_summary) = NULL

method_status_summary = do.call(rbind, lapply(
  dataset_result_list,
  `[[`,
  "method_status"
))
rownames(method_status_summary) = NULL

rm(
  repo_root_candidates, config_row, eval_configs, eval_config,
  dataset_info_all, cell_type_mapping_all
)


## export ----

dir.create(summary_dir, recursive = TRUE, showWarnings = FALSE)

saveRDS(
  truth_weighted_ordering_score_summary_list,
  file.path(summary_dir, truth_weighted_ordering_score_filename)
)
saveRDS(
  ccc_summary_list,
  file.path(summary_dir, ccc_filename)
)
write.table(
  intra_cell_type_cor_summary,
  file = file.path(summary_dir, intra_cell_type_cor_filename),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  na = "NA"
)
