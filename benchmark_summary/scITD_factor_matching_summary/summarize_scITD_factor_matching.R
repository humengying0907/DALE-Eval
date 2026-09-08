# Match independently fitted scITD factors to truth factors using loadings.
#
# Notebook-style usage:
# 1. Edit the dataset, config, or method settings if needed.
# 2. Run from the repository root or this summary directory.
# 3. Inspect factor_matching, loading_correlation_matrix_list, and
#    scITD_factor_matching_summary.
# 4. The export section writes one RDS containing the correlation matrices
#    and matched-factor tables.


## settings ----

dataset_include = c(
  "PBMC_Perez2022",
  "ROSMAP_AD430_Mathys2023"
)

config_include = c("config01", "config02")

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

truth_result_name = "truth_sumcount"
summary_rds_filename = "scITD_factor_matching_summary.RDS"


## repository inputs ----

repo_root_candidates = unique(normalizePath(
  c(getwd(), file.path(getwd(), ".."), file.path(getwd(), "../..")),
  mustWork = FALSE
))

repo_root = repo_root_candidates[file.exists(file.path(
  repo_root_candidates,
  "DALE_Eval",
  "configs",
  "scITD_configs.txt"
))][1]

if (is.na(repo_root)) {
  stop("Could not locate the ctse_benchmark repository root")
}

summary_dir = file.path(
  repo_root,
  "benchmark_summary",
  "scITD_factor_matching_summary"
)

source(file.path(repo_root, "DALE_Eval", "modules", "config_helpers.R"))

method_preprocessing = read.delim(
  file.path(
    repo_root,
    "DALE_Eval",
    "configs",
    "scITD_method_preprocessing.txt"
  ),
  stringsAsFactors = FALSE
)
libnorm_method_include = method_include[
  method_preprocessing$preprocess_mode[
    match(method_include, method_preprocessing$method)
  ] == "log1p"
]
result_include = c(
  method_include,
  paste0(libnorm_method_include, "_libnorm")
)

config_table = do.call(rbind, lapply(config_include, function(config_id) {
  config_row = read_deconv_config(config_id, repo_root = repo_root)
  data.frame(
    config_id = config_id,
    config_slug = config_slug(config_row),
    stringsAsFactors = FALSE
  )
}))

job_table = expand.grid(
  dataset = dataset_include,
  config_id = config_include,
  result_name = result_include,
  stringsAsFactors = FALSE
)
job_table$config_slug = config_table$config_slug[
  match(job_table$config_id, config_table$config_id)
]
job_table$truth_dir = file.path(
  repo_root,
  "Benchmarking_obj",
  job_table$dataset,
  "scITD_res",
  truth_result_name
)
job_table$method_dir = file.path(
  repo_root,
  "Benchmarking_obj",
  job_table$dataset,
  "scITD_res",
  job_table$config_slug,
  job_table$result_name
)


## helpers ----

read_scitd_matrix = function(path) {
  x = as.matrix(read.delim(
    path,
    row.names = 1,
    stringsAsFactors = FALSE,
    check.names = FALSE
  ))
  suppressWarnings(storage.mode(x) <- "double")
  x
}

scitd_loading_paths = function(result_dir) {
  paths = list.files(
    result_dir,
    pattern = "^gene_celltype_loading_Factor[0-9]+\\.txt$",
    full.names = TRUE
  )
  factor_number = as.integer(sub(
    ".*Factor([0-9]+)\\.txt$",
    "\\1",
    paths
  ))
  paths[order(factor_number)]
}

has_scitd_result = function(result_dir) {
  file.exists(file.path(result_dir, "sample_scores.txt")) &&
    length(scitd_loading_paths(result_dir)) > 0L
}

read_scitd_loadings = function(result_dir) {
  loading_paths = scitd_loading_paths(result_dir)
  factor_names = sub(
    "^gene_celltype_loading_(Factor[0-9]+)\\.txt$",
    "\\1",
    basename(loading_paths)
  )
  loadings = lapply(loading_paths, read_scitd_matrix)
  names(loadings) = factor_names
  loadings
}

factor_permutations = function(x) {
  if (length(x) == 1L) {
    return(matrix(x, nrow = 1L))
  }
  do.call(rbind, lapply(seq_along(x), function(index) {
    cbind(x[index], factor_permutations(x[-index]))
  }))
}

best_factor_assignment = function(correlation_matrix) {
  stopifnot(
    nrow(correlation_matrix) == ncol(correlation_matrix),
    all(is.finite(correlation_matrix))
  )

  permutations = factor_permutations(seq_len(ncol(correlation_matrix)))
  permutation_scores = apply(permutations, 1L, function(assignment) {
    sum(abs(correlation_matrix[cbind(
      seq_len(nrow(correlation_matrix)),
      assignment
    )]))
  })
  as.integer(permutations[which.max(permutation_scores), ])
}

compare_scitd_factors = function(job_row) {
  truth_loadings = read_scitd_loadings(job_row$truth_dir)
  method_loadings = read_scitd_loadings(job_row$method_dir)
  truth_scores = read_scitd_matrix(file.path(
    job_row$truth_dir,
    "sample_scores.txt"
  ))
  method_scores = read_scitd_matrix(file.path(
    job_row$method_dir,
    "sample_scores.txt"
  ))

  common_genes = intersect(
    rownames(truth_loadings[[1L]]),
    rownames(method_loadings[[1L]])
  )
  common_cell_types = intersect(
    colnames(truth_loadings[[1L]]),
    colnames(method_loadings[[1L]])
  )
  common_samples = intersect(rownames(truth_scores), rownames(method_scores))

  truth_loading_vectors = do.call(cbind, lapply(
    truth_loadings,
    function(x) as.vector(x[
      common_genes,
      common_cell_types,
      drop = FALSE
    ])
  ))
  method_loading_vectors = do.call(cbind, lapply(
    method_loadings,
    function(x) as.vector(x[
      common_genes,
      common_cell_types,
      drop = FALSE
    ])
  ))
  colnames(truth_loading_vectors) = names(truth_loadings)
  colnames(method_loading_vectors) = names(method_loadings)

  loading_correlation_matrix = cor(
    truth_loading_vectors,
    method_loading_vectors,
    method = "pearson"
  )
  matched_method_index = best_factor_assignment(
    loading_correlation_matrix
  )
  factor_index = seq_len(nrow(loading_correlation_matrix))
  matched_loading_cor = loading_correlation_matrix[cbind(
    factor_index,
    matched_method_index
  )]
  truth_factor = rownames(loading_correlation_matrix)
  matched_method_factor = colnames(loading_correlation_matrix)[
    matched_method_index
  ]
  sign_flip = ifelse(matched_loading_cor < 0, -1, 1)
  sample_score_pearson = vapply(factor_index, function(index) {
    cor(
      truth_scores[common_samples, truth_factor[index]],
      sign_flip[index] * method_scores[
        common_samples,
        matched_method_factor[index]
      ],
      method = "pearson"
    )
  }, numeric(1))

  factor_matches = data.frame(
    method = job_row$result_name,
    truth_factor = truth_factor,
    matched_method_factor = matched_method_factor,
    loading_pearson = matched_loading_cor,
    abs_loading_pearson = abs(matched_loading_cor),
    sign_flip = sign_flip,
    sample_score_pearson = sample_score_pearson,
    stringsAsFactors = FALSE
  )

  list(
    factor_matches = factor_matches,
    loading_correlation_matrix = loading_correlation_matrix
  )
}


## factor matching ----

available_job_table = job_table[
  vapply(job_table$truth_dir, has_scitd_result, logical(1)) &
    vapply(job_table$method_dir, has_scitd_result, logical(1)),
  ,
  drop = FALSE
]

comparison_list = lapply(seq_len(nrow(available_job_table)), function(index) {
  compare_scitd_factors(available_job_table[index, , drop = FALSE])
})
names(comparison_list) = paste(
  available_job_table$dataset,
  available_job_table$config_id,
  available_job_table$result_name,
  sep = "__"
)

loading_correlation_matrix_list = lapply(
  comparison_list,
  `[[`,
  "loading_correlation_matrix"
)

dataset_config_keys = unlist(lapply(dataset_include, function(dataset) {
  paste(dataset, config_include, sep = "__")
}))
available_dataset_config = paste(
  available_job_table$dataset,
  available_job_table$config_id,
  sep = "__"
)

factor_matching = setNames(lapply(dataset_config_keys, function(key) {
  matched_rows = do.call(rbind, lapply(
    comparison_list[available_dataset_config == key],
    `[[`,
    "factor_matches"
  ))
  rownames(matched_rows) = NULL
  matched_rows
}), dataset_config_keys)

scITD_factor_matching_summary = list(
  loading_correlation_matrix_list = loading_correlation_matrix_list,
  factor_matching = factor_matching
)


## export ----

dir.create(summary_dir, recursive = TRUE, showWarnings = FALSE)

saveRDS(
  scITD_factor_matching_summary,
  file.path(summary_dir, summary_rds_filename)
)
