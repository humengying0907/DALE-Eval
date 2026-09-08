# Summarize high-fraction sample retention and inferred cell-type abundance.
#
# Notebook-style workflow:
#   1. Edit dataset_include or config_include below.
#   2. Run from the repository root or
#      benchmark_summary/sample_filter_effect_summary/.
#   3. Inspect highFrac_sample_metadata_list and highFrac_sample_metadata.
#   4. The final section writes highFrac_sample_metadata_list.RDS.
#
# This is cell-type metadata shared by all CTSE methods and top-gene groups.
# It therefore has one row per dataset, config, truth type, and cell type.
#
# Sample definitions:
#   test samples with inferred fraction =
#     test samples intersect config-specific InstaPrism fraction samples
#
#   high-fraction samples =
#     test samples with inferred fraction
#     intersect config-specific InstaPrism fraction > instaprism_min_frac
#
#   all evaluation samples =
#     test samples
#     intersect truth CTSE samples
#     intersect the config-specific truth fraction > min_frac
#
#   filtered evaluation samples =
#     all evaluation samples intersect high-fraction samples
#
# inferred_ct_abundance_all_test_samples matches the abundance basis used for
# the archived adaptive top-N assignment. The explicitly named
# inferred_ct_abundance_after_filtering is the mean inferred fraction among
# retained high-fraction samples.


## settings ----

# Set to NULL to include all annotated benchmark datasets.
dataset_include = NULL

# High-fraction evaluation is currently defined for these independent-reference
# configurations. A dataset/config pair is skipped when its inputs are absent.
config_include = c("config01", "config02")

instaprism_min_frac = 0.1


## repository inputs ----

repo_root_candidates = unique(normalizePath(
  c(getwd(), file.path(getwd(), ".."), file.path(getwd(), "../..")),
  mustWork = FALSE
))

repo_root = repo_root_candidates[file.exists(file.path(
  repo_root_candidates,
  "DALE_Eval",
  "configs",
  "eval_configs.txt"
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

metadata_output_path = file.path(
  repo_root,
  "benchmark_summary",
  "sample_filter_effect_summary",
  "highFrac_sample_metadata_list.RDS"
)

deconv_configs = read.delim(
  file.path(repo_root, "DALE_Eval", "configs", "deconv_configs.txt"),
  sep = "\t",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

eval_configs = read.delim(
  file.path(repo_root, "DALE_Eval", "configs", "eval_configs.txt"),
  sep = "\t",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

ref_assignment = read.delim(
  file.path(
    repo_root,
    "DALE_Eval",
    "configs",
    "benchmark_ref_assignment.txt"
  ),
  sep = "\t",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

cell_type_mapping = read.delim(
  file.path(repo_root, "Indep_scReference", "cell_type_mapping.txt"),
  sep = "\t",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

dataset_info = read_benchmark_dataset_info(repo_root)


## selected datasets and config contracts ----

config_contract = merge(
  deconv_configs,
  eval_configs[, c(
    "config_id",
    "expected_truth_type",
    "filter_frac",
    "min_frac",
    "min_n_sample"
  )],
  by = "config_id",
  all = FALSE,
  sort = FALSE
)

config_contract = config_contract[
  config_contract$refType == "indep",
  ,
  drop = FALSE
]

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

datasets = unique(dataset_info$dataset)

if (!is.null(dataset_include)) {
  missing_datasets = setdiff(dataset_include, datasets)
  if (length(missing_datasets) > 0L) {
    stop(
      "Datasets missing from benchmark_dataset_info.txt: ",
      paste(missing_datasets, collapse = ", ")
    )
  }
  datasets = dataset_include
}


## small helpers ----

empty_metadata_table = function() {
  data.frame(
    dataset = character(),
    tissue = character(),
    cell_type_raw = character(),
    cell_type = character(),
    indep_ref_cell_type = character(),
    config_id = character(),
    truth_type = character(),
    instaprism_min_frac = numeric(),
    truth_min_frac = numeric(),
    min_n_sample = integer(),
    n_test_samples_with_inferred_fraction = integer(),
    n_highFrac_samples = integer(),
    sample_fraction_left = numeric(),
    sample_percent_left = numeric(),
    n_eval_samples_all = integer(),
    n_eval_samples_filtered = integer(),
    eval_sample_fraction_left = numeric(),
    eval_sample_percent_left = numeric(),
    inferred_ct_abundance_all_test_samples = numeric(),
    inferred_ct_abundance_after_filtering = numeric(),
    passes_min_n_sample = logical(),
    stringsAsFactors = FALSE
  )
}

read_numeric_matrix = function(path) {
  x = read.delim(
    path,
    sep = "\t",
    row.names = 1,
    check.names = FALSE
  )
  x = as.matrix(x)
  storage.mode(x) = "numeric"
  x
}

read_test_samples = function(path) {
  sample_split = read.delim(
    path,
    sep = "\t",
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  unique(sample_split$sampleIDs[
    sample_split$group == "test" &
      !is.na(sample_split$sampleIDs) &
      sample_split$sampleIDs != ""
  ])
}

read_truth_files = function(truth_dir) {
  truth_paths = list.files(
    truth_dir,
    pattern = "[.]txt([.]gz)?$",
    full.names = TRUE
  )

  setNames(
    truth_paths,
    sub("[.]txt([.]gz)?$", "", basename(truth_paths))
  )
}

read_truth_sample_names = function(path) {
  truth_stub = read.delim(
    path,
    sep = "\t",
    row.names = 1,
    check.names = FALSE,
    nrows = 1
  )
  colnames(truth_stub)
}

mean_or_na = function(values) {
  values = values[is.finite(values)]
  if (length(values) == 0L) {
    return(NA_real_)
  }
  mean(values)
}

ratio_or_na = function(numerator, denominator) {
  if (denominator == 0L) {
    return(NA_real_)
  }
  numerator / denominator
}

format_cutoff = function(value) {
  format(value, scientific = FALSE, trim = TRUE)
}


## build cell-type metadata ----

highFrac_sample_metadata_list = setNames(
  vector("list", length(datasets)),
  datasets
)

for (dataset_name in datasets) {
  message("Dataset: ", dataset_name)

  dataset_root = file.path(
    repo_root,
    "Benchmarking_obj",
    dataset_name
  )

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
  indep_ref = assignment_row$indep_ref[[1]]

  mapping_one = cell_type_mapping[
    cell_type_mapping$dataset == dataset_name &
      cell_type_mapping$indep_ref == indep_ref,
    c("target_cell_type", "indep_ref_cell_type"),
    drop = FALSE
  ]
  mapping_one = unique(mapping_one)
  if (anyDuplicated(mapping_one$target_cell_type)) {
    stop(dataset_name, ": target cell-type mapping is not unique")
  }

  annotation_one = dataset_info[
    dataset_info$dataset == dataset_name,
    c("cell_type_raw", "cell_type", "tissue"),
    drop = FALSE
  ]

  sample_split_path = file.path(
    dataset_root,
    "self_reference",
    "sample_split.txt"
  )
  if (!file.exists(sample_split_path)) {
    message("  skip dataset: missing sample_split.txt")
    highFrac_sample_metadata_list[[dataset_name]] = empty_metadata_table()
    next
  }
  test_samples = read_test_samples(sample_split_path)

  dataset_rows = list()

  for (config_index in seq_len(nrow(config_contract))) {
    config_row = config_contract[config_index, , drop = FALSE]
    config_id = config_row$config_id[[1]]
    truth_type = config_row$expected_truth_type[[1]]
    truth_min_frac = config_row$min_frac[[1]]
    min_n_sample = as.integer(config_row$min_n_sample[[1]])

    truth_dir = file.path(
      dataset_root,
      "ctse_truth",
      truth_type
    )
    truth_fraction_path = file.path(
      dataset_root,
      "frac_input",
      paste0(config_row$filter_frac[[1]], ".txt")
    )
    sample_mask_path = file.path(
      dataset_root,
      "evaluation_metadata",
      "sample_mask",
      paste0(
        "InstaPrismfrac_",
        config_id,
        "_gt",
        format_cutoff(instaprism_min_frac),
        ".txt"
      )
    )
    instaprism_fraction_path = file.path(
      dataset_root,
      "deconv_res",
      config_slug(config_row),
      "InstaPrism",
      "InstaPrismfrac.txt"
    )

    required_paths = c(
      truth_dir,
      truth_fraction_path,
      sample_mask_path,
      instaprism_fraction_path
    )
    if (!all(file.exists(required_paths))) {
      message("  ", config_id, " skip: incomplete high-fraction inputs")
      next
    }

    truth_files = read_truth_files(truth_dir)
    truth_fraction = read_numeric_matrix(truth_fraction_path)
    sample_mask = read_numeric_matrix(sample_mask_path)
    instaprism_fraction = read_numeric_matrix(instaprism_fraction_path)

    mapping_config = mapping_one[
      mapping_one$target_cell_type %in% names(truth_files) &
        mapping_one$target_cell_type %in% colnames(truth_fraction) &
        mapping_one$target_cell_type %in% colnames(sample_mask),
      ,
      drop = FALSE
    ]

    source_column = ifelse(
      mapping_config$indep_ref_cell_type %in%
        colnames(instaprism_fraction),
      mapping_config$indep_ref_cell_type,
      mapping_config$target_cell_type
    )
    mapping_config$instaprism_column = source_column
    mapping_config = mapping_config[
      mapping_config$instaprism_column %in%
        colnames(instaprism_fraction),
      ,
      drop = FALSE
    ]

    for (cell_type_index in seq_len(nrow(mapping_config))) {
      mapping_row = mapping_config[
        cell_type_index,
        ,
        drop = FALSE
      ]
      cell_type_raw = mapping_row$target_cell_type[[1]]
      instaprism_column = mapping_row$instaprism_column[[1]]

      annotation_row = annotation_one[
        annotation_one$cell_type_raw == cell_type_raw,
        ,
        drop = FALSE
      ]
      if (nrow(annotation_row) != 1L) {
        stop(
          dataset_name,
          " / ", cell_type_raw,
          ": expected one benchmark dataset annotation row"
        )
      }

      test_samples_with_inferred_fraction = Reduce(
        intersect,
        list(
          test_samples,
          rownames(instaprism_fraction),
          rownames(sample_mask)
        )
      )

      mask_values = sample_mask[, cell_type_raw]
      highFrac_samples = rownames(sample_mask)[
        is.finite(mask_values) & mask_values == 1
      ]
      highFrac_samples = intersect(
        test_samples_with_inferred_fraction,
        highFrac_samples
      )

      truth_samples = read_truth_sample_names(
        truth_files[[cell_type_raw]]
      )
      truth_fraction_values = truth_fraction[, cell_type_raw]
      truth_fraction_samples = rownames(truth_fraction)[
        is.finite(truth_fraction_values) &
          truth_fraction_values > truth_min_frac
      ]

      all_eval_samples = Reduce(
        intersect,
        list(
          test_samples,
          truth_samples,
          truth_fraction_samples
        )
      )
      filtered_eval_samples = intersect(
        all_eval_samples,
        highFrac_samples
      )

      n_test_samples_with_inferred_fraction = length(
        test_samples_with_inferred_fraction
      )
      n_highFrac_samples = length(highFrac_samples)
      n_eval_samples_all = length(all_eval_samples)
      n_eval_samples_filtered = length(filtered_eval_samples)

      sample_fraction_left = ratio_or_na(
        n_highFrac_samples,
        n_test_samples_with_inferred_fraction
      )
      eval_sample_fraction_left = ratio_or_na(
        n_eval_samples_filtered,
        n_eval_samples_all
      )

      inferred_ct_abundance_all_test_samples = mean_or_na(
        instaprism_fraction[
          test_samples_with_inferred_fraction,
          instaprism_column
        ]
      )
      inferred_ct_abundance_after_filtering = mean_or_na(
        instaprism_fraction[
          highFrac_samples,
          instaprism_column
        ]
      )

      dataset_rows[[length(dataset_rows) + 1L]] = data.frame(
        dataset = dataset_name,
        tissue = annotation_row$tissue[[1]],
        cell_type_raw = cell_type_raw,
        cell_type = annotation_row$cell_type[[1]],
        indep_ref_cell_type = mapping_row$indep_ref_cell_type[[1]],
        config_id = config_id,
        truth_type = truth_type,
        instaprism_min_frac = instaprism_min_frac,
        truth_min_frac = truth_min_frac,
        min_n_sample = min_n_sample,
        n_test_samples_with_inferred_fraction =
          n_test_samples_with_inferred_fraction,
        n_highFrac_samples = n_highFrac_samples,
        sample_fraction_left = sample_fraction_left,
        sample_percent_left = 100 * sample_fraction_left,
        n_eval_samples_all = n_eval_samples_all,
        n_eval_samples_filtered = n_eval_samples_filtered,
        eval_sample_fraction_left = eval_sample_fraction_left,
        eval_sample_percent_left = 100 * eval_sample_fraction_left,
        inferred_ct_abundance_all_test_samples =
          inferred_ct_abundance_all_test_samples,
        inferred_ct_abundance_after_filtering =
          inferred_ct_abundance_after_filtering,
        passes_min_n_sample = n_eval_samples_filtered >= min_n_sample,
        stringsAsFactors = FALSE
      )
    }
  }

  highFrac_sample_metadata_list[[dataset_name]] = if (
    length(dataset_rows) > 0L
  ) {
    result = do.call(rbind, dataset_rows)
    rownames(result) = NULL
    result
  } else {
    empty_metadata_table()
  }

  message(
    "  Metadata rows: ",
    nrow(highFrac_sample_metadata_list[[dataset_name]])
  )
}

highFrac_sample_metadata = do.call(
  rbind,
  highFrac_sample_metadata_list
)
rownames(highFrac_sample_metadata) = NULL


## export ----

saveRDS(
  highFrac_sample_metadata_list,
  file = metadata_output_path
)

message("Saved: ", metadata_output_path)
