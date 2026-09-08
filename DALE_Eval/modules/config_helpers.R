find_repo_root <- function(start = getwd()) {
  cur <- normalizePath(start, mustWork = TRUE)
  repeat {
    if (
      file.exists(file.path(cur, "DALE_Eval", "configs", "deconv_configs.txt")) &&
        dir.exists(file.path(cur, "Benchmarking_obj")) &&
        dir.exists(file.path(cur, "DALE_Eval"))
    ) {
      return(cur)
    }
    parent <- dirname(cur)
    if (identical(parent, cur)) {
      stop("Could not find repo root from: ", start)
    }
    cur <- parent
  }
}

script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[[1]]))))
  }
  getwd()
}

normalize_dir <- function(path) {
  sub("/+$", "", path)
}

normalize_config_value <- function(x) {
  x <- as.character(x[[1]])
  trimws(tolower(x))
}

config_bulk_normalization <- function(config_row) {
  normalize_config_value(config_row$bulk_normalization)
}

config_slug <- function(config_row) {
  bulk_normalization <- config_bulk_normalization(config_row)
  normalization_suffix <- if (bulk_normalization == "cpm") {
    ""
  } else {
    paste0("__norm-", bulk_normalization)
  }

  paste0(
    config_row$config_id,
    "_bulk-", config_row$bulk_input,
    "__frac-", config_row$frac_input,
    "__ref-", config_row$refType,
    normalization_suffix
  )
}

read_deconv_config <- function(config_id, repo_root = find_repo_root()) {
  config_path <- file.path(repo_root, "DALE_Eval", "configs", "deconv_configs.txt")
  configs <- read.delim(config_path, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
  required_cols <- c("config_id", "bulk_input", "bulk_scale", "bulk_normalization", "frac_input", "refType")
  missing_cols <- setdiff(required_cols, colnames(configs))
  if (length(missing_cols) > 0) {
    stop("deconv_configs.txt missing columns: ", paste(missing_cols, collapse = ", "))
  }
  config_row <- configs[configs$config_id == config_id, , drop = FALSE]
  if (nrow(config_row) != 1) {
    stop("Expected exactly one deconv config for config_id=", config_id, ", found ", nrow(config_row))
  }
  config_row
}

find_deconv_config <- function(
  bulk_input,
  frac_input,
  refType,
  bulk_normalization = "cpm",
  repo_root = find_repo_root()
) {
  config_path <- file.path(repo_root, "DALE_Eval", "configs", "deconv_configs.txt")
  configs <- read.delim(config_path, sep = "	", stringsAsFactors = FALSE, check.names = FALSE)
  required_cols <- c("config_id", "bulk_input", "bulk_scale", "bulk_normalization", "frac_input", "refType")
  missing_cols <- setdiff(required_cols, colnames(configs))
  if (length(missing_cols) > 0) {
    stop("deconv_configs.txt missing columns: ", paste(missing_cols, collapse = ", "))
  }

  config_row <- configs[
    configs$bulk_input == bulk_input &
      tolower(configs$bulk_normalization) == tolower(bulk_normalization) &
      configs$frac_input == frac_input &
      configs$refType == refType,
    ,
    drop = FALSE
  ]
  if (nrow(config_row) != 1) {
    stop(
      "Expected exactly one deconv config for bulk_input=", bulk_input,
      ", bulk_normalization=", bulk_normalization,
      ", frac_input=", frac_input,
      ", refType=", refType,
      "; found ", nrow(config_row)
    )
  }
  config_row
}

resolve_indep_ref_dir <- function(dataset, repo_root = find_repo_root()) {
  assignment_path <- file.path(repo_root, "DALE_Eval", "configs", "benchmark_ref_assignment.txt")
  assignments <- read.delim(assignment_path, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
  ref_row <- assignments[assignments$dataset == dataset, , drop = FALSE]
  if (nrow(ref_row) != 1) {
    stop("Expected exactly one independent reference assignment for dataset=", dataset)
  }
  file.path(repo_root, "Indep_scReference", ref_row$indep_ref[[1]])
}

resolve_deconv_paths <- function(dataset, config_id, repo_root = find_repo_root()) {
  config_row <- read_deconv_config(config_id, repo_root = repo_root)
  ref_type <- config_row$refType[[1]]
  if (!ref_type %in% c("indep", "self")) {
    stop("Invalid refType in deconv config ", config_id, ": ", ref_type)
  }

  obj_dir <- file.path(repo_root, "Benchmarking_obj", dataset)
  bulk_path <- file.path(obj_dir, "bulk_input", paste0(config_row$bulk_input[[1]], ".txt"))
  frac_path <- file.path(obj_dir, "frac_input", paste0(config_row$frac_input[[1]], ".txt"))
  ref_dir <- if (ref_type == "self") {
    file.path(obj_dir, "self_reference")
  } else {
    resolve_indep_ref_dir(dataset, repo_root = repo_root)
  }
  output_dir <- file.path(obj_dir, "deconv_res", config_slug(config_row))
  runtime_path <- file.path(obj_dir, "logs", "deconv_runs.txt")

  list(
    dataset = dataset,
    config_id = config_id,
    config = config_row,
    refType = ref_type,
    obj_dir = obj_dir,
    ref_dir = normalize_dir(ref_dir),
    bulk_path = bulk_path,
    frac_path = frac_path,
    output_dir = output_dir,
    runtime_path = runtime_path
  )
}

read_scitd_config <- function(
  scitd_config_id = "tutorial_v1",
  repo_root = find_repo_root()
) {
  config_path <- file.path(
    repo_root,
    "DALE_Eval",
    "configs",
    "scITD_configs.txt"
  )

  if (!file.exists(config_path)) {
    stop("Missing scITD config file: ", config_path)
  }

  configs <- read.delim(
    config_path,
    sep = "\t",
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  required_cols <- c(
    "scitd_config_id",
    "deconv_config_ids",
    "truth_type",
    "use_test_samples",
    "truth_gene_scope",
    "method_gene_scope",
    "method_sample_scope",
    "drop_zero_profiles",
    "zero_profile_tolerance",
    "truth_donor_min_cells",
    "truth_norm_method",
    "scale_factor",
    "vargenes_method",
    "vargenes_thresh",
    "scale_var",
    "var_scale_power",
    "donor_rank",
    "gene_rank",
    "tucker_type",
    "rotation_type",
    "inferred_negative_action",
    "n_cores",
    "random_seed"
  )

  missing_cols <- setdiff(required_cols, colnames(configs))
  if (length(missing_cols) > 0) {
    stop(
      "scITD_configs.txt missing columns: ",
      paste(missing_cols, collapse = ", ")
    )
  }

  config_row <- configs[
    configs$scitd_config_id == scitd_config_id,
    ,
    drop = FALSE
  ]

  if (nrow(config_row) != 1) {
    stop(
      "Expected exactly one scITD config for scitd_config_id=",
      scitd_config_id,
      ", found ",
      nrow(config_row)
    )
  }

  scalar_character <- function(field) {
    value <- trimws(as.character(config_row[[field]][[1]]))
    if (!nzchar(value)) {
      stop("Empty value for scITD config field: ", field)
    }
    value
  }

  scalar_boolean <- function(field) {
    value <- tolower(scalar_character(field))
    if (!value %in% c("true", "false")) {
      stop("scITD config field ", field, " must be true or false")
    }
    identical(value, "true")
  }

  scalar_numeric <- function(field) {
    raw_value <- scalar_character(field)
    value <- suppressWarnings(as.numeric(raw_value))
    if (!is.finite(value)) {
      stop("scITD config field ", field, " must be finite numeric")
    }
    value
  }

  scalar_integer <- function(field) {
    value <- scalar_numeric(field)
    if (value != floor(value)) {
      stop("scITD config field ", field, " must be an integer")
    }
    as.integer(value)
  }

  deconv_config_ids <- trimws(strsplit(
    scalar_character("deconv_config_ids"),
    ",",
    fixed = TRUE
  )[[1]])
  deconv_config_ids <- deconv_config_ids[nzchar(deconv_config_ids)]
  if (length(deconv_config_ids) == 0 || anyDuplicated(deconv_config_ids)) {
    stop("deconv_config_ids must contain unique, comma-separated IDs")
  }
  if (any(!deconv_config_ids %in% c("config01", "config02"))) {
    stop("The scITD benchmark is restricted to config01 and config02")
  }

  config <- list(
    scitd_config_id = scalar_character("scitd_config_id"),
    deconv_config_ids = deconv_config_ids,
    truth_type = scalar_character("truth_type"),
    use_test_samples = scalar_boolean("use_test_samples"),
    truth_gene_scope = scalar_character("truth_gene_scope"),
    method_gene_scope = scalar_character("method_gene_scope"),
    method_sample_scope = scalar_character("method_sample_scope"),
    drop_zero_profiles = scalar_boolean("drop_zero_profiles"),
    zero_profile_tolerance = scalar_numeric("zero_profile_tolerance"),
    truth_donor_min_cells = scalar_integer("truth_donor_min_cells"),
    truth_norm_method = scalar_character("truth_norm_method"),
    scale_factor = scalar_numeric("scale_factor"),
    vargenes_method = scalar_character("vargenes_method"),
    vargenes_thresh = scalar_numeric("vargenes_thresh"),
    scale_var = scalar_boolean("scale_var"),
    var_scale_power = scalar_numeric("var_scale_power"),
    donor_rank = scalar_integer("donor_rank"),
    gene_rank = scalar_integer("gene_rank"),
    tucker_type = scalar_character("tucker_type"),
    rotation_type = scalar_character("rotation_type"),
    inferred_negative_action = scalar_character("inferred_negative_action"),
    n_cores = scalar_integer("n_cores"),
    random_seed = scalar_integer("random_seed")
  )

  supported_policy <- c(
    truth_type = "sumcount",
    truth_gene_scope = "indep_ref_overlap",
    method_gene_scope = "truth_selected_intersection",
    method_sample_scope = "available_nonzero_by_result",
    inferred_negative_action = "gene_celltype_shift"
  )
  for (field in names(supported_policy)) {
    if (!identical(config[[field]], unname(supported_policy[[field]]))) {
      stop(
        "Unsupported ", field, " in scITD config: ", config[[field]],
        ". Expected ", supported_policy[[field]], "."
      )
    }
  }

  if (!config$use_test_samples) {
    stop("The current scITD benchmark requires use_test_samples=true")
  }
  if (!config$drop_zero_profiles) {
    stop("The current scITD benchmark requires drop_zero_profiles=true")
  }
  if (config$zero_profile_tolerance < 0) {
    stop("zero_profile_tolerance must be nonnegative")
  }
  if (config$truth_donor_min_cells != 0L) {
    stop(
      "truth_donor_min_cells must be 0 because each truth pseudobulk is ",
      "represented by one pseudo-cell"
    )
  }
  if (!config$truth_norm_method %in% c("trim", "regular")) {
    stop("truth_norm_method must be trim or regular")
  }
  if (!config$vargenes_method %in% c("norm_var", "norm_var_pvals", "anova")) {
    stop("Unsupported vargenes_method: ", config$vargenes_method)
  }
  if (config$scale_factor <= 0 || config$vargenes_thresh <= 0) {
    stop("scale_factor and vargenes_thresh must be positive")
  }
  if (config$var_scale_power < 0) {
    stop("var_scale_power must be nonnegative")
  }
  if (config$donor_rank < 1L || config$gene_rank < 1L) {
    stop("donor_rank and gene_rank must be positive integers")
  }
  if (!identical(config$tucker_type, "regular")) {
    stop("tucker_type must be regular")
  }
  if (!config$rotation_type %in% c("hybrid", "ica_dsc", "ica_lds")) {
    stop("Unsupported rotation_type: ", config$rotation_type)
  }
  if (config$n_cores < 1L) {
    stop("n_cores must be a positive integer")
  }
  if (config$random_seed < 0L) {
    stop("random_seed must be a nonnegative integer")
  }

  for (config_id in config$deconv_config_ids) {
    deconv_config <- read_deconv_config(config_id, repo_root = repo_root)
    if (!identical(as.character(deconv_config$refType[[1]]), "indep")) {
      stop("scITD deconvolution configs must use refType=indep: ", config_id)
    }
  }

  structure(config, class = c("scitd_benchmark_config", "list"))
}

read_scitd_method_preprocessing <- function(repo_root = find_repo_root()) {
  config_path <- file.path(
    repo_root,
    "DALE_Eval",
    "configs",
    "scITD_method_preprocessing.txt"
  )

  if (!file.exists(config_path)) {
    stop("Missing scITD method preprocessing config: ", config_path)
  }

  preprocessing <- read.delim(
    config_path,
    sep = "\t",
    stringsAsFactors = FALSE,
    check.names = FALSE,
    colClasses = "character"
  )

  required_cols <- c("method", "preprocess_mode")
  missing_cols <- setdiff(required_cols, colnames(preprocessing))
  unexpected_cols <- setdiff(colnames(preprocessing), required_cols)
  if (length(missing_cols) > 0) {
    stop(
      "scITD_method_preprocessing.txt missing columns: ",
      paste(missing_cols, collapse = ", ")
    )
  }
  if (length(unexpected_cols) > 0) {
    stop(
      "scITD_method_preprocessing.txt has unexpected columns: ",
      paste(unexpected_cols, collapse = ", ")
    )
  }
  if (nrow(preprocessing) == 0) {
    stop("scITD_method_preprocessing.txt contains no method mappings")
  }

  preprocessing <- preprocessing[, required_cols, drop = FALSE]
  preprocessing$method <- trimws(preprocessing$method)
  preprocessing$preprocess_mode <- tolower(trimws(
    preprocessing$preprocess_mode
  ))

  if (any(is.na(preprocessing$method) | !nzchar(preprocessing$method))) {
    stop("scITD method preprocessing contains an empty method name")
  }
  if (
    any(
      is.na(preprocessing$preprocess_mode) |
        !nzchar(preprocessing$preprocess_mode)
    )
  ) {
    stop("scITD method preprocessing contains an empty preprocess_mode")
  }
  if (anyDuplicated(preprocessing$method)) {
    duplicated_methods <- unique(preprocessing$method[
      duplicated(preprocessing$method)
    ])
    stop(
      "Duplicate scITD method preprocessing mappings: ",
      paste(duplicated_methods, collapse = ", ")
    )
  }

  supported_modes <- c("log1p", "already_transformed")
  unsupported_modes <- setdiff(
    unique(preprocessing$preprocess_mode),
    supported_modes
  )
  if (length(unsupported_modes) > 0) {
    stop(
      "Unsupported scITD preprocess_mode value(s): ",
      paste(unsupported_modes, collapse = ", "),
      ". Supported values: ",
      paste(supported_modes, collapse = ", ")
    )
  }

  rownames(preprocessing) <- NULL
  structure(
    preprocessing,
    class = c("scitd_method_preprocessing", "data.frame")
  )
}
