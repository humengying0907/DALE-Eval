# Shared scITD benchmark helpers.
#
# Truth sum-count pseudobulks use scITD's make_new_container() and
# form_tensor() workflow. Continuous inferred CTSE estimates bypass
# make_new_container(), enter scITD at the pseudobulk stage, and receive
# method-scale-appropriate preprocessing.

scitd_truth_result_dir <- function(dataset, repo_root = find_repo_root()) {
  file.path(
    repo_root,
    "Benchmarking_obj",
    dataset,
    "scITD_res",
    "truth_sumcount"
  )
}

scitd_method_result_dir <- function(
  dataset,
  config_id,
  method,
  lib_norm = FALSE,
  repo_root = find_repo_root()
) {
  if (
    length(lib_norm) != 1L ||
      !is.logical(lib_norm) ||
      is.na(lib_norm)
  ) {
    stop("lib_norm must be one non-missing logical value")
  }

  config_row <- read_deconv_config(config_id, repo_root = repo_root)
  method_output_label <- paste0(
    method,
    if (lib_norm) "_libnorm" else ""
  )
  file.path(
    repo_root,
    "Benchmarking_obj",
    dataset,
    "scITD_res",
    config_slug(config_row),
    method_output_label
  )
}

scitd_cell_type_map <- function(truth_cell_types) {
  truth_cell_types <- as.character(truth_cell_types)
  if (
    length(truth_cell_types) == 0 ||
      anyNA(truth_cell_types) ||
      any(!nzchar(truth_cell_types)) ||
      anyDuplicated(truth_cell_types)
  ) {
    stop("Truth cell-type names must be nonempty and unique")
  }

  internal_cell_types <- gsub("[._:]", "", truth_cell_types)
  if (any(!nzchar(internal_cell_types))) {
    stop("A truth cell-type name is empty after scITD name sanitization")
  }
  if (anyDuplicated(internal_cell_types)) {
    collisions <- unique(internal_cell_types[duplicated(internal_cell_types)])
    stop(
      "Truth cell-type names collide after removing '.', '_', and ':': ",
      paste(collisions, collapse = ", ")
    )
  }

  data.frame(
    truth_cell_type = truth_cell_types,
    scitd_cell_type = internal_cell_types,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

read_scitd_expression_matrix <- function(
  path,
  require_nonnegative = FALSE,
  require_integer = FALSE
) {
  if (!file.exists(path)) {
    stop("Missing expression matrix: ", path)
  }

  connection <- if (grepl("\\.gz$", path, ignore.case = TRUE)) {
    gzfile(path, open = "rt")
  } else {
    file(path, open = "rt")
  }
  on.exit(try(close(connection), silent = TRUE), add = TRUE)

  values <- read.delim(
    connection,
    header = TRUE,
    row.names = 1,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  if (nrow(values) == 0 || ncol(values) == 0) {
    stop("Expression matrix has zero rows or columns: ", path)
  }

  values <- as.matrix(values)
  suppressWarnings(storage.mode(values) <- "double")

  if (
    is.null(rownames(values)) ||
      anyNA(rownames(values)) ||
      any(!nzchar(rownames(values))) ||
      anyDuplicated(rownames(values))
  ) {
    stop("Gene names must be nonempty and unique in: ", path)
  }
  if (
    is.null(colnames(values)) ||
      anyNA(colnames(values)) ||
      any(!nzchar(colnames(values))) ||
      anyDuplicated(colnames(values))
  ) {
    stop("Sample names must be nonempty and unique in: ", path)
  }
  if (any(!is.finite(values))) {
    stop("Non-finite or nonnumeric values found in: ", path)
  }
  if (require_nonnegative && any(values < 0)) {
    stop("Negative values found in nonnegative matrix: ", path)
  }
  if (
    require_integer &&
      any(abs(values - round(values)) > sqrt(.Machine$double.eps))
  ) {
    stop("Non-integer values found in truth sum-count matrix: ", path)
  }

  values
}

scitd_truth_files <- function(
  dataset,
  truth_type = "sumcount",
  repo_root = find_repo_root()
) {
  truth_dir <- file.path(
    repo_root,
    "Benchmarking_obj",
    dataset,
    "ctse_truth",
    truth_type
  )
  if (!dir.exists(truth_dir)) {
    stop("Missing truth CTSE directory: ", truth_dir)
  }

  paths <- list.files(
    truth_dir,
    pattern = "\\.txt\\.gz$",
    full.names = TRUE
  )
  if (length(paths) == 0) {
    stop("No truth CTSE .txt.gz files found in: ", truth_dir)
  }

  truth_cell_types <- sub("\\.txt\\.gz$", "", basename(paths))
  if (anyDuplicated(truth_cell_types)) {
    stop("Duplicated truth cell-type filenames in: ", truth_dir)
  }

  stats::setNames(paths, truth_cell_types)
}

read_scitd_test_samples <- function(dataset, repo_root = find_repo_root()) {
  split_path <- file.path(
    repo_root,
    "Benchmarking_obj",
    dataset,
    "self_reference",
    "sample_split.txt"
  )
  if (!file.exists(split_path)) {
    stop("Missing sample split file: ", split_path)
  }

  sample_split <- read.delim(
    split_path,
    header = TRUE,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  required_cols <- c("group", "sampleIDs")
  missing_cols <- setdiff(required_cols, colnames(sample_split))
  if (length(missing_cols) > 0) {
    stop(
      "sample_split.txt missing columns: ",
      paste(missing_cols, collapse = ", ")
    )
  }

  group <- tolower(trimws(as.character(sample_split$group)))
  samples <- as.character(sample_split$sampleIDs[group == "test"])
  samples <- samples[!is.na(samples) & nzchar(samples)]
  samples <- unique(samples)
  if (length(samples) == 0) {
    stop("No test samples found in: ", split_path)
  }

  samples
}

scitd_reference_genes <- function(dataset, repo_root = find_repo_root()) {
  reference_dir <- resolve_indep_ref_dir(dataset, repo_root = repo_root)
  reference_path <- file.path(reference_dir, "refPhi.RDS")
  if (!file.exists(reference_path)) {
    stop("Missing independent-reference refPhi.RDS: ", reference_path)
  }

  reference_object <- readRDS(reference_path)

  # S4 slots are stored as attributes in the serialized object. Reading this
  # first avoids requiring the reference object's defining package merely to
  # recover its gene names.
  phi_cs <- attr(reference_object, "phi.cs", exact = TRUE)
  if (is.null(phi_cs) && is.environment(reference_object) && exists(
    "phi.cs",
    envir = reference_object,
    inherits = FALSE
  )) {
    phi_cs <- get("phi.cs", envir = reference_object, inherits = FALSE)
  } else if (
    is.null(phi_cs) &&
      is.list(reference_object) &&
      !is.null(reference_object$phi.cs)
  ) {
    phi_cs <- reference_object$phi.cs
  }

  reference_genes <- rownames(phi_cs)
  if (
    is.null(reference_genes) ||
      length(reference_genes) == 0 ||
      anyNA(reference_genes) ||
      any(!nzchar(reference_genes))
  ) {
    stop("Could not recover genes from the phi.cs component of: ", reference_path)
  }

  unique(reference_genes)
}

scitd_ordered_common_values <- function(values_by_source) {
  if (length(values_by_source) == 0) {
    return(character())
  }

  common_values <- as.character(values_by_source[[1]])
  if (length(values_by_source) > 1) {
    for (i in 2:length(values_by_source)) {
      common_values <- common_values[
        common_values %in% as.character(values_by_source[[i]])
      ]
    }
  }
  common_values
}

scitd_zero_profiles <- function(
  matrices,
  sample_order,
  tolerance = 1e-12
) {
  rows <- lapply(names(matrices), function(cell_type) {
    matrix_values <- matrices[[cell_type]][
      ,
      sample_order,
      drop = FALSE
    ]
    profile_sizes <- colSums(abs(matrix_values))
    zero_samples <- sample_order[
      !is.finite(profile_sizes) | profile_sizes <= tolerance
    ]
    data.frame(
      sample = zero_samples,
      cell_type = rep(cell_type, length(zero_samples)),
      stringsAsFactors = FALSE
    )
  })

  profiles <- do.call(rbind, rows)
  rownames(profiles) <- NULL
  zero_samples <- sample_order[sample_order %in% unique(profiles$sample)]

  list(samples = zero_samples, profiles = profiles)
}

validate_scitd_ranks <- function(n_samples, n_genes, config, label) {
  if (n_samples <= config$donor_rank) {
    stop(
      label, " has ", n_samples, " usable samples; donor_rank=",
      config$donor_rank, " requires more samples"
    )
  }
  if (n_genes < config$gene_rank) {
    stop(
      label, " has ", n_genes, " usable genes; gene_rank=",
      config$gene_rank, " cannot be fitted"
    )
  }
  invisible(TRUE)
}

validate_scitd_tensor <- function(container, samples, genes, cell_types) {
  tensor_data <- container$tensor_data
  if (is.null(tensor_data) || length(tensor_data) < 4) {
    stop("scITD did not create tensor_data")
  }
  if (!identical(tensor_data[[1]], samples)) {
    stop("scITD tensor sample order differs from the expected order")
  }
  if (!identical(tensor_data[[2]], genes)) {
    stop("scITD tensor gene order differs from the expected order")
  }
  if (!identical(tensor_data[[3]], cell_types)) {
    stop("scITD tensor cell-type order differs from the expected order")
  }
  if (any(!is.finite(tensor_data[[4]]))) {
    stop("scITD tensor contains non-finite values")
  }
  invisible(TRUE)
}

validate_scitd_cell_type_variation <- function(
  container,
  cell_type_labels = NULL,
  tolerance = 1e-12
) {
  if (!is.finite(tolerance) || tolerance < 0) {
    stop("Cell-type variation tolerance must be finite and nonnegative")
  }

  tensor_data <- container$tensor_data
  cell_types <- tensor_data[[3]]
  tensor <- tensor_data[[4]]
  if (is.null(cell_type_labels)) {
    cell_type_labels <- stats::setNames(cell_types, cell_types)
  }
  missing_labels <- setdiff(cell_types, names(cell_type_labels))
  if (length(missing_labels) > 0) {
    stop(
      "Missing display labels for scITD cell types: ",
      paste(missing_labels, collapse = ", ")
    )
  }

  for (cell_type_index in seq_along(cell_types)) {
    cell_type <- cell_types[[cell_type_index]]
    cell_type_slice <- tensor[, , cell_type_index, drop = FALSE]
    dim(cell_type_slice) <- dim(tensor)[1:2]
    gene_standard_deviations <- apply(
      cell_type_slice,
      MARGIN = 2,
      FUN = stats::sd
    )
    has_patient_variation <- any(
      is.finite(gene_standard_deviations) &
        gene_standard_deviations > tolerance
    )
    if (!has_patient_variation) {
      stop(
        unname(cell_type_labels[[cell_type]]),
        " has no patient-level variation after variance QC"
      )
    }
  }
  invisible(TRUE)
}

run_scitd_truth <- function(
  dataset,
  scitd_config,
  repo_root = find_repo_root()
) {
  if (!requireNamespace("scITD", quietly = TRUE)) {
    stop("The scITD package is not available on .libPaths()")
  }
  if (!requireNamespace("Matrix", quietly = TRUE)) {
    stop("The Matrix package is required")
  }
  if (
    identical(scitd_config$truth_norm_method, "trim") &&
      !requireNamespace("edgeR", quietly = TRUE)
  ) {
    stop("edgeR is required for scITD trim normalization of truth counts")
  }

  message("=== scITD truth: ", dataset, " ===")

  truth_paths <- scitd_truth_files(
    dataset,
    truth_type = scitd_config$truth_type,
    repo_root = repo_root
  )
  truth_cell_types <- names(truth_paths)
  cell_type_map <- scitd_cell_type_map(truth_cell_types)

  truth_counts <- lapply(
    truth_paths,
    read_scitd_expression_matrix,
    require_nonnegative = TRUE,
    require_integer = TRUE
  )

  test_samples <- read_scitd_test_samples(dataset, repo_root = repo_root)
  missing_test_samples <- lapply(truth_counts, function(values) {
    setdiff(test_samples, colnames(values))
  })
  missing_counts <- lengths(missing_test_samples)
  if (any(missing_counts > 0)) {
    problem_cell_types <- names(missing_counts)[missing_counts > 0]
    details <- paste0(
      problem_cell_types,
      " (",
      missing_counts[problem_cell_types],
      ")"
    )
    stop(
      "Test samples are missing from truth cell types: ",
      paste(details, collapse = ", ")
    )
  }

  shared_truth_genes <- scitd_ordered_common_values(
    lapply(truth_counts, rownames)
  )
  reference_genes <- scitd_reference_genes(dataset, repo_root = repo_root)
  overlapping_genes <- shared_truth_genes[
    shared_truth_genes %in% reference_genes
  ]
  if (length(overlapping_genes) == 0) {
    stop("No genes overlap truth CTSE and the independent reference")
  }

  truth_counts <- lapply(truth_counts, function(values) {
    values[overlapping_genes, test_samples, drop = FALSE]
  })

  zero_profiles <- scitd_zero_profiles(
    truth_counts,
    sample_order = test_samples,
    tolerance = scitd_config$zero_profile_tolerance
  )
  samples_used <- test_samples[
    !test_samples %in% zero_profiles$samples
  ]
  message(
    "Truth samples: ", length(test_samples),
    "; excluded zero profiles: ", length(zero_profiles$samples),
    "; retained: ", length(samples_used)
  )

  validate_scitd_ranks(
    n_samples = length(samples_used),
    n_genes = length(overlapping_genes),
    config = scitd_config,
    label = paste(dataset, "truth")
  )

  pseudo_count_blocks <- vector("list", length(truth_cell_types))
  pseudo_metadata_blocks <- vector("list", length(truth_cell_types))
  for (i in seq_along(truth_cell_types)) {
    truth_cell_type <- truth_cell_types[[i]]
    internal_cell_type <- cell_type_map$scitd_cell_type[[i]]
    count_block <- truth_counts[[truth_cell_type]][
      overlapping_genes,
      samples_used,
      drop = FALSE
    ]
    pseudo_cell_ids <- paste0(
      "pseudocell__",
      i,
      "__",
      seq_along(samples_used)
    )
    colnames(count_block) <- pseudo_cell_ids
    pseudo_count_blocks[[i]] <- count_block
    pseudo_metadata_blocks[[i]] <- data.frame(
      donors = samples_used,
      ctypes = rep(internal_cell_type, length(samples_used)),
      row.names = pseudo_cell_ids,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  }

  pseudo_count_matrix <- do.call(cbind, pseudo_count_blocks)
  pseudo_metadata <- do.call(rbind, pseudo_metadata_blocks)
  pseudo_metadata$ctypes <- factor(
    pseudo_metadata$ctypes,
    levels = cell_type_map$scitd_cell_type
  )
  if (!identical(colnames(pseudo_count_matrix), rownames(pseudo_metadata))) {
    stop("Truth pseudo-count columns do not match pseudo-metadata rows")
  }

  truth_params <- scITD::initialize_params(
    ctypes_use = cell_type_map$scitd_cell_type,
    ncores = scitd_config$n_cores,
    rand_seed = scitd_config$random_seed
  )
  truth_container <- scITD::make_new_container(
    params = truth_params,
    count_data = Matrix::Matrix(pseudo_count_matrix, sparse = TRUE),
    meta_data = pseudo_metadata
  )
  truth_container <- scITD::form_tensor(
    container = truth_container,
    donor_min_cells = scitd_config$truth_donor_min_cells,
    norm_method = scitd_config$truth_norm_method,
    scale_factor = scitd_config$scale_factor,
    vargenes_method = scitd_config$vargenes_method,
    vargenes_thresh = scitd_config$vargenes_thresh,
    batch_var = NULL,
    scale_var = scitd_config$scale_var,
    var_scale_power = scitd_config$var_scale_power,
    custom_genes = NULL,
    verbose = TRUE
  )

  genes_used <- truth_container$tensor_data[[2]]
  samples_in_tensor <- truth_container$tensor_data[[1]]
  if (length(genes_used) == 0) {
    stop("scITD selected no variable genes for truth")
  }
  if (!setequal(samples_in_tensor, samples_used)) {
    stop("scITD truth tensor did not retain exactly the expected samples")
  }
  if (!identical(genes_used, truth_container$all_vargenes)) {
    stop("scITD truth tensor genes differ from all_vargenes")
  }
  validate_scitd_ranks(
    n_samples = length(samples_in_tensor),
    n_genes = length(genes_used),
    config = scitd_config,
    label = paste(dataset, "truth tensor")
  )
  validate_scitd_tensor(
    truth_container,
    samples = samples_in_tensor,
    genes = genes_used,
    cell_types = cell_type_map$scitd_cell_type
  )

  set.seed(scitd_config$random_seed)
  truth_container <- scITD::run_tucker_ica(
    container = truth_container,
    ranks = c(scitd_config$donor_rank, scitd_config$gene_rank),
    tucker_type = scitd_config$tucker_type,
    rotation_type = scitd_config$rotation_type
  )

  message(
    "Truth tensor: ",
    paste(dim(truth_container$tensor_data[[4]]), collapse = " x ")
  )

  list(
    container = truth_container,
    cell_type_map = cell_type_map,
    negative_shifts = NULL
  )
}

scitd_truth_spec_from_result <- function(truth_result) {
  container <- truth_result$container
  cell_type_map <- truth_result$cell_type_map
  internal_cell_types <- container$tensor_data[[3]]

  if (!identical(internal_cell_types, cell_type_map$scitd_cell_type)) {
    stop("Truth result cell-type mapping does not match its tensor")
  }

  list(
    samples = container$tensor_data[[1]],
    genes = container$tensor_data[[2]],
    truth_cell_types = cell_type_map$truth_cell_type,
    scitd_cell_types = cell_type_map$scitd_cell_type
  )
}

read_scitd_truth_spec <- function(truth_result_dir) {
  sample_path <- file.path(truth_result_dir, "sample_scores.txt")
  loading_path <- file.path(
    truth_result_dir,
    "gene_celltype_loading_Factor1.txt"
  )
  if (!file.exists(sample_path) || !file.exists(loading_path)) {
    stop(
      "Saved truth scITD outputs require sample_scores.txt and ",
      "gene_celltype_loading_Factor1.txt in: ",
      truth_result_dir
    )
  }

  sample_scores <- read_scitd_expression_matrix(sample_path)
  factor_one <- read_scitd_expression_matrix(loading_path)
  cell_type_map <- scitd_cell_type_map(colnames(factor_one))

  list(
    samples = rownames(sample_scores),
    genes = rownames(factor_one),
    truth_cell_types = cell_type_map$truth_cell_type,
    scitd_cell_types = cell_type_map$scitd_cell_type
  )
}

discover_scitd_method_jobs <- function(
  dataset,
  deconv_config_ids,
  methods = NULL,
  lib_norm_methods = character(),
  repo_root = find_repo_root()
) {
  empty_jobs <- data.frame(
    dataset = character(),
    config_id = character(),
    config_slug = character(),
    method = character(),
    lib_norm = logical(),
    input_dir = character(),
    output_dir = character(),
    stringsAsFactors = FALSE
  )
  jobs <- list()

  for (config_id in deconv_config_ids) {
    paths <- tryCatch(
      resolve_deconv_paths(dataset, config_id, repo_root = repo_root),
      error = function(error) {
        message(
          "[SKIP] ", dataset, " / ", config_id, ": ",
          conditionMessage(error)
        )
        NULL
      }
    )
    if (is.null(paths)) {
      next
    }

    config_input_dir <- paths$output_dir
    if (!dir.exists(config_input_dir)) {
      message(
        "[SKIP] ", dataset, " / ", config_id,
        ": deconvolution result directory is unavailable: ",
        config_input_dir
      )
      next
    }

    method_dirs <- list.dirs(
      config_input_dir,
      recursive = FALSE,
      full.names = TRUE
    )
    if (length(method_dirs) == 0) {
      message(
        "[SKIP] ", dataset, " / ", config_id,
        ": no method result directories"
      )
      next
    }

    for (method_dir in method_dirs) {
      ctse_files <- list.files(
        method_dir,
        pattern = "\\.txt\\.gz$",
        full.names = TRUE
      )
      if (length(ctse_files) == 0) {
        message(
          "[SKIP] ", dataset, " / ", config_id, " / ",
          basename(method_dir), ": no CTSE .txt.gz files"
        )
        next
      }

      method <- basename(method_dir)
      if (
        !is.null(methods) &&
          !method %in% methods
      ) {
        next
      }
      lib_norm <- method %in% lib_norm_methods
      jobs[[length(jobs) + 1L]] <- data.frame(
        dataset = dataset,
        config_id = config_id,
        config_slug = basename(config_input_dir),
        method = method,
        lib_norm = lib_norm,
        input_dir = method_dir,
        output_dir = scitd_method_result_dir(
          dataset,
          config_id,
          method,
          lib_norm = lib_norm,
          repo_root = repo_root
        ),
        stringsAsFactors = FALSE
      )
    }
  }

  if (length(jobs) == 0) {
    return(empty_jobs)
  }
  do.call(rbind, jobs)
}

scitd_shift_negative_values <- function(
  matrices,
  tolerance = 1e-12
) {
  shift_rows <- list()

  for (cell_type in names(matrices)) {
    values <- matrices[[cell_type]]
    original_minimum <- apply(values, 1, min)
    shift_added <- pmax(0, -original_minimum)
    shifted_genes <- which(shift_added > 0)

    if (length(shifted_genes) > 0) {
      shift_rows[[length(shift_rows) + 1L]] <- data.frame(
        cell_type = rep(cell_type, length(shifted_genes)),
        gene = rownames(values)[shifted_genes],
        original_minimum = unname(original_minimum[shifted_genes]),
        shift_added = unname(shift_added[shifted_genes]),
        stringsAsFactors = FALSE,
        check.names = FALSE
      )
    }

    values <- sweep(values, 1, shift_added, FUN = "+")
    tiny_negative <- values < 0 & values >= -tolerance
    values[tiny_negative] <- 0
    if (any(values < 0)) {
      stop("Negative values remain after shifting cell type: ", cell_type)
    }
    matrices[[cell_type]] <- values
  }

  shifts <- if (length(shift_rows) == 0) {
    data.frame(
      cell_type = character(),
      gene = character(),
      original_minimum = numeric(),
      shift_added = numeric(),
      stringsAsFactors = FALSE
    )
  } else {
    do.call(rbind, shift_rows)
  }
  rownames(shifts) <- NULL

  list(matrices = matrices, shifts = shifts)
}

scitd_method_preprocess_mode <- function(method, method_preprocessing) {
  required_columns <- c("method", "preprocess_mode")
  missing_columns <- setdiff(required_columns, colnames(method_preprocessing))
  if (length(missing_columns) > 0) {
    stop(
      "scITD method preprocessing is missing columns: ",
      paste(missing_columns, collapse = ", ")
    )
  }

  mode_rows <- which(method_preprocessing$method == method)
  if (length(mode_rows) == 0L) {
    stop("No scITD preprocessing mapping for method: ", method)
  }
  if (length(mode_rows) > 1L) {
    stop(
      "Multiple scITD preprocessing mappings found for method: ",
      method
    )
  }

  preprocess_mode <- as.character(
    method_preprocessing$preprocess_mode[[mode_rows]]
  )
  supported_modes <- c("log1p", "already_transformed")
  if (!preprocess_mode %in% supported_modes) {
    stop(
      "Unsupported scITD preprocessing mode for ", method, ": ",
      preprocess_mode
    )
  }
  preprocess_mode
}

scitd_preprocess_inferred_container <- function(
  container,
  preprocess_mode,
  scale_factor,
  lib_norm = FALSE
) {
  if (
    length(lib_norm) != 1L ||
      !is.logical(lib_norm) ||
      is.na(lib_norm)
  ) {
    stop("lib_norm must be one non-missing logical value")
  }

  if (lib_norm && !identical(preprocess_mode, "log1p")) {
    stop(
      "lib_norm=true is supported only for preprocess_mode=log1p"
    )
  }

  if (identical(preprocess_mode, "log1p")) {
    if (lib_norm) {
      return(scITD::normalize_pseudobulk(
        container = container,
        method = "regular",
        scale_factor = scale_factor
      ))
    }

    for (cell_type in names(container$scMinimal_ctype)) {
      cell_type_object <- container$scMinimal_ctype[[cell_type]]
      transformed <- log1p(cell_type_object$pseudobulk)
      cell_type_object$pseudobulk <- Matrix::Matrix(
        Matrix::t(transformed),
        sparse = TRUE
      )
    }
    return(container)
  }

  if (identical(preprocess_mode, "already_transformed")) {
    for (cell_type in names(container$scMinimal_ctype)) {
      cell_type_object <- container$scMinimal_ctype[[cell_type]]
      cell_type_object$pseudobulk <- Matrix::Matrix(
        Matrix::t(cell_type_object$pseudobulk),
        sparse = TRUE
      )
    }
    return(container)
  }

  stop("Unsupported inferred CTSE preprocessing mode: ", preprocess_mode)
}

scitd_normalized_variance_qc <- function(container, genes, cell_types) {
  invalid_by_cell_type <- stats::setNames(
    vector("list", length(cell_types)),
    cell_types
  )

  for (cell_type in cell_types) {
    cell_type_object <- container$scMinimal_ctype[[cell_type]]
    normalized_variances <- cell_type_object$norm_variances
    if (is.null(normalized_variances) || length(normalized_variances) == 0) {
      selected_variances <- rep(NA_real_, length(genes))
      names(selected_variances) <- genes
    } else {
      if (!is.numeric(normalized_variances)) {
        stop("Normalized variances are not numeric in: ", cell_type)
      }
      selected_variances <- normalized_variances[genes]
    }

    invalid <- !is.finite(selected_variances) | selected_variances < 0
    invalid_by_cell_type[[cell_type]] <- genes[invalid]
  }

  invalid_values <- unique(unlist(
    invalid_by_cell_type,
    use.names = FALSE
  ))
  invalid_genes <- genes[genes %in% invalid_values]

  list(
    invalid_by_cell_type = invalid_by_cell_type,
    invalid_genes = invalid_genes,
    retained_genes = genes[!genes %in% invalid_genes]
  )
}

run_scitd_inferred <- function(
  dataset,
  config_id,
  method,
  truth_spec,
  scitd_config,
  method_preprocessing,
  repo_root = find_repo_root(),
  lib_norm = FALSE,
  input_dir = NULL
) {
  if (!requireNamespace("scITD", quietly = TRUE)) {
    stop("The scITD package is not available on .libPaths()")
  }
  if (!requireNamespace("Matrix", quietly = TRUE)) {
    stop("The Matrix package is required")
  }
  if (!config_id %in% scitd_config$deconv_config_ids) {
    stop("Config is not enabled for scITD: ", config_id)
  }

  preprocess_mode <- scitd_method_preprocess_mode(
    method,
    method_preprocessing
  )

  if (is.null(input_dir)) {
    paths <- resolve_deconv_paths(dataset, config_id, repo_root = repo_root)
    input_dir <- file.path(paths$output_dir, method)
  }
  if (!dir.exists(input_dir)) {
    stop("Method result directory is unavailable: ", input_dir)
  }

  message("=== scITD inferred: ", dataset, " / ", config_id, " / ", method, " ===")

  cell_type_map <- scitd_cell_type_map(truth_spec$truth_cell_types)
  if (!identical(cell_type_map$scitd_cell_type, truth_spec$scitd_cell_types)) {
    stop("Truth specification has an inconsistent cell-type mapping")
  }

  method_paths <- file.path(
    input_dir,
    paste0(cell_type_map$truth_cell_type, ".txt.gz")
  )
  names(method_paths) <- cell_type_map$truth_cell_type
  missing_paths <- method_paths[!file.exists(method_paths)]
  if (length(missing_paths) > 0) {
    stop(
      "Required truth-space cell types are unavailable: ",
      paste(names(missing_paths), collapse = ", ")
    )
  }

  inferred_values <- lapply(
    method_paths,
    read_scitd_expression_matrix,
    require_nonnegative = FALSE,
    require_integer = FALSE
  )

  common_method_genes <- scitd_ordered_common_values(
    lapply(inferred_values, rownames)
  )
  candidate_genes <- truth_spec$genes[
    truth_spec$genes %in% common_method_genes
  ]
  if (length(candidate_genes) == 0) {
    stop("No truth-selected genes are shared across method cell types")
  }

  samples_used <- scitd_ordered_common_values(c(
    list(truth_spec$samples),
    lapply(inferred_values, colnames)
  ))
  if (length(samples_used) == 0) {
    stop("No truth samples are shared across method cell types")
  }

  inferred_values <- lapply(inferred_values, function(values) {
    values[common_method_genes, samples_used, drop = FALSE]
  })

  # Zero profiles and negative shifts are defined over the common method gene
  # universe used for preprocessing and normalized-variance fitting. The tensor
  # is restricted to the valid truth-selected intersection afterward.
  zero_before_shift <- scitd_zero_profiles(
    inferred_values,
    sample_order = samples_used,
    tolerance = scitd_config$zero_profile_tolerance
  )
  samples_used <- samples_used[
    !samples_used %in% zero_before_shift$samples
  ]
  if (length(samples_used) == 0) {
    stop("All shared samples have a zero profile before negative-value shifting")
  }
  inferred_values <- lapply(inferred_values, function(values) {
    values[, samples_used, drop = FALSE]
  })

  shifted <- scitd_shift_negative_values(
    inferred_values,
    tolerance = scitd_config$zero_profile_tolerance
  )
  inferred_values <- shifted$matrices
  negative_shifts <- shifted$shifts

  zero_after_shift <- scitd_zero_profiles(
    inferred_values,
    sample_order = samples_used,
    tolerance = scitd_config$zero_profile_tolerance
  )
  samples_used <- samples_used[
    !samples_used %in% zero_after_shift$samples
  ]
  if (length(samples_used) == 0) {
    stop("All shared samples have a zero profile after negative-value shifting")
  }
  inferred_values <- lapply(inferred_values, function(values) {
    values[, samples_used, drop = FALSE]
  })

  message(
    "Shared truth samples: ", length(truth_spec$samples),
    "; excluded before shift: ", length(zero_before_shift$samples),
    "; excluded after shift: ", length(zero_after_shift$samples),
    "; retained: ", length(samples_used)
  )
  message(
    "Truth-selected genes: ", length(truth_spec$genes),
    "; method intersection: ", length(candidate_genes),
    "; preprocessing/variance universe: ", length(common_method_genes)
  )
  if (nrow(negative_shifts) > 0) {
    message(
      "Shifted negative gene/cell-type profiles: ",
      nrow(negative_shifts)
    )
  }

  internal_values <- stats::setNames(
    inferred_values[cell_type_map$truth_cell_type],
    cell_type_map$scitd_cell_type
  )
  inferred_params <- scITD::initialize_params(
    ctypes_use = cell_type_map$scitd_cell_type,
    ncores = scitd_config$n_cores,
    rand_seed = scitd_config$random_seed
  )
  inferred_container <- new.env(parent = emptyenv())
  inferred_container$experiment_params <- inferred_params
  inferred_container$scMinimal_ctype <- stats::setNames(
    vector("list", length(internal_values)),
    names(internal_values)
  )
  inferred_container$gn_convert <- NULL

  for (cell_type in names(internal_values)) {
    cell_type_object <- new.env(parent = emptyenv())
    cell_type_object$pseudobulk <- Matrix::Matrix(
      internal_values[[cell_type]],
      sparse = TRUE
    )
    inferred_container$scMinimal_ctype[[cell_type]] <- cell_type_object
  }

  inferred_container <- scitd_preprocess_inferred_container(
    container = inferred_container,
    preprocess_mode = preprocess_mode,
    scale_factor = scitd_config$scale_factor,
    lib_norm = lib_norm
  )
  if (identical(preprocess_mode, "log1p")) {
    if (lib_norm) {
      message(
        "Preprocessing mode: log1p; applied regular library normalization ",
        "to ", scitd_config$scale_factor, " followed by log1p"
      )
    } else {
      message(
        "Preprocessing mode: log1p; applied log1p without library ",
        "normalization or scale-factor rescaling"
      )
    }
  } else {
    message(
      "Preprocessing mode: already transformed; used shifted values ",
      "directly (transpose only, no normalization or additional log)"
    )
  }
  for (cell_type in cell_type_map$scitd_cell_type) {
    processed_values <- inferred_container$scMinimal_ctype[[cell_type]]$pseudobulk
    if (!identical(rownames(processed_values), samples_used)) {
      stop("Sample order differs after preprocessing in: ", cell_type)
    }
    if (!identical(colnames(processed_values), common_method_genes)) {
      stop("Gene order differs after preprocessing in: ", cell_type)
    }
    if (any(!is.finite(processed_values))) {
      stop("Non-finite values after preprocessing in: ", cell_type)
    }
  }

  inferred_container <- scITD::get_normalized_variance(inferred_container)
  for (cell_type in cell_type_map$scitd_cell_type) {
    cell_type_object <- inferred_container$scMinimal_ctype[[cell_type]]
    missing_genes <- setdiff(
      candidate_genes,
      colnames(cell_type_object$pseudobulk)
    )
    if (length(missing_genes) > 0) {
      stop("Candidate genes are missing after preprocessing in: ", cell_type)
    }
  }

  variance_qc <- scitd_normalized_variance_qc(
    container = inferred_container,
    genes = candidate_genes,
    cell_types = cell_type_map$scitd_cell_type
  )
  genes_used <- variance_qc$retained_genes
  invalid_counts <- lengths(variance_qc$invalid_by_cell_type)
  names(invalid_counts) <- cell_type_map$truth_cell_type[
    match(names(invalid_counts), cell_type_map$scitd_cell_type)
  ]
  message(
    "Normalized-variance QC: candidates=", length(candidate_genes),
    "; dropped=", length(variance_qc$invalid_genes),
    "; retained=", length(genes_used)
  )
  if (length(variance_qc$invalid_genes) > 0) {
    message(
      "Invalid normalized variance by cell type: ",
      paste0(names(invalid_counts), "=", invalid_counts, collapse = "; ")
    )
    example_genes <- head(variance_qc$invalid_genes, 10L)
    message(
      "Dropped invalid-variance genes",
      if (length(variance_qc$invalid_genes) > 10L) " (first 10)" else "",
      ": ", paste(example_genes, collapse = ", ")
    )
  }

  validate_scitd_ranks(
    n_samples = length(samples_used),
    n_genes = length(genes_used),
    config = scitd_config,
    label = paste(dataset, config_id, method, "after variance QC")
  )

  for (cell_type in cell_type_map$scitd_cell_type) {
    cell_type_object <- inferred_container$scMinimal_ctype[[cell_type]]
    cell_type_object$pseudobulk <- cell_type_object$pseudobulk[
      samples_used,
      genes_used,
      drop = FALSE
    ]
  }
  inferred_container$all_vargenes <- genes_used

  if (scitd_config$scale_var) {
    inferred_container <- scITD::scale_variance(
      container = inferred_container,
      var_scale_power = scitd_config$var_scale_power
    )
  }
  inferred_container <- scITD::stack_tensor(inferred_container)
  validate_scitd_tensor(
    inferred_container,
    samples = samples_used,
    genes = genes_used,
    cell_types = cell_type_map$scitd_cell_type
  )
  validate_scitd_cell_type_variation(
    inferred_container,
    cell_type_labels = stats::setNames(
      cell_type_map$truth_cell_type,
      cell_type_map$scitd_cell_type
    ),
    tolerance = scitd_config$zero_profile_tolerance
  )

  set.seed(scitd_config$random_seed)
  inferred_container <- scITD::run_tucker_ica(
    container = inferred_container,
    ranks = c(scitd_config$donor_rank, scitd_config$gene_rank),
    tucker_type = scitd_config$tucker_type,
    rotation_type = scitd_config$rotation_type
  )

  message(
    "Inferred tensor: ",
    paste(dim(inferred_container$tensor_data[[4]]), collapse = " x ")
  )

  list(
    container = inferred_container,
    cell_type_map = cell_type_map,
    negative_shifts = negative_shifts
  )
}

scitd_managed_output_files <- function(output_dir) {
  if (!dir.exists(output_dir)) {
    return(character())
  }
  c(
    file.path(output_dir, "sample_scores.txt"),
    file.path(output_dir, "negative_gene_shifts.txt"),
    list.files(
      output_dir,
      pattern = "^gene_celltype_loading_Factor[0-9]+\\.txt$",
      full.names = TRUE
    )
  )
}

write_scitd_table <- function(
  values,
  path,
  row.names = TRUE,
  col.names = NA
) {
  old_options <- options(digits = 17)
  on.exit(options(old_options), add = TRUE)
  write.table(
    values,
    file = path,
    sep = "\t",
    quote = FALSE,
    row.names = row.names,
    col.names = col.names
  )
}

write_scitd_outputs <- function(
  result,
  output_dir,
  overwrite = FALSE
) {
  container <- result$container
  cell_type_map <- result$cell_type_map
  samples <- container$tensor_data[[1]]
  genes <- container$tensor_data[[2]]
  internal_cell_types <- container$tensor_data[[3]]
  sample_scores <- as.matrix(container$tucker_results[[1]])
  factor_loadings <- as.matrix(container$tucker_results[[2]])

  if (!identical(rownames(sample_scores), samples)) {
    stop("Sample-score row names do not match tensor samples")
  }
  if (!identical(internal_cell_types, cell_type_map$scitd_cell_type)) {
    stop("Output cell-type mapping does not match tensor cell types")
  }
  if (nrow(factor_loadings) != ncol(sample_scores)) {
    stop("The numbers of score factors and loading factors differ")
  }
  if (ncol(factor_loadings) != length(genes) * length(internal_cell_types)) {
    stop("Loading dimensions do not match tensor genes and cell types")
  }
  expected_loading_names <- unlist(
    lapply(
      internal_cell_types,
      function(cell_type) paste(cell_type, genes, sep = ":")
    ),
    use.names = FALSE
  )
  if (!identical(colnames(factor_loadings), expected_loading_names)) {
    stop("scITD loading column order does not match cell-type-major tensor order")
  }
  if (any(!is.finite(sample_scores)) || any(!is.finite(factor_loadings))) {
    stop("scITD scores or loadings contain non-finite values")
  }

  factor_names <- paste0("Factor", seq_len(ncol(sample_scores)))
  colnames(sample_scores) <- factor_names
  rownames(factor_loadings) <- factor_names

  managed_files <- scitd_managed_output_files(output_dir)
  existing_files <- managed_files[file.exists(managed_files)]
  if (length(existing_files) > 0 && !overwrite) {
    stop(
      "Managed scITD outputs already exist in ", output_dir,
      ". Set overwrite_outputs <- TRUE to replace them."
    )
  }
  if (overwrite && length(existing_files) > 0) {
    unlink(existing_files, force = TRUE)
  }
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  write_scitd_table(
    sample_scores,
    file.path(output_dir, "sample_scores.txt"),
    row.names = TRUE,
    col.names = NA
  )

  for (factor_index in seq_len(nrow(factor_loadings))) {
    loading_matrix <- matrix(
      factor_loadings[factor_index, ],
      nrow = length(genes),
      ncol = length(internal_cell_types),
      byrow = FALSE,
      dimnames = list(genes, cell_type_map$truth_cell_type)
    )
    write_scitd_table(
      loading_matrix,
      file.path(
        output_dir,
        paste0(
          "gene_celltype_loading_Factor",
          factor_index,
          ".txt"
        )
      ),
      row.names = TRUE,
      col.names = NA
    )
  }

  negative_shifts <- result$negative_shifts
  if (!is.null(negative_shifts) && nrow(negative_shifts) > 0) {
    write_scitd_table(
      negative_shifts,
      file.path(output_dir, "negative_gene_shifts.txt"),
      row.names = FALSE,
      col.names = TRUE
    )
  }

  message("Saved scITD outputs: ", output_dir)
  invisible(output_dir)
}
