# Helpers for CTSE sample-mean and scITD normalized-variance evaluation.

sample_stats_enable_vendored_scitd <- function(repo_root = find_repo_root()) {
  library_root <- file.path(repo_root, "DALE_Eval", "external_modules")
  if (dir.exists(library_root)) {
    .libPaths(unique(c(normalizePath(library_root), .libPaths())))
  }
  invisible(.libPaths())
}

sample_stats_output_label <- function(method, lib_norm = FALSE) {
  if (
    length(lib_norm) != 1L ||
      !is.logical(lib_norm) ||
      is.na(lib_norm)
  ) {
    stop("lib_norm must be one non-missing logical value")
  }
  paste0(method, if (lib_norm) "_libnorm" else "")
}

sample_stats_method_performance_slug <- function(paths, expected_truth_type) {
  expected_truth_type <- trimws(as.character(expected_truth_type))
  if (
    length(expected_truth_type) != 1L ||
      is.na(expected_truth_type) ||
      !nzchar(expected_truth_type) ||
      !grepl("^[A-Za-z0-9_]+$", expected_truth_type)
  ) {
    stop("expected_truth_type must be one non-empty path-safe value")
  }
  paste0(
    config_slug(paths$config),
    "__truth-",
    expected_truth_type
  )
}

sample_stats_method_output_paths <- function(
  paths,
  method,
  expected_truth_type,
  lib_norm = FALSE
) {
  label <- sample_stats_output_label(method, lib_norm = lib_norm)
  root <- file.path(
    paths$obj_dir,
    "deconv_performance",
    sample_stats_method_performance_slug(paths, expected_truth_type)
  )
  list(
    sample_mean = file.path(root, "sample_mean", paste0(label, ".txt")),
    sample_logmean = file.path(
      root,
      "sample_logmean",
      paste0(label, ".txt")
    ),
    normalized_variance = file.path(
      root,
      "normalized_variance",
      paste0(label, ".txt")
    )
  )
}

sample_stats_truth_output_paths <- function(
  dataset,
  truth_type,
  repo_root = find_repo_root()
) {
  output_labels <- c(
    sumcount = "truth_sumcount",
    meancpm = "truth_meancpm",
    sumcount_cpm = "truth_sumcount_cpm"
  )
  if (!truth_type %in% names(output_labels)) {
    stop(
      "Unsupported truth_type=", truth_type,
      ". Supported truth types: ",
      paste(names(output_labels), collapse = ", ")
    )
  }
  root <- file.path(
    repo_root,
    "Benchmarking_obj",
    dataset,
    "deconv_performance",
    unname(output_labels[[truth_type]])
  )
  list(
    sample_mean = file.path(root, "sample_mean", "Z_truth.txt"),
    normalized_variance = file.path(
      root,
      "normalized_variance",
      "Z_truth.txt"
    )
  )
}

write_sample_stats_matrix <- function(x, path, digits = 6L) {
  if (
    length(digits) != 1L ||
      is.na(digits) ||
      digits < 0L ||
      digits != floor(digits)
  ) {
    stop("digits must be one non-negative integer")
  }
  x <- as.matrix(x)
  storage.mode(x) <- "double"
  if (
    nrow(x) == 0L ||
      ncol(x) == 0L ||
      is.null(rownames(x)) ||
      is.null(colnames(x))
  ) {
    stop("Sample-statistic output must have gene rows and cell-type columns")
  }
  finite <- is.finite(x)
  x[finite] <- round(x[finite], digits = as.integer(digits))
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  if (file.exists(path)) {
    warning("Overwriting existing sample-statistic file: ", path)
  }
  write.table(
    x,
    file = path,
    sep = "\t",
    quote = FALSE,
    col.names = NA
  )
}

read_sample_stats_matrices <- function(
  files,
  cell_types,
  require_nonnegative = FALSE,
  require_integer = FALSE,
  label = "CTSE"
) {
  missing <- setdiff(cell_types, names(files))
  if (length(missing) > 0L) {
    stop(label, " is missing cell-type files: ", paste(missing, collapse = ", "))
  }
  stats::setNames(
    lapply(cell_types, function(cell_type) {
      read_scitd_expression_matrix(
        files[[cell_type]],
        require_nonnegative = require_nonnegative,
        require_integer = require_integer
      )
    }),
    cell_types
  )
}

sample_stats_common_genes <- function(matrices, label = "CTSE") {
  genes <- sort(Reduce(intersect, lapply(matrices, rownames)))
  if (length(genes) == 0L) {
    stop(label, " has no genes shared across selected cell types")
  }
  genes
}

prepare_truth_sample_stats_input <- function(
  files,
  cell_types,
  test_samples,
  truth_type,
  label = paste("truth", truth_type)
) {
  supported_truth_types <- c("sumcount", "meancpm", "sumcount_cpm")
  if (!truth_type %in% supported_truth_types) {
    stop(
      "Unsupported truth_type=", truth_type,
      ". Supported truth types: ",
      paste(supported_truth_types, collapse = ", ")
    )
  }
  matrices <- read_sample_stats_matrices(
    files,
    cell_types,
    require_nonnegative = TRUE,
    require_integer = identical(truth_type, "sumcount"),
    label = label
  )
  missing_samples <- lapply(matrices, function(x) {
    setdiff(test_samples, colnames(x))
  })
  if (any(lengths(missing_samples) > 0L)) {
    problem_cell_types <- names(missing_samples)[lengths(missing_samples) > 0L]
    details <- paste0(
      problem_cell_types,
      " (",
      lengths(missing_samples)[problem_cell_types],
      ")"
    )
    stop(
      "Testing samples are missing from truth cell types: ",
      paste(details, collapse = ", ")
    )
  }
  genes <- sample_stats_common_genes(matrices, label = label)
  matrices <- lapply(matrices, function(x) {
    x[genes, test_samples, drop = FALSE]
  })
  list(matrices = matrices, genes = genes, samples = test_samples)
}

prepare_method_sample_stats_input <- function(
  files,
  cell_types,
  test_samples,
  label = "inferred CTSE"
) {
  matrices <- read_sample_stats_matrices(
    files,
    cell_types,
    require_nonnegative = FALSE,
    require_integer = FALSE,
    label = label
  )
  reference_samples <- colnames(matrices[[1L]])
  mismatched <- names(matrices)[!vapply(
    matrices,
    function(x) setequal(colnames(x), reference_samples),
    logical(1)
  )]
  if (length(mismatched) > 0L) {
    stop(
      label,
      " has non-identical sample sets across cell types: ",
      paste(mismatched, collapse = ", ")
    )
  }
  unexpected_samples <- setdiff(reference_samples, test_samples)
  if (length(unexpected_samples) > 0L) {
    stop(
      label,
      " contains samples not marked as testing samples: ",
      paste(unexpected_samples, collapse = ", ")
    )
  }
  samples <- test_samples[test_samples %in% reference_samples]
  if (length(samples) < 2L) {
    stop(label, " has fewer than two testing samples")
  }
  if (length(samples) != length(reference_samples)) {
    stop(label, " sample alignment unexpectedly changed the sample count")
  }
  genes <- sample_stats_common_genes(matrices, label = label)
  matrices <- lapply(matrices, function(x) {
    x[genes, samples, drop = FALSE]
  })
  list(matrices = matrices, genes = genes, samples = samples)
}

sample_stats_library_log1p <- function(
  values,
  scale_factor,
  norm_method = c("regular", "trim"),
  zero_tolerance = 1e-12
) {
  norm_method <- match.arg(norm_method)
  if (any(values < 0)) {
    stop("Library normalization requires nonnegative values")
  }
  library_sizes <- colSums(values)
  positive <- is.finite(library_sizes) & library_sizes > zero_tolerance
  transformed <- matrix(
    0,
    nrow = nrow(values),
    ncol = ncol(values),
    dimnames = dimnames(values)
  )
  if (any(positive)) {
    positive_values <- values[, positive, drop = FALSE]
    norm_factors <- rep(1, sum(positive))
    if (identical(norm_method, "trim") && sum(positive) > 1L) {
      if (!requireNamespace("edgeR", quietly = TRUE)) {
        stop("edgeR is required for TMM normalization of truth sum counts")
      }
      norm_factors <- edgeR::calcNormFactors(positive_values)
    }
    effective_sizes <- library_sizes[positive] * norm_factors
    if (any(!is.finite(effective_sizes) | effective_sizes <= 0)) {
      stop("Library normalization produced invalid effective library sizes")
    }
    normalized <- sweep(
      positive_values,
      MARGIN = 2L,
      STATS = effective_sizes,
      FUN = "/"
    )
    transformed[, positive] <- log1p(normalized * scale_factor)
  }
  transformed
}

preprocess_truth_sample_stats <- function(
  matrices,
  truth_type,
  scale_factor,
  norm_method = "trim",
  zero_tolerance = 1e-12
) {
  supported_truth_types <- c("sumcount", "meancpm", "sumcount_cpm")
  if (!truth_type %in% supported_truth_types) {
    stop(
      "Unsupported truth_type=", truth_type,
      ". Supported truth types: ",
      paste(supported_truth_types, collapse = ", ")
    )
  }
  zero_profiles <- scitd_zero_profiles(
    matrices,
    sample_order = colnames(matrices[[1L]]),
    tolerance = zero_tolerance
  )
  y <- lapply(matrices, function(values) {
    transformed <- if (identical(truth_type, "sumcount")) {
      sample_stats_library_log1p(
        values,
        scale_factor = scale_factor,
        norm_method = norm_method,
        zero_tolerance = zero_tolerance
      )
    } else {
      log1p(values)
    }
    Matrix::Matrix(Matrix::t(transformed), sparse = TRUE)
  })
  list(y = y, zero_profiles = zero_profiles$profiles)
}

preprocess_method_sample_stats <- function(
  matrices,
  preprocess_mode,
  scale_factor,
  lib_norm = FALSE,
  zero_tolerance = 1e-12
) {
  if (lib_norm && !identical(preprocess_mode, "log1p")) {
    stop("lib_norm=true is supported only for preprocess_mode=log1p")
  }
  shifted <- scitd_shift_negative_values(
    matrices,
    tolerance = zero_tolerance
  )
  zero_profiles <- scitd_zero_profiles(
    shifted$matrices,
    sample_order = colnames(shifted$matrices[[1L]]),
    tolerance = zero_tolerance
  )
  y <- lapply(shifted$matrices, function(values) {
    transformed <- if (identical(preprocess_mode, "log1p")) {
      if (lib_norm) {
        sample_stats_library_log1p(
          values,
          scale_factor = scale_factor,
          norm_method = "regular",
          zero_tolerance = zero_tolerance
        )
      } else {
        log1p(values)
      }
    } else if (identical(preprocess_mode, "already_transformed")) {
      values
    } else {
      stop("Unsupported inferred CTSE preprocessing mode: ", preprocess_mode)
    }
    Matrix::Matrix(Matrix::t(transformed), sparse = TRUE)
  })
  list(
    y = y,
    negative_shifts = shifted$shifts,
    zero_profiles = zero_profiles$profiles
  )
}

calculate_sample_logmean <- function(
  matrices,
  genes,
  cell_types,
  preprocess_mode,
  epsilon = 1e-10
) {
  if (
    length(epsilon) != 1L ||
      is.na(epsilon) ||
      !is.finite(epsilon) ||
      epsilon <= 0
  ) {
    stop("epsilon must be one positive finite value")
  }
  if (!preprocess_mode %in% c("log1p", "already_transformed")) {
    stop("Unsupported inferred CTSE preprocessing mode: ", preprocess_mode)
  }
  if (!identical(names(matrices), cell_types)) {
    stop("Raw CTSE cell-type order does not match the requested output order")
  }
  for (cell_type in cell_types) {
    values <- matrices[[cell_type]]
    if (
      ncol(values) < 2L ||
        !identical(rownames(values), genes) ||
        any(!is.finite(values))
    ) {
      stop("Invalid raw CTSE matrix for cell type: ", cell_type)
    }
  }

  raw_sample_mean <- do.call(cbind, lapply(cell_types, function(cell_type) {
    rowMeans(matrices[[cell_type]])[genes]
  }))
  rownames(raw_sample_mean) <- genes
  colnames(raw_sample_mean) <- cell_types

  mean_shifts <- stats::setNames(rep(0, length(genes)), genes)
  sample_logmean <- if (identical(preprocess_mode, "log1p")) {
    mean_shifts <- pmax(0, -apply(raw_sample_mean, 1L, min))
    names(mean_shifts) <- genes
    shifted_mean <- sweep(
      raw_sample_mean,
      MARGIN = 1L,
      STATS = mean_shifts,
      FUN = "+"
    )
    log2(pmax(shifted_mean, epsilon))
  } else {
    # These methods already export CTSE on a transformed/log-like scale, so
    # their arithmetic sample mean is retained without another logarithm.
    raw_sample_mean
  }
  dimnames(sample_logmean) <- dimnames(raw_sample_mean)

  shifted_genes <- names(mean_shifts)[mean_shifts > 0]
  list(
    sample_logmean = sample_logmean,
    negative_mean_shifts = data.frame(
      gene = shifted_genes,
      shift = unname(mean_shifts[shifted_genes]),
      stringsAsFactors = FALSE
    )
  )
}

sample_stats_nv_gam_diagnostic <- function(values, gam_k = 5L) {
  if (
    length(gam_k) != 1L ||
      is.na(gam_k) ||
      gam_k < 1L ||
      gam_k != floor(gam_k)
  ) {
    stop("gam_k must be one positive integer")
  }

  mean_variance <- scITD::colMeanVars(values, rowSel = NULL)
  log_mean <- suppressWarnings(log(mean_variance$m))
  log_variance <- suppressWarnings(log(mean_variance$v))
  usable <- is.finite(log_mean) &
    is.finite(log_variance) &
    is.finite(mean_variance$nobs) &
    mean_variance$nobs >= 0
  n_usable_genes <- sum(usable)
  n_unique_log_means <- length(unique(log_mean[usable]))

  reason <- NA_character_
  if (n_usable_genes < gam_k) {
    reason <- paste0(
      "fewer than ", gam_k,
      " genes have finite log-mean and log-variance"
    )
  } else if (n_unique_log_means < gam_k) {
    reason <- paste0(
      "fewer than ", gam_k,
      " distinct finite log-means are available"
    )
  }

  list(
    can_fit = is.na(reason),
    n_usable_genes = n_usable_genes,
    n_unique_log_means = n_unique_log_means,
    gam_k = as.integer(gam_k),
    reason = reason
  )
}

calculate_sample_mean_nv <- function(y, genes, cell_types) {
  if (!requireNamespace("scITD", quietly = TRUE)) {
    stop(
      "The scITD package is unavailable. Use an R environment with the ",
      "vendored scITD dependencies installed."
    )
  }
  if (!requireNamespace("Matrix", quietly = TRUE)) {
    stop("The Matrix package is required")
  }
  if (!identical(names(y), cell_types)) {
    stop("Y cell-type order does not match the requested output order")
  }
  for (cell_type in cell_types) {
    values <- y[[cell_type]]
    if (
      nrow(values) < 2L ||
        !identical(colnames(values), genes) ||
        any(!is.finite(values))
    ) {
      stop("Invalid Y matrix for cell type: ", cell_type)
    }
  }

  sample_mean <- do.call(cbind, lapply(cell_types, function(cell_type) {
    Matrix::colMeans(y[[cell_type]])[genes]
  }))
  rownames(sample_mean) <- genes
  colnames(sample_mean) <- cell_types

  normalized_variance <- matrix(
    NA_real_,
    nrow = length(genes),
    ncol = length(cell_types),
    dimnames = list(genes, cell_types)
  )
  nv_gam_skipped <- vector("list", length(cell_types))
  names(nv_gam_skipped) <- cell_types

  # scITD::norm_var_helper() uses gam.k = 5. Fit each cell type
  # independently so a degenerate profile produces an NA column without
  # preventing valid cell types or sample means from being returned.
  for (cell_type in cell_types) {
    diagnostic <- sample_stats_nv_gam_diagnostic(
      y[[cell_type]],
      gam_k = 5L
    )
    if (!diagnostic$can_fit) {
      nv_gam_skipped[[cell_type]] <- data.frame(
        cell_type = cell_type,
        n_usable_genes = diagnostic$n_usable_genes,
        n_unique_log_means = diagnostic$n_unique_log_means,
        gam_k = diagnostic$gam_k,
        reason = diagnostic$reason,
        stringsAsFactors = FALSE
      )
      next
    }

    container <- new.env(parent = emptyenv())
    container$experiment_params <- list(ctypes_use = cell_type)
    cell_type_object <- new.env(parent = emptyenv())
    cell_type_object$pseudobulk <- y[[cell_type]]
    container$scMinimal_ctype <- stats::setNames(
      list(cell_type_object),
      cell_type
    )
    container <- tryCatch(
      scITD::get_normalized_variance(container),
      error = function(e) {
        stop(
          "Normalized-variance calculation failed for cell type ",
          cell_type,
          ": ",
          conditionMessage(e),
          call. = FALSE
        )
      }
    )
    values <- container$scMinimal_ctype[[cell_type]]$norm_variances
    normalized_variance[, cell_type] <- values[genes]
  }

  nv_gam_skipped <- nv_gam_skipped[!vapply(
    nv_gam_skipped,
    is.null,
    logical(1)
  )]
  nv_gam_skipped <- if (length(nv_gam_skipped) == 0L) {
    data.frame(
      cell_type = character(0),
      n_usable_genes = integer(0),
      n_unique_log_means = integer(0),
      gam_k = integer(0),
      reason = character(0),
      stringsAsFactors = FALSE
    )
  } else {
    do.call(rbind, nv_gam_skipped)
  }
  rownames(nv_gam_skipped) <- NULL

  invalid_nv <- !is.finite(normalized_variance)
  normalized_variance[invalid_nv] <- NA_real_

  list(
    sample_mean = sample_mean,
    normalized_variance = normalized_variance,
    invalid_nv_by_cell_type = colSums(invalid_nv),
    nv_gam_skipped = nv_gam_skipped
  )
}
