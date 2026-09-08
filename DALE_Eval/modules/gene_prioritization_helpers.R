# Helpers for truth-independent, config-level gene-prioritization scores.

gene_prioritization_score_names <- function() {
  c(
    "meanlog_margin",
    "meanlog_mean_contrast",
    "logmean_margin",
    "logmean_contrast",
    "partial_r2",
    "genesigtest_adapted_fdr"
  )
}

gene_prioritization_output_root <- function(paths) {
  file.path(
    paths$obj_dir,
    "gene_prioritization",
    config_slug(paths$config)
  )
}

gene_prioritization_method_output_paths <- function(
  paths,
  method,
  lib_norm = FALSE
) {
  label <- sample_stats_output_label(method, lib_norm = lib_norm)
  root <- gene_prioritization_output_root(paths)
  list(
    meanlog_margin = file.path(
      root,
      "post_hoc_meanlog_margin",
      paste0(label, ".txt.gz")
    ),
    meanlog_mean_contrast = file.path(
      root,
      "post_hoc_meanlog_mean_contrast",
      paste0(label, ".txt.gz")
    ),
    logmean_margin = file.path(
      root,
      "post_hoc_logmean_margin",
      paste0(label, ".txt.gz")
    ),
    logmean_contrast = file.path(
      root,
      "post_hoc_logmean_contrast",
      paste0(label, ".txt.gz")
    ),
    partial_r2 = file.path(
      root,
      "post_hoc_partial_r2",
      paste0(label, ".txt.gz")
    )
  )
}

gene_prioritization_genesigtest_output_path <- function(paths) {
  file.path(
    gene_prioritization_output_root(paths),
    "genesigtest_adapted_fdr",
    "GeneSigTest_adapted.txt.gz"
  )
}

read_gene_prioritization_matrix <- function(path, label) {
  if (!file.exists(path)) {
    stop(label, " not found: ", path)
  }
  x <- read.delim(
    path,
    sep = "\t",
    row.names = 1,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  x <- as.matrix(x)
  suppressWarnings(storage.mode(x) <- "double")
  if (
    nrow(x) == 0L ||
      ncol(x) == 0L ||
      is.null(rownames(x)) ||
      is.null(colnames(x))
  ) {
    stop(label, " must have gene rows and cell-type columns: ", path)
  }
  if (
    anyNA(rownames(x)) ||
      anyNA(colnames(x)) ||
      any(!nzchar(rownames(x))) ||
      any(!nzchar(colnames(x))) ||
      anyDuplicated(rownames(x)) ||
      anyDuplicated(colnames(x))
  ) {
    stop(label, " has invalid or duplicated dimension names: ", path)
  }
  if (any(!is.finite(x))) {
    stop(label, " contains non-finite values: ", path)
  }
  x
}

write_gene_prioritization_matrix <- function(
  x,
  path,
  digits = 6L,
  digit_mode = c("decimal", "significant")
) {
  digit_mode <- match.arg(digit_mode)
  if (
    length(digits) != 1L ||
      is.na(digits) ||
      digits < 1L ||
      digits != floor(digits)
  ) {
    stop("digits must be one positive integer")
  }
  x <- as.matrix(x)
  suppressWarnings(storage.mode(x) <- "double")
  if (
    nrow(x) == 0L ||
      ncol(x) == 0L ||
      is.null(rownames(x)) ||
      is.null(colnames(x))
  ) {
    stop("Gene-prioritization output must have gene rows and cell-type columns")
  }
  finite <- is.finite(x)
  if (identical(digit_mode, "significant")) {
    x[finite] <- signif(x[finite], digits = as.integer(digits))
  } else {
    x[finite] <- round(x[finite], digits = as.integer(digits))
  }
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  if (file.exists(path)) {
    warning("Overwriting existing gene-prioritization file: ", path)
  }
  old_options <- options(digits = 17L)
  on.exit(options(old_options), add = TRUE)
  connection <- if (grepl("\\.gz$", path, ignore.case = TRUE)) {
    gzfile(path, open = "wt")
  } else {
    file(path, open = "wt")
  }
  on.exit(try(close(connection), silent = TRUE), add = TRUE)
  write.table(
    x,
    file = connection,
    sep = "\t",
    quote = FALSE,
    col.names = NA
  )
}

calculate_meanlog_specificity <- function(sample_mean) {
  sample_mean <- as.matrix(sample_mean)
  suppressWarnings(storage.mode(sample_mean) <- "double")
  if (
    nrow(sample_mean) == 0L ||
      ncol(sample_mean) < 2L ||
      is.null(rownames(sample_mean)) ||
      is.null(colnames(sample_mean))
  ) {
    stop("sample_mean must contain genes and at least two cell types")
  }
  if (any(!is.finite(sample_mean))) {
    stop("sample_mean contains non-finite values")
  }

  margin <- matrix(
    NA_real_,
    nrow = nrow(sample_mean),
    ncol = ncol(sample_mean),
    dimnames = dimnames(sample_mean)
  )
  mean_contrast <- margin

  for (j in seq_len(ncol(sample_mean))) {
    other <- sample_mean[, -j, drop = FALSE]
    margin[, j] <- sample_mean[, j] - apply(other, 1L, max)
    mean_contrast[, j] <- sample_mean[, j] - rowMeans(other)
  }

  list(
    meanlog_margin = margin,
    meanlog_mean_contrast = mean_contrast
  )
}

calculate_logmean_specificity <- function(sample_logmean) {
  sample_logmean <- as.matrix(sample_logmean)
  suppressWarnings(storage.mode(sample_logmean) <- "double")
  if (
    nrow(sample_logmean) == 0L ||
      ncol(sample_logmean) < 2L ||
      is.null(rownames(sample_logmean)) ||
      is.null(colnames(sample_logmean))
  ) {
    stop("sample_logmean must contain genes and at least two cell types")
  }
  if (any(!is.finite(sample_logmean))) {
    stop("sample_logmean contains non-finite values")
  }

  margin <- matrix(
    NA_real_,
    nrow = nrow(sample_logmean),
    ncol = ncol(sample_logmean),
    dimnames = dimnames(sample_logmean)
  )
  contrast <- margin

  for (j in seq_len(ncol(sample_logmean))) {
    other <- sample_logmean[, -j, drop = FALSE]
    margin[, j] <- sample_logmean[, j] - apply(other, 1L, max)
    contrast[, j] <- sample_logmean[, j] - rowMeans(other)
  }

  list(
    logmean_margin = margin,
    logmean_contrast = contrast
  )
}

partial_r2_one_gene <- function(response, predictors) {
  response <- as.numeric(response)
  predictors <- as.matrix(predictors)
  suppressWarnings(storage.mode(predictors) <- "double")
  keep <- is.finite(response) & rowSums(!is.finite(predictors)) == 0L
  response <- response[keep]
  predictors <- predictors[keep, , drop = FALSE]
  n_predictors <- ncol(predictors)
  result <- rep(NA_real_, n_predictors)

  if (length(response) <= n_predictors + 1L) {
    return(result)
  }

  design <- cbind(`(Intercept)` = 1, predictors)
  fit <- lm.fit(design, response)
  residual_df <- length(response) - fit$rank
  if (residual_df <= 0L || fit$rank == 0L) {
    return(result)
  }

  rss <- sum(fit$residuals ^ 2)
  sigma_squared <- rss / residual_df
  rank_indices <- seq_len(fit$rank)
  r_matrix <- qr.R(fit$qr)[rank_indices, rank_indices, drop = FALSE]
  covariance_pivot <- tryCatch(
    chol2inv(r_matrix),
    error = function(e) NULL
  )
  if (is.null(covariance_pivot)) {
    return(result)
  }

  standard_error <- rep(NA_real_, ncol(design))
  estimable <- fit$qr$pivot[rank_indices]
  standard_error[estimable] <- sqrt(
    pmax(0, diag(covariance_pivot) * sigma_squared)
  )
  coefficients <- fit$coefficients

  for (j in seq_len(n_predictors)) {
    coefficient_index <- j + 1L
    coefficient <- coefficients[[coefficient_index]]
    se <- standard_error[[coefficient_index]]
    if (!is.finite(coefficient) || !is.finite(se)) {
      next
    }
    if (se == 0) {
      if (coefficient != 0) {
        result[[j]] <- 1
      }
      next
    }
    t_squared <- (coefficient / se) ^ 2
    result[[j]] <- t_squared / (t_squared + residual_df)
  }

  pmin(1, pmax(0, result))
}

calculate_post_hoc_partial_r2 <- function(
  y,
  bulk,
  n_core = 1L,
  min_positive_samples = 2L
) {
  if (
    length(n_core) != 1L ||
      is.na(n_core) ||
      n_core < 1L ||
      n_core != floor(n_core)
  ) {
    stop("n_core must be one positive integer")
  }
  if (
    length(min_positive_samples) != 1L ||
      is.na(min_positive_samples) ||
      min_positive_samples < 1L ||
      min_positive_samples != floor(min_positive_samples)
  ) {
    stop("min_positive_samples must be one positive integer")
  }
  if (length(y) < 2L || is.null(names(y)) || anyDuplicated(names(y))) {
    stop("y must be a named list containing at least two cell types")
  }

  bulk <- as.matrix(bulk)
  suppressWarnings(storage.mode(bulk) <- "double")
  if (
    nrow(bulk) == 0L ||
      ncol(bulk) == 0L ||
      is.null(rownames(bulk)) ||
      is.null(colnames(bulk)) ||
      any(!is.finite(bulk))
  ) {
    stop("bulk must be one finite gene-by-sample matrix")
  }

  reference_samples <- rownames(y[[1L]])
  reference_genes <- colnames(y[[1L]])
  invalid_y <- names(y)[!vapply(
    y,
    function(x) {
      identical(rownames(x), reference_samples) &&
        identical(colnames(x), reference_genes) &&
        all(is.finite(x))
    },
    logical(1)
  )]
  if (length(invalid_y) > 0L) {
    stop(
      "Preprocessed CTSE matrices are not finite and identically aligned: ",
      paste(invalid_y, collapse = ", ")
    )
  }

  samples <- reference_samples[reference_samples %in% colnames(bulk)]
  if (length(samples) <= length(y) + 1L) {
    stop(
      "Too few aligned samples for partial R2: n_samples=",
      length(samples),
      ", n_cell_types=",
      length(y)
    )
  }
  genes <- reference_genes[reference_genes %in% rownames(bulk)]
  bulk_use <- bulk[genes, samples, drop = FALSE]
  expressed <- rowSums(bulk_use > 0) >= as.integer(min_positive_samples)
  genes <- genes[expressed]
  bulk_use <- bulk_use[expressed, , drop = FALSE]
  if (length(genes) == 0L) {
    stop("No shared genes pass the bulk-expression filter for partial R2")
  }
  if (any(bulk_use < 0)) {
    stop("Partial R2 log1p bulk transformation requires nonnegative values")
  }
  bulk_use <- log1p(bulk_use)

  y_use <- lapply(y, function(x) x[samples, genes, drop = FALSE])
  cell_types <- names(y_use)
  calculate_indices <- function(indices) {
    result <- matrix(
      NA_real_,
      nrow = length(indices),
      ncol = length(cell_types),
      dimnames = list(genes[indices], cell_types)
    )
    for (position in seq_along(indices)) {
      gene_index <- indices[[position]]
      predictors <- do.call(cbind, lapply(
        y_use,
        function(x) as.numeric(x[, gene_index])
      ))
      result[position, ] <- partial_r2_one_gene(
        bulk_use[gene_index, ],
        predictors
      )
    }
    result
  }

  workers <- min(as.integer(n_core), length(genes))
  index_groups <- split(
    seq_along(genes),
    rep(seq_len(workers), length.out = length(genes))
  )
  result_parts <- if (workers > 1L && .Platform$OS.type != "windows") {
    parallel::mclapply(
      index_groups,
      calculate_indices,
      mc.cores = workers,
      mc.preschedule = TRUE
    )
  } else {
    lapply(index_groups, calculate_indices)
  }
  result <- do.call(rbind, result_parts)
  result[genes, cell_types, drop = FALSE]
}

gene_sig_test_adapted <- function(Bulk, frac, p_threshold = 0.05) {
  Bulk <- as.matrix(Bulk)
  frac <- as.matrix(frac)
  suppressWarnings(storage.mode(Bulk) <- "double")
  suppressWarnings(storage.mode(frac) <- "double")

  if (
    nrow(Bulk) == 0L ||
      ncol(Bulk) == 0L ||
      is.null(rownames(Bulk)) ||
      is.null(colnames(Bulk)) ||
      any(!is.finite(Bulk))
  ) {
    stop("Bulk must be one finite gene-by-sample matrix")
  }
  if (
    nrow(frac) == 0L ||
      ncol(frac) < 2L ||
      is.null(rownames(frac)) ||
      is.null(colnames(frac)) ||
      any(!is.finite(frac))
  ) {
    stop("frac must be one finite sample-by-cell-type matrix")
  }
  if (!identical(colnames(Bulk), rownames(frac))) {
    stop("Bulk columns and fraction rows must be identically aligned")
  }
  if (nrow(frac) <= ncol(frac)) {
    stop(
      "GeneSigTest_adapted requires more samples than cell types: n_samples=",
      nrow(frac),
      ", n_cell_types=",
      ncol(frac)
    )
  }
  if (
    length(p_threshold) != 1L ||
      !is.finite(p_threshold) ||
      p_threshold <= 0 ||
      p_threshold >= 1
  ) {
    stop("p_threshold must be strictly between zero and one")
  }

  design_qr <- qr(frac)
  if (design_qr$rank != ncol(frac)) {
    stop("GeneSigTest_adapted fraction design matrix is rank deficient")
  }
  response <- t(Bulk)
  fit <- lm.fit(frac, response)
  coefficients <- as.matrix(fit$coefficients)
  residuals <- as.matrix(fit$residuals)
  residual_df <- nrow(frac) - ncol(frac)

  rank_indices <- seq_len(design_qr$rank)
  r_matrix <- qr.R(design_qr)[rank_indices, rank_indices, drop = FALSE]
  covariance_pivot <- chol2inv(r_matrix)
  covariance_diagonal <- rep(NA_real_, ncol(frac))
  covariance_diagonal[design_qr$pivot[rank_indices]] <- diag(covariance_pivot)

  residual_variance <- colSums(residuals ^ 2) / residual_df
  standard_errors <- sqrt(outer(covariance_diagonal, residual_variance))
  t_statistics <- coefficients / standard_errors
  p_values <- stats::pt(
    t_statistics,
    df = residual_df,
    lower.tail = FALSE
  )
  p_values[is.nan(p_values)] <- 1
  p_values[!is.finite(p_values)] <- 1
  p_values <- t(p_values)
  rownames(p_values) <- rownames(Bulk)
  colnames(p_values) <- colnames(frac)

  adjusted <- apply(p_values, 2L, stats::p.adjust, method = "BH")
  adjusted <- as.matrix(adjusted)
  rownames(adjusted) <- rownames(Bulk)
  colnames(adjusted) <- colnames(frac)
  calls <- 1L * (adjusted < p_threshold)
  dimnames(calls) <- dimnames(adjusted)

  list(call = calls, pval = adjusted)
}

map_gene_prioritization_indep_columns <- function(
  x,
  paths,
  repo_root = find_repo_root()
) {
  x <- as.matrix(x)
  if (paths$refType != "indep") {
    return(x)
  }
  mapping <- read_indep_ref_cell_type_mapping(paths, repo_root)
  if (nrow(mapping) == 0L) {
    return(x)
  }
  duplicated_reference <- unique(
    mapping$indep_ref_cell_type[duplicated(mapping$indep_ref_cell_type)]
  )
  if (length(duplicated_reference) > 0L) {
    stop(
      "Duplicate independent-reference cell-type mappings: ",
      paste(duplicated_reference, collapse = ", ")
    )
  }
  target_by_reference <- stats::setNames(
    mapping$target_cell_type,
    mapping$indep_ref_cell_type
  )
  mapped <- ifelse(
    colnames(x) %in% names(target_by_reference),
    target_by_reference[colnames(x)],
    colnames(x)
  )
  if (anyDuplicated(mapped)) {
    stop(
      "Independent-reference mapping creates duplicated output names: ",
      paste(unique(mapped[duplicated(mapped)]), collapse = ", ")
    )
  }
  colnames(x) <- unname(mapped)
  x
}
