# Helpers for config02 CTSE differential-expression benchmarking.
#
# The notebook-style step6 runner owns the dataset, method, and cell-type loops.
# Helpers stop on a missing or malformed input; step6 records the error and skips
# only that method/cell-type job.

ctse_de_find_repo_root <- function(start = getwd()) {
  current <- normalizePath(start, mustWork = TRUE)
  repeat {
    if (
      dir.exists(file.path(current, "DALE_Eval")) &&
        dir.exists(file.path(current, "Benchmarking_obj"))
    ) {
      return(current)
    }
    parent <- dirname(current)
    if (identical(parent, current)) {
      stop("Could not find the CTSE benchmark repository from: ", start)
    }
    current <- parent
  }
}

ctse_de_read_tsv <- function(path, col_classes = NA) {
  if (!file.exists(path)) {
    stop("Missing required table: ", path)
  }
  read.delim(
    path,
    sep = "\t",
    header = TRUE,
    stringsAsFactors = FALSE,
    check.names = FALSE,
    colClasses = col_classes
  )
}

ctse_de_require_columns <- function(values, required, label) {
  missing <- setdiff(required, colnames(values))
  if (length(missing) > 0L) {
    stop(label, " missing columns: ", paste(missing, collapse = ", "))
  }
  invisible(values)
}

read_ctse_de_configs <- function(
  repo_root = ctse_de_find_repo_root(),
  config_path = file.path(repo_root, "DALE_Eval", "configs", "ctse_de_configs.txt")
) {
  configs <- ctse_de_read_tsv(config_path, col_classes = "character")
  ctse_de_require_columns(
    configs,
    c(
      "dataset", "phenotype_col", "control_level", "case_level",
      "ctse_meta_cell_type_col", "bmind_de_chunk_size"
    ),
    "ctse_de_configs.txt"
  )
  configs
}

read_ctse_de_config <- function(
  dataset,
  repo_root = ctse_de_find_repo_root()
) {
  configs <- read_ctse_de_configs(repo_root)
  row <- configs[configs$dataset == dataset, , drop = FALSE]
  if (nrow(row) != 1L) {
    stop("Expected exactly one CTSE DE config row for: ", dataset)
  }
  config <- lapply(as.list(row[1L, , drop = FALSE]), function(value) value[[1L]])

  # Shared, fixed config02 choices are kept here instead of repeated in the TSV.
  config$config_id <- "config02"
  config$sample_split_relpath <- file.path("self_reference", "sample_split.txt")
  config$sample_split_group <- "test"
  config$ctse_meta_sample_col <- "sample"
  config$ctse_meta_ncells_col <- "nCells"
  config$min_cells <- 20L
  config$truth_type <- "sumcount"
  config$truth_fraction <- "truth_transcriptfrac"
  config$bmind_de_chunk_size <- as.integer(config$bmind_de_chunk_size)
  if (is.na(config$bmind_de_chunk_size) || config$bmind_de_chunk_size < 1L) {
    stop("bmind_de_chunk_size must be a positive integer for: ", dataset)
  }
  config$direct_folder <- "direct_ctsDEG_config02_bulk-sumcount"
  config
}

ctse_de_deconv_slug <- function(
  config_id = "config02",
  repo_root = ctse_de_find_repo_root()
) {
  configs <- ctse_de_read_tsv(
    file.path(repo_root, "DALE_Eval", "configs", "deconv_configs.txt"),
    col_classes = "character"
  )
  row <- configs[configs$config_id == config_id, , drop = FALSE]
  if (
    nrow(row) != 1L ||
      row$bulk_input[[1L]] != "sumcount" ||
      row$bulk_scale[[1L]] != "counts"
  ) {
    stop("CTSE DE requires config02 with sumcount/counts input")
  }
  normalization_suffix <- if (tolower(row$bulk_normalization[[1L]]) == "cpm") {
    ""
  } else {
    paste0("__norm-", tolower(row$bulk_normalization[[1L]]))
  }
  paste0(
    row$config_id[[1L]], "_bulk-", row$bulk_input[[1L]], "__frac-",
    row$frac_input[[1L]], "__ref-", row$refType[[1L]],
    normalization_suffix
  )
}

ctse_de_resolve_paths <- function(
  dataset,
  config,
  repo_root = ctse_de_find_repo_root()
) {
  object_dir <- file.path(repo_root, "Benchmarking_obj", dataset)
  if (!dir.exists(object_dir)) {
    stop("Missing benchmark object directory: ", object_dir)
  }
  deconv_slug <- ctse_de_deconv_slug(config$config_id, repo_root)
  list(
    dataset = dataset,
    object_dir = object_dir,
    sample_split = file.path(object_dir, config$sample_split_relpath),
    sample_meta = file.path(object_dir, "metadata", "sample_meta.txt"),
    ctse_meta = file.path(object_dir, "metadata", "ctse_meta.txt"),
    truth_dir = file.path(object_dir, "ctse_truth", config$truth_type),
    deconv_slug = deconv_slug,
    deconv_dir = file.path(object_dir, "deconv_res", deconv_slug),
    truth_output_dir = file.path(object_dir, "DE_res", "truth_sumcount"),
    inferred_output_dir = file.path(object_dir, "DE_res", deconv_slug),
    direct_output_dir = file.path(object_dir, "DE_res", config$direct_folder),
    baseline_output_dir = file.path(object_dir, "DE_res", "baseline")
  )
}

ctse_de_read_matrix <- function(
  path,
  require_nonnegative = FALSE,
  require_integer = FALSE
) {
  if (!file.exists(path)) {
    stop("Missing matrix: ", path)
  }
  connection <- if (grepl("\\.gz$", path, ignore.case = TRUE)) {
    gzfile(path, open = "rt")
  } else {
    file(path, open = "rt")
  }
  on.exit(close(connection), add = TRUE)
  values <- read.delim(
    connection,
    header = TRUE,
    row.names = 1L,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  values <- as.matrix(values)
  suppressWarnings(storage.mode(values) <- "double")
  if (nrow(values) == 0L || ncol(values) == 0L || any(!is.finite(values))) {
    stop("Matrix is empty, nonnumeric, or non-finite: ", path)
  }
  if (require_nonnegative && any(values < 0)) {
    stop("Matrix contains negative values: ", path)
  }
  if (require_integer && any(values != round(values))) {
    stop("Count matrix contains non-integer values: ", path)
  }
  if (require_integer) {
    values <- round(values)
  }
  values
}

ctse_de_truth_files <- function(
  dataset,
  config,
  repo_root = ctse_de_find_repo_root()
) {
  truth_dir <- ctse_de_resolve_paths(dataset, config, repo_root)$truth_dir
  files <- list.files(truth_dir, pattern = "\\.txt\\.gz$", full.names = TRUE)
  if (length(files) == 0L) {
    stop("No truth CTSE files in: ", truth_dir)
  }
  stats::setNames(files, sub("\\.txt\\.gz$", "", basename(files)))
}

ctse_de_method_files <- function(method_dir) {
  if (!dir.exists(method_dir)) {
    stop("Missing inferred CTSE method directory: ", method_dir)
  }
  files <- list.files(method_dir, pattern = "\\.txt\\.gz$", full.names = TRUE)
  if (length(files) == 0L) {
    stop("No inferred CTSE cell-type files in: ", method_dir)
  }
  stats::setNames(files, sub("\\.txt\\.gz$", "", basename(files)))
}

ctse_de_read_sample_context <- function(
  dataset,
  config,
  repo_root = ctse_de_find_repo_root()
) {
  paths <- ctse_de_resolve_paths(dataset, config, repo_root)
  split <- ctse_de_read_tsv(paths$sample_split, col_classes = "character")
  ctse_de_require_columns(split, c("group", "sampleIDs"), "sample_split.txt")
  test_samples <- unique(split$sampleIDs[
    tolower(split$group) == tolower(config$sample_split_group)
  ])
  test_samples <- test_samples[!is.na(test_samples) & nzchar(test_samples)]
  if (length(test_samples) == 0L) {
    stop("No testing samples for: ", dataset)
  }

  sample_meta <- ctse_de_read_tsv(paths$sample_meta, col_classes = "character")
  ctse_de_require_columns(
    sample_meta,
    c("sample", config$phenotype_col),
    "sample_meta.txt"
  )
  groups <- vapply(test_samples, function(sample) {
    observed <- unique(sample_meta[
      sample_meta$sample == sample,
      config$phenotype_col,
      drop = TRUE
    ])
    observed <- observed[!is.na(observed) & nzchar(observed)]
    if (length(observed) != 1L) {
      stop("Test sample does not have exactly one phenotype: ", sample)
    }
    observed
  }, character(1L))
  names(groups) <- test_samples
  if (!all(groups %in% c(config$control_level, config$case_level))) {
    stop("Testing samples contain phenotype levels outside the requested contrast")
  }

  list(
    dataset = dataset,
    config = config,
    paths = paths,
    test_samples = test_samples,
    groups = groups
  )
}

ctse_de_read_sample_covariates <- function(context, columns) {
  columns <- as.character(columns)
  if (length(columns) == 0L || anyNA(columns) || any(!nzchar(columns))) {
    stop("DE covariate columns must be a non-empty character vector")
  }
  if (anyDuplicated(columns)) {
    stop("Duplicated DE covariate columns: ", paste(columns, collapse = ", "))
  }

  sample_meta <- ctse_de_read_tsv(
    context$paths$sample_meta,
    col_classes = "character"
  )
  ctse_de_require_columns(
    sample_meta,
    c("sample", columns),
    "sample_meta.txt"
  )

  covariates <- data.frame(row.names = context$test_samples)
  for (column in columns) {
    values <- vapply(context$test_samples, function(sample) {
      observed <- unique(sample_meta[
        sample_meta$sample == sample,
        column,
        drop = TRUE
      ])
      observed <- observed[!is.na(observed) & nzchar(trimws(observed))]
      if (length(observed) != 1L) {
        stop(
          "Test sample does not have exactly one nonmissing ", column,
          " value: ", sample
        )
      }
      observed
    }, character(1L))
    covariates[[column]] <- factor(values)
  }
  covariates
}

ctse_de_covariates_for_samples <- function(context, samples) {
  if (is.null(context$covariates)) {
    return(data.frame(row.names = samples))
  }
  if (!all(samples %in% rownames(context$covariates))) {
    stop("DE covariates are missing retained analysis samples")
  }

  covariates <- context$covariates[samples, , drop = FALSE]
  for (column in colnames(covariates)) {
    covariates[[column]] <- droplevels(factor(covariates[[column]]))
    if (nlevels(covariates[[column]]) < 2L) {
      stop(
        "Configured DE covariate has fewer than two levels in retained samples: ",
        column
      )
    }
  }
  covariates
}

ctse_de_read_truth_meta <- function(context) {
  config <- context$config
  meta <- ctse_de_read_tsv(context$paths$ctse_meta, col_classes = "character")
  ctse_de_require_columns(
    meta,
    c(
      config$ctse_meta_cell_type_col,
      config$ctse_meta_sample_col,
      config$ctse_meta_ncells_col
    ),
    "ctse_meta.txt"
  )
  meta[[config$ctse_meta_ncells_col]] <- as.numeric(
    meta[[config$ctse_meta_ncells_col]]
  )
  if (anyNA(meta[[config$ctse_meta_ncells_col]])) {
    stop("ctse_meta.txt has invalid nCells values")
  }
  meta
}

ctse_de_group_info <- function(samples, context) {
  config <- context$config
  groups <- factor(
    context$groups[samples],
    levels = c(config$control_level, config$case_level)
  )
  counts <- table(groups)
  if (any(counts < 2L)) {
    stop(
      "Each contrast group needs at least two samples: ",
      paste(names(counts), as.integer(counts), sep = "=", collapse = ", ")
    )
  }
  list(
    groups = groups,
    n_control = unname(counts[[config$control_level]]),
    n_case = unname(counts[[config$case_level]])
  )
}

ctse_de_truth_eligibility <- function(cell_type, counts, context, ctse_meta) {
  config <- context$config
  cell_meta <- ctse_meta[
    ctse_meta[[config$ctse_meta_cell_type_col]] == cell_type,
    ,
    drop = FALSE
  ]
  if (anyDuplicated(cell_meta[[config$ctse_meta_sample_col]])) {
    stop("Duplicated ctse_meta rows for cell type: ", cell_type)
  }
  n_cells <- stats::setNames(
    cell_meta[[config$ctse_meta_ncells_col]],
    cell_meta[[config$ctse_meta_sample_col]]
  )
  samples <- context$test_samples
  samples <- samples[samples %in% names(n_cells)]
  samples <- samples[n_cells[samples] >= config$min_cells]
  samples <- samples[samples %in% colnames(counts)]
  samples <- samples[colSums(counts[, samples, drop = FALSE]) > 0]
  if (length(samples) == 0L) {
    stop("No eligible truth samples for cell type: ", cell_type)
  }
  c(list(samples = samples), ctse_de_group_info(samples, context))
}

ctse_de_design <- function(groups, covariates = NULL) {
  if (is.null(covariates)) {
    covariates <- data.frame(row.names = names(groups))
  }
  if (nrow(covariates) != length(groups)) {
    stop("DE covariate rows do not match the disease-group vector")
  }
  if (
    !is.null(names(groups)) &&
      !is.null(rownames(covariates)) &&
      !identical(names(groups), rownames(covariates))
  ) {
    stop("DE covariates and disease groups are not in the same sample order")
  }

  model_data <- data.frame(
    disease_case_vs_control = as.integer(groups == levels(groups)[[2L]]),
    check.names = FALSE
  )
  if (ncol(covariates) > 0L) {
    model_data <- cbind(model_data, covariates)
  }
  rownames(model_data) <- names(groups)
  design <- stats::model.matrix(~ ., data = model_data)
  if (qr(design)$rank < ncol(design)) {
    stop(
      "Disease plus covariate design matrix is rank deficient: ",
      paste(colnames(design), collapse = ", ")
    )
  }
  list(
    design = design,
    coefficient = "disease_case_vs_control",
    covariates = colnames(covariates)
  )
}

ctse_de_empty_limma_result <- function(genes) {
  data.frame(
    gene = genes,
    effect_size = NA_real_,
    t = NA_real_,
    p_value = NA_real_,
    q_value = NA_real_,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

ctse_de_expand_limma_result <- function(top, all_genes) {
  result <- ctse_de_empty_limma_result(all_genes)
  matched <- match(rownames(top), all_genes)
  result$effect_size[matched] <- top$logFC
  result$t[matched] <- top$t
  result$p_value[matched] <- top$P.Value
  result$q_value[matched] <- top$adj.P.Val
  result
}

ctse_de_fit_limma <- function(expression, groups, covariates = NULL) {
  if (!requireNamespace("limma", quietly = TRUE)) {
    stop("Package 'limma' is required for CTSE DE fitting")
  }
  testable <- apply(expression, 1L, stats::var) > 0
  if (!any(testable)) {
    stop("No genes have nonzero variance")
  }
  model <- ctse_de_design(groups, covariates)
  fit <- limma::lmFit(expression[testable, , drop = FALSE], model$design)
  fit <- limma::eBayes(fit, trend = TRUE, robust = TRUE)
  top <- limma::topTable(
    fit,
    coef = model$coefficient,
    number = Inf,
    sort.by = "none",
    adjust.method = "BH"
  )
  list(
    result = ctse_de_expand_limma_result(top, rownames(expression)),
    n_tested = sum(testable),
    covariates = model$covariates
  )
}

ctse_de_run_truth_cell_type <- function(
  dataset,
  cell_type,
  config,
  context,
  ctse_meta
) {
  if (!requireNamespace("edgeR", quietly = TRUE)) {
    stop("Package 'edgeR' is required for truth CTSE DE")
  }
  if (!requireNamespace("limma", quietly = TRUE)) {
    stop("Package 'limma' is required for truth CTSE DE")
  }
  counts <- ctse_de_read_matrix(
    file.path(context$paths$truth_dir, paste0(cell_type, ".txt.gz")),
    require_nonnegative = TRUE,
    require_integer = TRUE
  )
  eligibility <- ctse_de_truth_eligibility(
    cell_type, counts, context, ctse_meta
  )
  counts_used <- counts[, eligibility$samples, drop = FALSE]
  covariates <- ctse_de_covariates_for_samples(
    context,
    eligibility$samples
  )
  model <- ctse_de_design(eligibility$groups, covariates)

  dge <- edgeR::DGEList(counts = counts_used, group = eligibility$groups)
  keep <- edgeR::filterByExpr(dge, group = eligibility$groups)
  if (!any(keep)) {
    stop("No truth genes pass filterByExpr for: ", cell_type)
  }
  dge <- dge[keep, , keep.lib.sizes = FALSE]
  dge <- edgeR::calcNormFactors(dge, method = "TMM")
  voom_values <- limma::voom(dge, design = model$design, plot = FALSE)
  fit <- limma::lmFit(voom_values, model$design)
  fit <- limma::eBayes(fit, robust = TRUE)
  top <- limma::topTable(
    fit,
    coef = model$coefficient,
    number = Inf,
    sort.by = "none",
    adjust.method = "BH"
  )
  list(
    result = ctse_de_expand_limma_result(top, rownames(counts)),
    summary = ctse_de_summary_row(
      dataset = dataset,
      source = "truth",
      method = "truth",
      cell_type = cell_type,
      status = "completed",
      n_input_samples = ncol(counts),
      n_signal_samples = length(eligibility$samples),
      n_control = eligibility$n_control,
      n_case = eligibility$n_case,
      n_input_genes = nrow(counts),
      n_tested_genes = sum(keep),
      preprocessing = "TMM + voom",
      effect_scale = "voom log2-expression difference"
    )
  )
}

ctse_de_inferred_preprocessing_rules <- function() {
  c(
    InstaPrism = "log2p1_normalized",
    CIBERSORTx = "log2p1_normalized",
    BLUE = "native",
    scTAPE = "native",
    bMIND = "native",
    EPICunmix = "native",
    TCA = "already_transformed",
    ENIGMAL2 = "signed_native",
    ENIGMAtrace = "signed_native",
    Unico = "signed_native"
  )
}

ctse_de_detect_method_preprocessing <- function(method, method_dir) {
  files <- ctse_de_method_files(method_dir)
  contains_negative <- FALSE
  max_ctse_value <- -Inf

  for (cell_type in names(files)) {
    values <- ctse_de_read_matrix(files[[cell_type]])
    contains_negative <- contains_negative || any(values < 0)
    max_ctse_value <- max(max_ctse_value, max(values))
  }

  rules <- ctse_de_inferred_preprocessing_rules()
  if (!method %in% names(rules)) {
    stop(
      "No inferred DE preprocessing rule is configured for method '", method,
      "'. Add an explicit rule before running DE."
    )
  }
  preprocessing_id <- unname(rules[[method]])

  if (identical(preprocessing_id, "log2p1_normalized")) {
    if (contains_negative) {
      stop(
        method, " is configured for log2(1+x), but its inferred CTSE contains negative values"
      )
    }
    preprocessing <- "log2(1+x) + limma (already normalized)"
    effect_scale <- "log2(normalized CTSE + 1) difference"
  } else if (identical(preprocessing_id, "already_transformed")) {
    preprocessing <- "direct limma (TCA already transformed)"
    effect_scale <- "TCA native transformed-scale difference"
  } else if (identical(preprocessing_id, "signed_native")) {
    preprocessing <- "direct limma (signed native scale)"
    effect_scale <- "native signed-scale difference"
  } else if (identical(preprocessing_id, "native")) {
    preprocessing <- "direct limma (native scale)"
    effect_scale <- "native-scale difference"
  } else {
    stop("Unsupported inferred DE preprocessing rule for ", method, ": ", preprocessing_id)
  }

  list(
    method = method,
    files = files,
    contains_negative = contains_negative,
    max_ctse_value = max_ctse_value,
    preprocessing_id = preprocessing_id,
    preprocessing = preprocessing,
    effect_scale = effect_scale
  )
}

ctse_de_prepare_inferred_values <- function(values, specification) {
  if (identical(specification$preprocessing_id, "log2p1_normalized")) {
    if (any(values < 0)) {
      stop("log2(1+x) preprocessing received negative inferred CTSE values")
    }
    return(log2(1 + values))
  }
  values
}

ctse_de_run_inferred_cell_type <- function(
  dataset,
  method,
  cell_type,
  specification,
  context
) {
  input_path <- specification$files[[cell_type]]
  if (is.null(input_path) || is.na(input_path)) {
    stop("Missing inferred CTSE file for ", method, " / ", cell_type)
  }
  values <- ctse_de_read_matrix(input_path)
  n_input_samples <- ncol(values)
  n_input_genes <- nrow(values)
  samples <- context$test_samples[context$test_samples %in% colnames(values)]
  if (length(samples) == 0L) {
    stop("No test samples are present in inferred CTSE")
  }
  values <- values[, samples, drop = FALSE]

  signal <- if (specification$contains_negative) {
    colSums(abs(values)) > 0
  } else {
    colSums(values) > 0
  }
  values <- values[, signal, drop = FALSE]
  samples <- samples[signal]
  if (length(samples) == 0L) {
    stop("No inferred samples have cell-type signal")
  }
  group_info <- ctse_de_group_info(samples, context)
  covariates <- ctse_de_covariates_for_samples(context, samples)
  processed <- ctse_de_prepare_inferred_values(values, specification)
  fitted <- ctse_de_fit_limma(
    processed,
    group_info$groups,
    covariates
  )

  list(
    result = fitted$result,
    summary = ctse_de_summary_row(
      dataset = dataset,
      source = "inferred",
      method = method,
      cell_type = cell_type,
      status = "completed",
      n_input_samples = n_input_samples,
      n_signal_samples = length(samples),
      n_control = group_info$n_control,
      n_case = group_info$n_case,
      n_input_genes = n_input_genes,
      n_tested_genes = fitted$n_tested,
      contains_negative = specification$contains_negative,
      max_ctse_value = specification$max_ctse_value,
      preprocessing = specification$preprocessing,
      effect_scale = specification$effect_scale
    )
  )
}

ctse_de_summary_row <- function(
  dataset,
  source,
  method,
  cell_type = NA_character_,
  status,
  n_input_samples = NA_integer_,
  n_signal_samples = NA_integer_,
  n_control = NA_integer_,
  n_case = NA_integer_,
  n_input_genes = NA_integer_,
  n_tested_genes = NA_integer_,
  contains_negative = NA,
  max_ctse_value = NA_real_,
  preprocessing = NA_character_,
  effect_scale = NA_character_,
  problem = "",
  output_file = NA_character_
) {
  data.frame(
    dataset = dataset,
    source = source,
    method = method,
    cell_type = cell_type,
    status = status,
    n_input_samples = as.integer(n_input_samples),
    n_signal_samples = as.integer(n_signal_samples),
    n_control = as.integer(n_control),
    n_case = as.integer(n_case),
    n_input_genes = as.integer(n_input_genes),
    n_tested_genes = as.integer(n_tested_genes),
    contains_negative = as.logical(contains_negative),
    max_ctse_value = as.numeric(max_ctse_value),
    preprocessing = preprocessing,
    effect_scale = effect_scale,
    problem = problem,
    output_file = output_file,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

ctse_de_bind_summaries <- function(rows) {
  if (length(rows) == 0L) {
    return(data.frame())
  }
  result <- do.call(rbind, rows)
  rownames(result) <- NULL
  result
}

ctse_de_assert_output_path <- function(path, overwrite) {
  target <- normalizePath(path, mustWork = FALSE)
  benchmark_root <- normalizePath(
    file.path(ctse_de_find_repo_root(), "Benchmarking_obj"),
    mustWork = TRUE
  )
  if (!startsWith(target, paste0(benchmark_root, .Platform$file.sep))) {
    stop("Refusing output outside Benchmarking_obj: ", path)
  }
  if (file.exists(path) && !overwrite) {
    stop("Protected CTSE DE output already exists: ", path)
  }
  invisible(path)
}

ctse_de_write_table <- function(values, path, overwrite = FALSE) {
  ctse_de_assert_output_path(path, overwrite)
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temporary <- tempfile(pattern = ".ctse_de_", tmpdir = dirname(path))
  on.exit(unlink(temporary), add = TRUE)
  connection <- if (grepl("\\.gz$", path, ignore.case = TRUE)) {
    gzfile(temporary, open = "wt")
  } else {
    file(temporary, open = "wt")
  }
  write.table(
    values,
    file = connection,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE,
    na = "NA"
  )
  close(connection)
  if (file.exists(path)) {
    unlink(path)
  }
  if (!file.rename(temporary, path)) {
    stop("Could not move CTSE DE output into place: ", path)
  }
  invisible(path)
}

ctse_de_write_limma_result <- function(result, output_path, overwrite = FALSE) {
  columns <- c("gene", "effect_size", "t", "p_value", "q_value")
  ctse_de_require_columns(result, columns, "limma CTSE DE result")
  ctse_de_write_table(result[, columns, drop = FALSE], output_path, overwrite)
}

ctse_de_write_native_result <- function(result, output_path, overwrite = FALSE) {
  columns <- c("gene", "effect_size", "p_value", "q_value")
  ctse_de_require_columns(result, columns, "native CTSE DE result")
  ctse_de_write_table(result[, columns, drop = FALSE], output_path, overwrite)
}

ctse_de_safe_cell_type_filename <- function(cell_type) {
  paste0(gsub("/", "_", cell_type, fixed = TRUE), ".txt.gz")
}
