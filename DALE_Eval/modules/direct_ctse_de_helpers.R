# Native cell type-specific DE adapters, resumable bMIND chunks, runtime
# logging, and the non-cell-type-specific bulk limma baseline for config02.
# Direct TCA is intentionally absent: the assessment paper and the public TCA
# vignette do not specify a reproducible native CTS-DE call for that comparison.

direct_ctse_de_prepare_inputs <- function(dataset, config, context) {
  bulk_path <- file.path(context$paths$object_dir, "bulk_input", "sumcount.txt")
  fraction_path <- file.path(
    context$paths$object_dir,
    "frac_input",
    paste0(config$truth_fraction, ".txt")
  )
  bulk_counts <- ctse_de_read_matrix(
    bulk_path,
    require_nonnegative = TRUE,
    require_integer = TRUE
  )
  fraction <- ctse_de_read_matrix(
    fraction_path,
    require_nonnegative = TRUE
  )

  samples <- context$test_samples
  samples <- samples[samples %in% colnames(bulk_counts)]
  samples <- samples[samples %in% rownames(fraction)]
  if (length(samples) == 0L) {
    stop("No test samples overlap bulk counts and truth transcript fractions")
  }
  bulk_counts <- bulk_counts[, samples, drop = FALSE]
  fraction <- fraction[samples, , drop = FALSE]

  keep <- colSums(bulk_counts) > 0 & rowSums(abs(fraction)) > 0
  bulk_counts <- bulk_counts[, keep, drop = FALSE]
  fraction <- fraction[keep, , drop = FALSE]
  samples <- samples[keep]
  if (length(samples) == 0L) {
    stop("No direct-method samples have both bulk and fraction signal")
  }
  group_info <- ctse_de_group_info(samples, context)

  library_sizes <- colSums(bulk_counts)
  bulk_cpm <- sweep(bulk_counts, 2L, library_sizes, "/") * 1e6

  list(
    dataset = dataset,
    samples = samples,
    groups = group_info$groups,
    y = as.integer(group_info$groups == config$case_level),
    n_control = group_info$n_control,
    n_case = group_info$n_case,
    bulk_counts = bulk_counts,
    bulk_cpm = bulk_cpm,
    bulk_log_cpm = log2(1 + bulk_cpm),
    fraction = fraction,
    genes = rownames(bulk_counts),
    cell_types = colnames(fraction),
    n_test_samples = length(context$test_samples),
    bulk_path = bulk_path,
    fraction_path = fraction_path,
    bulk_md5 = unname(tools::md5sum(bulk_path)),
    fraction_md5 = unname(tools::md5sum(fraction_path))
  )
}

direct_bulk_de_prepare_inputs <- function(dataset, context) {
  bulk_path <- file.path(context$paths$object_dir, "bulk_input", "sumcount.txt")
  counts <- ctse_de_read_matrix(
    bulk_path,
    require_nonnegative = TRUE,
    require_integer = TRUE
  )
  samples <- context$test_samples[context$test_samples %in% colnames(counts)]
  if (length(samples) == 0L) {
    stop("No config02 test samples are present in bulk sum counts")
  }
  counts <- counts[, samples, drop = FALSE]
  keep <- colSums(counts) > 0
  counts <- counts[, keep, drop = FALSE]
  samples <- samples[keep]
  if (length(samples) == 0L) {
    stop("No config02 test samples have positive bulk library size")
  }
  group_info <- ctse_de_group_info(samples, context)
  covariates <- ctse_de_covariates_for_samples(context, samples)
  list(
    dataset = dataset,
    counts = counts,
    genes = rownames(counts),
    samples = samples,
    groups = group_info$groups,
    covariates = covariates,
    n_control = group_info$n_control,
    n_case = group_info$n_case,
    n_test_samples = length(context$test_samples),
    bulk_path = bulk_path
  )
}

# Prepare the config02 truth-fraction-adjusted bulk baseline from test samples.
direct_bulk_truthfrac_de_prepare_inputs <- function(dataset, config, context) {
  prepared_bulk <- direct_bulk_de_prepare_inputs(dataset, context)
  fraction_path <- file.path(
    context$paths$object_dir,
    "frac_input",
    paste0(config$truth_fraction, ".txt")
  )
  fraction <- ctse_de_read_matrix(
    fraction_path,
    require_nonnegative = TRUE
  )

  samples <- prepared_bulk$samples[
    prepared_bulk$samples %in% rownames(fraction)
  ]
  if (length(samples) == 0L) {
    stop("No bulk test samples overlap truth transcript fractions")
  }
  counts <- prepared_bulk$counts[, samples, drop = FALSE]
  fraction <- fraction[samples, , drop = FALSE]

  keep <- rowSums(abs(fraction)) > 0
  counts <- counts[, keep, drop = FALSE]
  fraction <- fraction[keep, , drop = FALSE]
  samples <- samples[keep]
  if (length(samples) == 0L) {
    stop("No bulk test samples have truth transcript-fraction signal")
  }
  group_info <- ctse_de_group_info(samples, context)
  covariates <- ctse_de_covariates_for_samples(context, samples)

  list(
    dataset = dataset,
    counts = counts,
    genes = rownames(counts),
    samples = samples,
    groups = group_info$groups,
    covariates = covariates,
    fraction = fraction,
    fraction_path = fraction_path,
    fraction_cell_types = colnames(fraction),
    n_control = group_info$n_control,
    n_case = group_info$n_case,
    n_test_samples = length(context$test_samples),
    bulk_path = prepared_bulk$bulk_path
  )
}

# Fractions sum to one, so the joint model uses an intercept and K-1 fraction
# columns. The omitted final cell type is the deterministic reference fraction.
direct_bulk_truthfrac_design <- function(groups, fraction, covariates) {
  if (length(groups) != nrow(fraction)) {
    stop("Truth-fraction rows do not match the disease-group vector")
  }
  if (nrow(covariates) != nrow(fraction)) {
    stop("DE covariate rows do not match truth-fraction rows")
  }
  if (!identical(rownames(covariates), rownames(fraction))) {
    stop("DE covariates and truth fractions are not in the same sample order")
  }
  if (is.null(colnames(fraction)) || anyDuplicated(colnames(fraction))) {
    stop("Truth fractions require unique cell-type column names")
  }

  fraction_columns <- if (ncol(fraction) > 1L) {
    seq_len(ncol(fraction) - 1L)
  } else {
    integer()
  }
  fraction_design <- fraction[, fraction_columns, drop = FALSE]
  if (ncol(fraction_design) > 0L) {
    colnames(fraction_design) <- make.unique(paste0(
      "truthfrac_",
      make.names(colnames(fraction_design))
    ))
  }
  model_data <- data.frame(
    fraction_design,
    covariates,
    disease_case_vs_control = as.integer(groups == levels(groups)[[2L]]),
    check.names = FALSE
  )
  rownames(model_data) <- rownames(fraction)
  design <- stats::model.matrix(~ ., data = model_data)
  if (qr(design)$rank < ncol(design)) {
    stop(
      "Truth-fraction, covariate, and disease design matrix is rank deficient: ",
      paste(colnames(design), collapse = ", ")
    )
  }

  list(
    design = design,
    coefficient = "disease_case_vs_control",
    reference_fraction = colnames(fraction)[[ncol(fraction)]],
    fraction_covariates = colnames(fraction_design),
    de_covariates = colnames(covariates)
  )
}

run_bulk_truthfrac_adjusted_limma <- function(prepared) {
  if (!requireNamespace("edgeR", quietly = TRUE)) {
    stop("Package 'edgeR' is required for truth-fraction-adjusted bulk limma")
  }
  if (!requireNamespace("limma", quietly = TRUE)) {
    stop("Package 'limma' is required for truth-fraction-adjusted bulk limma")
  }
  if (!identical(colnames(prepared$counts), rownames(prepared$fraction))) {
    stop("Bulk counts and truth fractions are not in the same sample order")
  }

  model <- direct_bulk_truthfrac_design(
    prepared$groups,
    prepared$fraction,
    prepared$covariates
  )
  dge <- edgeR::DGEList(
    counts = prepared$counts,
    group = prepared$groups
  )
  keep <- edgeR::filterByExpr(dge, group = prepared$groups)
  if (!any(keep)) {
    stop("No bulk genes pass filterByExpr for truth-fraction-adjusted limma")
  }
  dge <- dge[keep, , keep.lib.sizes = FALSE]
  dge <- edgeR::calcNormFactors(dge, method = "TMM")
  voom_values <- limma::voom(
    dge,
    design = model$design,
    plot = FALSE
  )
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
    result = ctse_de_expand_limma_result(top, prepared$genes),
    n_tested_genes = sum(keep),
    n_input_samples = length(prepared$samples),
    n_control = prepared$n_control,
    n_case = prepared$n_case,
    fraction_path = prepared$fraction_path,
    fraction_cell_types = prepared$fraction_cell_types,
    reference_fraction = model$reference_fraction,
    fraction_covariates = model$fraction_covariates,
    de_covariates = model$de_covariates,
    preprocessing = paste(
      "raw sum counts -> filterByExpr -> TMM -> voom ->",
      "intercept + K-1 truth transcript fractions + DE covariates + disease"
    ),
    effect_scale = "truth-fraction-adjusted voom log2-expression difference"
  )
}

direct_ctse_de_empty_result <- function(genes) {
  data.frame(
    gene = genes,
    effect_size = NA_real_,
    p_value = NA_real_,
    q_value = NA_real_,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

direct_ctse_de_expand <- function(effect, p_value, q_value, source_genes, all_genes) {
  result <- direct_ctse_de_empty_result(all_genes)
  matched <- match(source_genes, all_genes)
  if (anyNA(matched)) {
    stop("Native DE result genes do not match the bulk genes")
  }
  result$effect_size[matched] <- as.numeric(effect)
  result$p_value[matched] <- as.numeric(p_value)
  result$q_value[matched] <- as.numeric(q_value)
  result
}

direct_ctse_de_native_bmind_names <- function(cell_types) {
  names <- gsub("[^0-9A-Za-z///' ]", "", cell_types)
  names <- gsub(" ", "_", names, fixed = TRUE)
  if (anyDuplicated(names)) {
    stop("bMIND cell-type name cleanup creates duplicated names")
  }
  names
}

# -----------------------------------------------------------------------------
# Permanent runtime table for native cell type-specific DE methods only
# -----------------------------------------------------------------------------

direct_ctse_de_runtime_columns <- function() {
  c(
    "run_id", "dataset", "method", "config_id", "status", "start_time",
    "end_time", "elapsed_seconds", "n_core", "n_input_samples",
    "n_input_genes", "n_tested_genes", "chunk_size", "n_chunks",
    "n_completed_chunks", "n_failed_chunks", "problem"
  )
}

direct_ctse_de_read_runtime <- function(path) {
  columns <- direct_ctse_de_runtime_columns()
  if (!file.exists(path)) {
    empty <- matrix(character(), nrow = 0L, ncol = length(columns))
    colnames(empty) <- columns
    return(as.data.frame(empty, stringsAsFactors = FALSE))
  }
  runtime <- read.delim(
    path,
    sep = "\t",
    stringsAsFactors = FALSE,
    check.names = FALSE,
    colClasses = "character",
    na.strings = character()
  )
  ctse_de_require_columns(runtime, columns, "direct CTS-DE runtime table")
  runtime[, columns, drop = FALSE]
}

direct_ctse_de_write_runtime <- function(runtime, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temporary <- tempfile(pattern = ".runtime_", tmpdir = dirname(path))
  on.exit(unlink(temporary), add = TRUE)
  write.table(
    runtime[, direct_ctse_de_runtime_columns(), drop = FALSE],
    file = temporary,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE,
    na = "NA"
  )
  if (!file.rename(temporary, path)) {
    stop("Could not update direct CTS-DE runtime table: ", path)
  }
  invisible(path)
}

direct_ctse_de_start_runtime <- function(
  paths,
  dataset,
  method,
  config_id,
  n_core,
  n_input_samples,
  n_input_genes,
  chunk_size = NA_integer_,
  n_chunks = NA_integer_
) {
  start_clock <- Sys.time()
  stamp <- gsub(
    "[^0-9]", "",
    format(start_clock, "%Y%m%d%H%M%OS3", tz = "UTC")
  )
  run_id <- paste(method, config_id, stamp, Sys.getpid(), sep = "_")
  runtime_path <- file.path(paths$direct_output_dir, "runtime.txt")
  runtime <- direct_ctse_de_read_runtime(runtime_path)
  row <- data.frame(
    run_id = run_id,
    dataset = dataset,
    method = method,
    config_id = config_id,
    status = "running",
    start_time = format(start_clock, "%Y-%m-%d %H:%M:%OS3 %Z"),
    end_time = "",
    elapsed_seconds = NA_character_,
    n_core = as.character(as.integer(n_core)),
    n_input_samples = as.character(as.integer(n_input_samples)),
    n_input_genes = as.character(as.integer(n_input_genes)),
    n_tested_genes = NA_character_,
    chunk_size = if (is.na(chunk_size)) NA_character_ else as.character(as.integer(chunk_size)),
    n_chunks = if (is.na(n_chunks)) NA_character_ else as.character(as.integer(n_chunks)),
    n_completed_chunks = "0",
    n_failed_chunks = "0",
    problem = "",
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  runtime <- rbind(runtime, row)
  direct_ctse_de_write_runtime(runtime, runtime_path)
  list(run_id = run_id, path = runtime_path, start_clock = start_clock)
}

direct_ctse_de_finish_runtime <- function(
  entry,
  status,
  n_tested_genes = NA_integer_,
  n_completed_chunks = NA_integer_,
  n_failed_chunks = NA_integer_,
  problem = ""
) {
  runtime <- direct_ctse_de_read_runtime(entry$path)
  index <- which(runtime$run_id == entry$run_id)
  if (length(index) != 1L) {
    stop("Could not identify direct CTS-DE runtime row: ", entry$run_id)
  }
  end_clock <- Sys.time()
  runtime$status[index] <- status
  runtime$end_time[index] <- format(end_clock, "%Y-%m-%d %H:%M:%OS3 %Z")
  runtime$elapsed_seconds[index] <- as.character(round(
    as.numeric(difftime(end_clock, entry$start_clock, units = "secs")),
    3
  ))
  runtime$n_tested_genes[index] <- if (is.na(n_tested_genes)) {
    NA_character_
  } else {
    as.character(as.integer(n_tested_genes))
  }
  runtime$n_completed_chunks[index] <- if (is.na(n_completed_chunks)) {
    NA_character_
  } else {
    as.character(as.integer(n_completed_chunks))
  }
  runtime$n_failed_chunks[index] <- if (is.na(n_failed_chunks)) {
    NA_character_
  } else {
    as.character(as.integer(n_failed_chunks))
  }
  runtime$problem[index] <- gsub("[\t\r\n]+", " ", problem)
  direct_ctse_de_write_runtime(runtime, entry$path)
  invisible(runtime[index, , drop = FALSE])
}

# -----------------------------------------------------------------------------
# Resumable bMIND DE chunks
# -----------------------------------------------------------------------------

direct_ctse_de_prepare_bmind_plan <- function(prepared, chunk_size) {
  chunk_size <- as.integer(chunk_size)
  if (is.na(chunk_size) || chunk_size < 1L) {
    stop("bMIND DE chunk size must be a positive integer")
  }
  variances <- apply(prepared$bulk_log_cpm, 1L, stats::var)
  genes <- rownames(prepared$bulk_log_cpm)[!is.na(variances) & variances > 0]
  if (length(genes) == 0L) {
    stop("No bulk genes have nonzero variance for bMIND DE")
  }
  chunk_number <- ceiling(seq_len(length(genes)) / chunk_size)
  gene_chunks <- split(genes, chunk_number)
  names(gene_chunks) <- sprintf("chunk_%04d", seq_len(length(gene_chunks)))
  list(
    genes = genes,
    gene_chunks = gene_chunks,
    chunk_size = chunk_size,
    n_chunks = length(gene_chunks)
  )
}

direct_ctse_de_bmind_manifest <- function(prepared, plan, config, n_core) {
  data.frame(
    field = c(
      "dataset", "config_id", "control_level", "case_level", "chunk_size",
      "n_core", "n_genes", "n_samples", "n_cell_types", "bulk_md5",
      "fraction_md5", "MIND_version", "np"
    ),
    value = as.character(c(
      prepared$dataset,
      config$config_id,
      config$control_level,
      config$case_level,
      plan$chunk_size,
      n_core,
      length(plan$genes),
      length(prepared$samples),
      length(prepared$cell_types),
      prepared$bulk_md5,
      prepared$fraction_md5,
      as.character(utils::packageVersion("MIND")),
      "TRUE"
    )),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

direct_ctse_de_write_temp_table <- function(values, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temporary <- tempfile(pattern = ".bmind_", tmpdir = dirname(path))
  on.exit(unlink(temporary), add = TRUE)
  write.table(
    values,
    file = temporary,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE,
    na = "NA"
  )
  if (!file.rename(temporary, path)) {
    stop("Could not update bMIND temporary table: ", path)
  }
  invisible(path)
}

direct_ctse_de_bmind_status_template <- function(plan) {
  chunks <- names(plan$gene_chunks)
  rows <- lapply(chunks, function(chunk) {
    genes <- plan$gene_chunks[[chunk]]
    data.frame(
      chunk = chunk,
      first_gene = genes[[1L]],
      last_gene = genes[[length(genes)]],
      gene_count = length(genes),
      n_returned_genes = NA_integer_,
      n_missing_genes = NA_integer_,
      status = "pending",
      message = "",
      result_path = file.path("chunks", paste0(chunk, ".rds")),
      start_time = "",
      end_time = "",
      elapsed_seconds = NA_real_,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  })
  result <- do.call(rbind, rows)
  rownames(result) <- NULL
  result
}

direct_ctse_de_write_bmind_status <- function(status, chunk_dir) {
  direct_ctse_de_write_temp_table(
    status,
    file.path(chunk_dir, "chunk_status.txt")
  )
}

direct_ctse_de_read_bmind_table <- function(path) {
  read.delim(
    path,
    sep = "\t",
    stringsAsFactors = FALSE,
    check.names = FALSE,
    colClasses = "character",
    na.strings = "NA"
  )
}

direct_ctse_de_bmind_run_is_compatible <- function(
  chunk_dir,
  manifest,
  prepared,
  plan
) {
  required <- file.path(
    chunk_dir,
    c(
      "manifest.txt", "genes.txt", "samples.txt", "phenotype_y.txt",
      "cell_types.txt", "chunk_status.txt"
    )
  )
  if (!all(file.exists(required)) || !dir.exists(file.path(chunk_dir, "chunks"))) {
    return(FALSE)
  }

  observed_manifest <- direct_ctse_de_read_bmind_table(required[[1L]])
  if (
    !identical(observed_manifest$field, manifest$field) ||
      !identical(observed_manifest$value, manifest$value)
  ) {
    return(FALSE)
  }
  observed_genes <- direct_ctse_de_read_bmind_table(required[[2L]])
  observed_samples <- direct_ctse_de_read_bmind_table(required[[3L]])
  observed_y <- direct_ctse_de_read_bmind_table(required[[4L]])
  observed_cell_types <- direct_ctse_de_read_bmind_table(required[[5L]])
  observed_status <- direct_ctse_de_read_bmind_table(required[[6L]])
  expected_status <- direct_ctse_de_bmind_status_template(plan)

  identical(observed_genes$gene, plan$genes) &&
    identical(observed_samples$sample, prepared$samples) &&
    identical(observed_y$sample, prepared$samples) &&
    identical(observed_y$y, as.character(prepared$y)) &&
    identical(observed_cell_types$cell_type, prepared$cell_types) &&
    identical(observed_status$chunk, expected_status$chunk) &&
    identical(observed_status$first_gene, expected_status$first_gene) &&
    identical(observed_status$last_gene, expected_status$last_gene) &&
    identical(observed_status$gene_count, as.character(expected_status$gene_count)) &&
    identical(observed_status$result_path, expected_status$result_path)
}

direct_ctse_de_create_bmind_run <- function(
  chunk_dir,
  manifest,
  prepared,
  plan
) {
  dir.create(file.path(chunk_dir, "chunks"), recursive = TRUE, showWarnings = FALSE)
  direct_ctse_de_write_temp_table(manifest, file.path(chunk_dir, "manifest.txt"))
  direct_ctse_de_write_temp_table(
    data.frame(gene = plan$genes, stringsAsFactors = FALSE),
    file.path(chunk_dir, "genes.txt")
  )
  direct_ctse_de_write_temp_table(
    data.frame(sample = prepared$samples, stringsAsFactors = FALSE),
    file.path(chunk_dir, "samples.txt")
  )
  direct_ctse_de_write_temp_table(
    data.frame(
      sample = prepared$samples,
      y = prepared$y,
      stringsAsFactors = FALSE
    ),
    file.path(chunk_dir, "phenotype_y.txt")
  )
  direct_ctse_de_write_temp_table(
    data.frame(cell_type = prepared$cell_types, stringsAsFactors = FALSE),
    file.path(chunk_dir, "cell_types.txt")
  )
  direct_ctse_de_write_bmind_status(
    direct_ctse_de_bmind_status_template(plan),
    chunk_dir
  )
  invisible(chunk_dir)
}

direct_ctse_de_find_or_create_bmind_run <- function(
  paths,
  run_id,
  manifest,
  prepared,
  plan,
  resume = TRUE
) {
  temporary_root <- file.path(paths$direct_output_dir, "tmp_bMIND")
  dir.create(temporary_root, recursive = TRUE, showWarnings = FALSE)

  if (resume) {
    candidates <- list.dirs(temporary_root, recursive = FALSE, full.names = TRUE)
    if (length(candidates) > 0L) {
      modified <- file.info(candidates)$mtime
      candidates <- candidates[order(modified, decreasing = TRUE, na.last = TRUE)]
      for (candidate in candidates) {
        compatible <- tryCatch(
          direct_ctse_de_bmind_run_is_compatible(
            candidate, manifest, prepared, plan
          ),
          error = function(error) FALSE
        )
        if (compatible) {
          return(list(chunk_dir = candidate, resumed = TRUE))
        }
      }
    }
  }

  chunk_dir <- file.path(temporary_root, run_id)
  if (dir.exists(chunk_dir)) {
    stop("New bMIND temporary run directory already exists: ", chunk_dir)
  }
  direct_ctse_de_create_bmind_run(
    chunk_dir, manifest, prepared, plan
  )
  list(chunk_dir = chunk_dir, resumed = FALSE)
}

direct_ctse_de_read_bmind_chunk <- function(path, expected_genes) {
  if (!file.exists(path)) {
    stop("Missing bMIND chunk result: ", path)
  }
  result <- readRDS(path)
  if (!all(c("genes", "coef", "pval") %in% names(result))) {
    stop("Malformed bMIND chunk result: ", path)
  }
  result$coef <- as.matrix(result$coef)
  result$pval <- as.matrix(result$pval)
  if (
    !identical(as.character(result$genes), expected_genes) ||
      !identical(rownames(result$coef), expected_genes) ||
      !identical(rownames(result$pval), expected_genes)
  ) {
    stop("bMIND chunk genes do not match its manifest: ", path)
  }
  result
}

direct_ctse_de_save_bmind_chunk <- function(result, path) {
  temporary <- tempfile(pattern = ".chunk_", tmpdir = dirname(path))
  on.exit(unlink(temporary), add = TRUE)
  saveRDS(result, temporary)
  if (!file.rename(temporary, path)) {
    stop("Could not save bMIND chunk result: ", path)
  }
  invisible(path)
}

direct_ctse_de_bmind_error <- function(message, chunk_dir, status) {
  structure(
    list(
      message = message,
      call = NULL,
      chunk_dir = chunk_dir,
      n_tested_genes = sum(
        as.numeric(status$n_returned_genes[status$status == "completed"]),
        na.rm = TRUE
      ),
      n_completed_chunks = sum(status$status == "completed"),
      n_failed_chunks = sum(status$status == "failed")
    ),
    class = c("direct_ctse_de_bmind_error", "error", "condition")
  )
}

direct_ctse_de_run_bmind_chunk <- function(prepared, genes, n_core) {
  native <- MIND::bmind_de(
    bulk = prepared$bulk_log_cpm[genes, , drop = FALSE],
    frac = prepared$fraction,
    y = prepared$y,
    np = TRUE,
    ncore = n_core
  )
  if (!all(c("coef", "pval") %in% names(native))) {
    stop("bMIND DE chunk did not return coef and pval")
  }
  coefficient <- as.matrix(native$coef)
  p_value <- as.matrix(native$pval)
  if (is.null(rownames(coefficient)) || is.null(rownames(p_value))) {
    stop("bMIND DE chunk returned matrices without gene names")
  }
  common <- genes[
    genes %in% rownames(coefficient) & genes %in% rownames(p_value)
  ]
  expanded_coefficient <- matrix(
    NA_real_,
    nrow = length(genes),
    ncol = ncol(coefficient),
    dimnames = list(genes, colnames(coefficient))
  )
  expanded_p_value <- matrix(
    NA_real_,
    nrow = length(genes),
    ncol = ncol(p_value),
    dimnames = list(genes, colnames(p_value))
  )
  expanded_coefficient[common, ] <- coefficient[common, , drop = FALSE]
  expanded_p_value[common, ] <- p_value[common, , drop = FALSE]
  failed_genes <- rowSums(!is.na(expanded_p_value)) == 0
  expanded_coefficient[failed_genes, ] <- NA_real_
  list(
    genes = genes,
    coef = expanded_coefficient,
    pval = expanded_p_value,
    n_returned_genes = sum(rowSums(!is.na(expanded_p_value)) > 0)
  )
}

run_bmind_ctse_de <- function(
  prepared,
  plan,
  config,
  paths,
  n_core = 4L,
  resume = TRUE,
  run_id
) {
  if (!requireNamespace("MIND", quietly = TRUE)) {
    stop("Package 'MIND' is required for native bMIND DE")
  }
  manifest <- direct_ctse_de_bmind_manifest(prepared, plan, config, n_core)
  chunk_state <- direct_ctse_de_find_or_create_bmind_run(
    paths = paths,
    run_id = run_id,
    manifest = manifest,
    prepared = prepared,
    plan = plan,
    resume = resume
  )
  chunk_dir <- chunk_state$chunk_dir
  status <- direct_ctse_de_read_bmind_table(
    file.path(chunk_dir, "chunk_status.txt")
  )

  for (chunk in names(plan$gene_chunks)) {
    genes <- plan$gene_chunks[[chunk]]
    status_row <- match(chunk, status$chunk)
    result_path <- file.path(chunk_dir, status$result_path[[status_row]])

    if (identical(status$status[[status_row]], "completed")) {
      valid <- tryCatch(
        {
          direct_ctse_de_read_bmind_chunk(result_path, genes)
          TRUE
        },
        error = function(error) FALSE
      )
      if (valid) {
        message("[bMIND resume] keeping completed ", chunk)
        next
      }
    }

    start_clock <- Sys.time()
    status$status[[status_row]] <- "running"
    status$message[[status_row]] <- ""
    status$start_time[[status_row]] <- format(
      start_clock, "%Y-%m-%d %H:%M:%OS3 %Z"
    )
    status$end_time[[status_row]] <- ""
    status$elapsed_seconds[[status_row]] <- NA_character_
    direct_ctse_de_write_bmind_status(status, chunk_dir)
    message(
      "[bMIND chunk] ", chunk, " (", length(genes), " genes)"
    )

    chunk_error <- tryCatch(
      {
        result <- direct_ctse_de_run_bmind_chunk(prepared, genes, n_core)
        direct_ctse_de_save_bmind_chunk(result, result_path)
        status$n_returned_genes[[status_row]] <- result$n_returned_genes
        status$n_missing_genes[[status_row]] <-
          length(genes) - result$n_returned_genes
        status$status[[status_row]] <- "completed"
        status$message[[status_row]] <- ""
        NULL
      },
      error = function(error) error
    )
    end_clock <- Sys.time()
    status$end_time[[status_row]] <- format(
      end_clock, "%Y-%m-%d %H:%M:%OS3 %Z"
    )
    status$elapsed_seconds[[status_row]] <- round(
      as.numeric(difftime(end_clock, start_clock, units = "secs")),
      3
    )
    if (!is.null(chunk_error)) {
      status$status[[status_row]] <- "failed"
      status$message[[status_row]] <- gsub(
        "[\t\r\n]+", " ", conditionMessage(chunk_error)
      )
    }
    direct_ctse_de_write_bmind_status(status, chunk_dir)
    if (!is.null(chunk_error)) {
      stop(direct_ctse_de_bmind_error(
        paste0(chunk, " failed: ", conditionMessage(chunk_error)),
        chunk_dir,
        status
      ))
    }
  }

  status <- direct_ctse_de_read_bmind_table(
    file.path(chunk_dir, "chunk_status.txt")
  )
  if (any(status$status != "completed")) {
    stop(direct_ctse_de_bmind_error(
      "Not all bMIND chunks completed",
      chunk_dir,
      status
    ))
  }

  pieces <- lapply(names(plan$gene_chunks), function(chunk) {
    row <- match(chunk, status$chunk)
    direct_ctse_de_read_bmind_chunk(
      file.path(chunk_dir, status$result_path[[row]]),
      plan$gene_chunks[[chunk]]
    )
  })
  coefficient <- do.call(rbind, lapply(pieces, function(piece) piece$coef))
  p_value <- do.call(rbind, lapply(pieces, function(piece) piece$pval))
  coefficient <- coefficient[plan$genes, , drop = FALSE]
  p_value <- p_value[plan$genes, , drop = FALSE]
  q_value <- matrix(
    NA_real_,
    nrow = nrow(p_value),
    ncol = ncol(p_value),
    dimnames = dimnames(p_value)
  )
  for (cell_type_column in colnames(p_value)) {
    q_value[, cell_type_column] <- stats::p.adjust(
      p_value[, cell_type_column],
      method = "BH"
    )
  }

  native_names <- direct_ctse_de_native_bmind_names(prepared$cell_types)
  results <- list()
  summaries <- list()
  for (cell_type in prepared$cell_types) {
    native_cell_type <- native_names[prepared$cell_types == cell_type][[1L]]
    control_column <- paste0(native_cell_type, ":co")
    case_column <- paste0(native_cell_type, ":ca")
    if (
      !control_column %in% colnames(coefficient) ||
        !case_column %in% colnames(coefficient) ||
        !native_cell_type %in% colnames(p_value)
    ) {
      stop("Could not identify bMIND output columns for cell type: ", cell_type)
    }
    effect <- coefficient[, case_column] - coefficient[, control_column]
    result <- direct_ctse_de_expand(
      effect = effect,
      p_value = p_value[, native_cell_type],
      q_value = q_value[, native_cell_type],
      source_genes = plan$genes,
      all_genes = prepared$genes
    )
    results[[cell_type]] <- result
    summaries[[cell_type]] <- ctse_de_summary_row(
      dataset = prepared$dataset,
      source = "direct",
      method = "bMIND",
      cell_type = cell_type,
      status = "completed",
      n_input_samples = length(prepared$samples),
      n_signal_samples = length(prepared$samples),
      n_control = prepared$n_control,
      n_case = prepared$n_case,
      n_input_genes = length(prepared$genes),
      n_tested_genes = sum(!is.na(result$p_value)),
      preprocessing = "sum counts -> CPM -> log2(1+x); native bmind_de(np=TRUE)",
      effect_scale = "bMIND case coefficient - control coefficient"
    )
  }
  list(
    results = results,
    summaries = summaries,
    chunk_dir = chunk_dir,
    chunk_status = status,
    n_chunks = nrow(status),
    n_completed_chunks = sum(status$status == "completed"),
    n_failed_chunks = sum(status$status == "failed"),
    n_tested_genes = sum(rowSums(!is.na(p_value)) > 0),
    resumed = chunk_state$resumed
  )
}

direct_ctse_de_cleanup_bmind_chunks <- function(chunk_dir) {
  expected_root <- normalizePath(
    file.path(dirname(dirname(chunk_dir)), "tmp_bMIND"),
    mustWork = FALSE
  )
  target <- normalizePath(chunk_dir, mustWork = TRUE)
  if (!identical(dirname(target), expected_root)) {
    stop("Refusing to remove unexpected bMIND temporary directory: ", target)
  }
  unlink(target, recursive = TRUE, force = TRUE)
  if (dir.exists(target)) {
    stop("Could not remove completed bMIND temporary directory: ", target)
  }
  if (dir.exists(expected_root) && length(list.files(expected_root)) == 0L) {
    unlink(expected_root, recursive = TRUE, force = TRUE)
  }
  invisible(target)
}

# -----------------------------------------------------------------------------
# ENIGMA input preparation and native DE
# -----------------------------------------------------------------------------

direct_ctse_de_read_enigma_reference <- function(paths) {
  candidates <- file.path(paths$ref_dir, c("rowMeans_sig.csv", "rowMeans_sig.txt"))
  existing <- candidates[file.exists(candidates)]
  if (length(existing) != 1L) {
    stop("Expected one ENIGMA reference in: ", paths$ref_dir)
  }
  separator <- if (grepl("\\.csv$", existing)) "," else "\t"
  values <- read.delim(
    existing,
    sep = separator,
    row.names = 1L,
    check.names = FALSE
  )
  values <- as.matrix(values)
  storage.mode(values) <- "double"
  if (any(!is.finite(values))) {
    stop("ENIGMA reference contains non-finite values: ", existing)
  }
  values
}

direct_ctse_de_map_truth_fraction <- function(fraction, reference, paths, repo_root) {
  if (paths$refType != "indep") {
    missing <- setdiff(colnames(reference), colnames(fraction))
    if (length(missing) > 0L) {
      stop("Truth fractions missing ENIGMA cell types: ", paste(missing, collapse = ", "))
    }
    cell_types <- colnames(reference)
    return(list(
      fraction = fraction[, cell_types, drop = FALSE],
      reference = reference[, cell_types, drop = FALSE]
    ))
  }
  aligned <- map_truth_fraction_to_indep_ref(
    fraction,
    colnames(reference),
    paths,
    repo_root
  )
  list(
    fraction = aligned$frac,
    reference = reference[, aligned$reference_cell_types, drop = FALSE]
  )
}

direct_ctse_de_enigma_output_names <- function(cell_types, paths, repo_root) {
  if (paths$refType != "indep") {
    return(stats::setNames(cell_types, cell_types))
  }
  mapping <- read_indep_ref_cell_type_mapping(paths, repo_root)
  target_by_reference <- stats::setNames(
    mapping$target_cell_type,
    mapping$indep_ref_cell_type
  )
  output <- ifelse(
    cell_types %in% names(target_by_reference),
    target_by_reference[cell_types],
    cell_types
  )
  output <- unname(output)
  if (anyDuplicated(output)) {
    stop("ENIGMA cell-type mapping creates duplicated output names")
  }
  stats::setNames(output, cell_types)
}

direct_ctse_de_prepare_enigma_inputs <- function(
  prepared,
  config,
  repo_root = ctse_de_find_repo_root()
) {
  paths <- resolve_deconv_paths(prepared$dataset, config$config_id, repo_root)
  reference <- direct_ctse_de_read_enigma_reference(paths)
  aligned <- direct_ctse_de_map_truth_fraction(
    prepared$fraction,
    reference,
    paths,
    repo_root
  )
  genes <- intersect(prepared$genes, rownames(aligned$reference))
  if (length(genes) == 0L) {
    stop("No genes overlap sum-count bulk and ENIGMA reference")
  }
  list(
    common = prepared,
    genes = genes,
    bulk = prepared$bulk_cpm[genes, , drop = FALSE],
    fraction = aligned$fraction,
    reference = aligned$reference[genes, , drop = FALSE],
    output_names = direct_ctse_de_enigma_output_names(
      colnames(aligned$reference), paths, repo_root
    )
  )
}

run_enigma_ctse_de <- function(
  prepared,
  mode = c("L2", "trace")
) {
  if (!requireNamespace("ENIGMA", quietly = TRUE)) {
    stop("Package 'ENIGMA' is required for native ENIGMA DE")
  }
  mode <- match.arg(mode)
  method <- if (mode == "L2") "ENIGMAL2" else "ENIGMAtrace"
  common <- prepared$common

  object <- ENIGMA::create_ENIGMA(
    bulk = prepared$bulk,
    ref = prepared$reference,
    ref_type = "aggre"
  )
  object@result_cell_proportion <- prepared$fraction
  if (mode == "L2") {
    object <- ENIGMA::ENIGMA_L2_max_norm(
      object,
      alpha = 0.1,
      model_tracker = FALSE,
      model_name = "log",
      preprocess = "log"
    )
  } else {
    object <- ENIGMA::ENIGMA_trace_norm(
      object,
      alpha = 0.1,
      model_tracker = FALSE,
      model_name = "log",
      preprocess = "log"
    )
  }
  native <- ENIGMA::FindCSE_DEG(
    object,
    y = common$y,
    FDR_control = TRUE,
    FoldChange = FALSE
  )

  results <- list()
  summaries <- list()
  for (reference_cell_type in names(native)) {
    cell_type <- prepared$output_names[[reference_cell_type]]
    if (is.null(cell_type)) {
      stop("Unexpected ENIGMA result cell type: ", reference_cell_type)
    }
    table <- as.data.frame(native[[reference_cell_type]], check.names = FALSE)
    ctse_de_require_columns(
      table,
      c("ExpressionDifference", "pvalue", "qvalue"),
      paste(method, reference_cell_type, "DE result")
    )
    result <- direct_ctse_de_expand(
      effect = table$ExpressionDifference,
      p_value = table$pvalue,
      q_value = table$qvalue,
      source_genes = rownames(table),
      all_genes = common$genes
    )
    results[[cell_type]] <- result
    summaries[[cell_type]] <- ctse_de_summary_row(
      dataset = common$dataset,
      source = "direct",
      method = method,
      cell_type = cell_type,
      status = "completed",
      n_input_samples = length(common$samples),
      n_signal_samples = length(common$samples),
      n_control = common$n_control,
      n_case = common$n_case,
      n_input_genes = length(common$genes),
      n_tested_genes = sum(!is.na(result$p_value)),
      preprocessing = paste0(
        "sum counts -> CPM; ENIGMA ", mode,
        " preprocess=log; native FindCSE_DEG(FDR_control=TRUE)"
      ),
      effect_scale = "ENIGMA normalized CSE case mean - control mean"
    )
  }
  tested_genes <- unique(unlist(lapply(results, function(result) {
    result$gene[!is.na(result$p_value)]
  })))
  list(
    results = results,
    summaries = summaries,
    n_tested_genes = length(tested_genes)
  )
}

# -----------------------------------------------------------------------------
# Non-cell-type-specific bulk limma baseline (no runtime record)
# -----------------------------------------------------------------------------

run_bulk_limma_de <- function(prepared) {
  if (!requireNamespace("edgeR", quietly = TRUE)) {
    stop("Package 'edgeR' is required for the bulk limma baseline")
  }
  if (!requireNamespace("limma", quietly = TRUE)) {
    stop("Package 'limma' is required for the bulk limma baseline")
  }
  model <- ctse_de_design(prepared$groups, prepared$covariates)
  dge <- edgeR::DGEList(counts = prepared$counts, group = prepared$groups)
  keep <- edgeR::filterByExpr(dge, group = prepared$groups)
  if (!any(keep)) {
    stop("No bulk genes pass filterByExpr")
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
    result = ctse_de_expand_limma_result(top, prepared$genes),
    n_tested_genes = sum(keep),
    n_input_samples = length(prepared$samples),
    n_control = prepared$n_control,
    n_case = prepared$n_case,
    covariates = model$covariates
  )
}
