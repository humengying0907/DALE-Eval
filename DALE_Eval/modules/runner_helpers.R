############# Matrix readers #############

read_tsv_matrix <- function(path) {
  x <- read.delim(path, sep = "\t", check.names = FALSE, row.names = 1)
  x <- as.matrix(x)
  storage.mode(x) <- "numeric"
  x[is.na(x)] <- 0
  x
}

read_delim_matrix <- function(path, sep, label, row_names = TRUE) {
  if (!file.exists(path)) {
    stop(label, " not found: ", path)
  }
  x <- read.delim(
    path,
    sep = sep,
    check.names = FALSE,
    row.names = if (isTRUE(row_names)) 1 else NULL
  )
  x <- as.matrix(x)
  storage.mode(x) <- "numeric"
  x[is.na(x)] <- 0
  x
}

read_reference_matrix <- function(ref_dir, basenames, label) {
  path <- ref_file(ref_dir, basenames, label)
  sep <- if (grepl("\\.csv$", path)) "," else "\t"
  read_delim_matrix(path, sep = sep, label = label, row_names = TRUE)
}

read_bulk <- function(paths) {
  if (!file.exists(paths$bulk_path)) {
    stop("Bulk input not found: ", paths$bulk_path)
  }
  read_tsv_matrix(paths$bulk_path)
}

normalize_counts_to_cpm <- function(counts) {
  rows_keep <- rowSums(counts, na.rm = TRUE) > 0
  cols_keep <- colSums(counts, na.rm = TRUE) > 0
  out <- matrix(0, nrow = nrow(counts), ncol = ncol(counts), dimnames = dimnames(counts))
  if (!any(rows_keep) || !any(cols_keep)) {
    return(out)
  }
  counts_use <- counts[rows_keep, cols_keep, drop = FALSE]
  lib_sizes <- colSums(counts_use, na.rm = TRUE)
  out[rows_keep, cols_keep] <- sweep(counts_use, 2, lib_sizes, "/") * 1e6
  out
}

normalize_counts_with_edgeR <- function(counts, method, config_id) {
  if (!requireNamespace("edgeR", quietly = TRUE)) {
    stop("bulk_normalization=", method, " requires the edgeR package for config_id=", config_id)
  }

  rows_keep <- rowSums(counts, na.rm = TRUE) > 0
  cols_keep <- colSums(counts, na.rm = TRUE) > 0
  out <- matrix(0, nrow = nrow(counts), ncol = ncol(counts), dimnames = dimnames(counts))
  if (!any(rows_keep) || !any(cols_keep)) {
    return(out)
  }

  counts_use <- counts[rows_keep, cols_keep, drop = FALSE]
  dge <- edgeR::DGEList(counts = counts_use)
  edgeR_method <- if (method == "uq") "upperquartile" else "TMM"
  dge <- edgeR::calcNormFactors(dge, method = edgeR_method)
  if (any(is.na(dge$samples$norm.factors))) {
    stop("bulk_normalization=", method, " produced NA normalization factors for config_id=", config_id)
  }
  out[rows_keep, cols_keep] <- edgeR::cpm(dge, normalized.lib.sizes = TRUE)
  out
}

prepare_bulk_for_deconv <- function(bulk_expr, paths, method = "This method") {
  config <- paths$config
  bulk_scale <- trimws(tolower(config$bulk_scale[[1]]))
  bulk_normalization <- trimws(tolower(config$bulk_normalization[[1]]))
  config_id <- paths$config_id

  if (!bulk_scale %in% c("counts", "cpm")) {
    stop("Unsupported bulk_scale=", bulk_scale, " for config_id=", config_id)
  }
  if (!bulk_normalization %in% c("cpm", "tmm", "uq", "none")) {
    stop("Unsupported bulk_normalization=", bulk_normalization, " for config_id=", config_id)
  }

  action <- "used input as configured"
  out <- bulk_expr

  if (bulk_normalization == "none") {
    action <- "used input without runner normalization"
  } else if (bulk_normalization == "cpm") {
    if (bulk_scale == "counts") {
      out <- normalize_counts_to_cpm(bulk_expr)
      action <- "converted counts to CPM"
    } else {
      action <- "used CPM input without runner normalization"
    }
  } else {
    if (bulk_scale != "counts") {
      stop(
        method, " requires bulk_scale=counts for bulk_normalization=", bulk_normalization,
        " in config_id=", config_id
      )
    }
    out <- normalize_counts_with_edgeR(bulk_expr, bulk_normalization, config_id)
    action <- paste0("converted counts to ", bulk_normalization, "-normalized CPM")
  }

  list(
    bulk_expr = out,
    bulk_scale = bulk_scale,
    bulk_normalization = bulk_normalization,
    action = action,
    note = paste0(
      method, " bulk preparation; config_id=", config_id,
      " bulk_input=", config$bulk_input[[1]],
      " bulk_scale=", bulk_scale,
      " bulk_normalization=", bulk_normalization,
      "; ", action
    )
  )
}



############# Small formatting / validation helpers #############

format_deconv_dim <- function(bulk_expr) {
  paste0(nrow(bulk_expr), " genes x ", ncol(bulk_expr), " samples")
}

validate_fraction_cell_types <- function(frac, method = "This method") {
  cell_types <- colnames(frac)
  if (is.null(cell_types) || any(!nzchar(cell_types))) {
    stop(method, " fraction input must have non-empty cell-type column names")
  }
  duplicated_cell_types <- unique(cell_types[duplicated(cell_types)])
  if (length(duplicated_cell_types) > 0) {
    stop(method, " fraction input has duplicated cell types: ", paste(duplicated_cell_types, collapse = ", "))
  }
  invisible(cell_types)
}

validate_non_log_bulk_input <- function(bulk_expr, bulk_input, config_id, method = "This method") {
  finite_values <- bulk_expr[is.finite(bulk_expr)]
  if (length(finite_values) == 0) {
    stop(method, " bulk input has no finite values for config_id=", config_id)
  }

  max_bulk <- max(finite_values)
  if (max_bulk < 50) {
    stop(
      method, " only takes non-log bulk input. ",
      "config_id=", config_id,
      " uses bulk_input=", bulk_input,
      " with max(bulk_expr)=", round(max_bulk, 4),
      ". Use a non-log bulk_input such as wcpm or a normalized sumcount input."
    )
  }
}

prepare_log_bulk_input <- function(bulk_expr, bulk_input, config_id, method = "This method") {
  action <- "applied log2(x + 1)"
  log_bulk <- log2(bulk_expr + 1)

  list(
    bulk_expr = log_bulk,
    was_logged = FALSE,
    action = action,
    note = paste0(
      method, " takes logged bulk input; config_id=", config_id,
      " bulk_input=", bulk_input, "; ", action
    )
  )
}

validate_bulk_frac <- function(bulk, frac) {
  missing_in_frac <- setdiff(colnames(bulk), rownames(frac))
  if (length(missing_in_frac) > 0) {
    stop("Bulk samples missing from fraction matrix: ", paste(missing_in_frac, collapse = ", "))
  }
  frac[colnames(bulk), , drop = FALSE]
}



############# Fraction path resolution #############

resolve_instaprism_fraction_path <- function(paths, repo_root = find_repo_root()) {
  default_paths <- paths
  if (config_bulk_normalization(paths$config) != "cpm") {
    default_config <- find_deconv_config(
      bulk_input = paths$config$bulk_input[[1]],
      frac_input = "InstaPrismfrac",
      refType = paths$refType,
      bulk_normalization = "cpm",
      repo_root = repo_root
    )
    default_paths <- resolve_deconv_paths(paths$dataset, default_config$config_id[[1]], repo_root = repo_root)
  }

  same_config_path <- file.path(default_paths$output_dir, "InstaPrism", "InstaPrismfrac.txt")
  if (file.exists(same_config_path)) {
    return(same_config_path)
  }

  source_bulk_input <- default_paths$config$bulk_input[[1]]
  config_path <- file.path(repo_root, "DALE_Eval", "configs", "deconv_configs.txt")
  configs <- read.delim(config_path, sep = "	", stringsAsFactors = FALSE, check.names = FALSE)
  candidate_configs <- configs[
    configs$bulk_input == source_bulk_input &
      tolower(configs$bulk_normalization) == "cpm" &
      configs$refType == default_paths$refType,
    ,
    drop = FALSE
  ]
  candidate_paths <- vapply(
    candidate_configs$config_id,
    function(config_id) {
      source_paths <- resolve_deconv_paths(default_paths$dataset, config_id, repo_root = repo_root)
      file.path(source_paths$output_dir, "InstaPrism", "InstaPrismfrac.txt")
    },
    character(1)
  )
  existing_paths <- candidate_paths[file.exists(candidate_paths)]
  if (length(existing_paths) != 1) {
    stop(
      "Expected exactly one existing default-CPM InstaPrism fraction for ",
      "dataset=", default_paths$dataset,
      ", source_bulk_input=", source_bulk_input,
      ", refType=", default_paths$refType,
      "; found ", length(existing_paths),
      ". Checked: ", paste(candidate_paths, collapse = ", ")
    )
  }
  existing_paths[[1]]
}

resolve_frac_input_path <- function(paths, repo_root = find_repo_root()) {
  frac_input <- paths$config$frac_input[[1]]
  if (frac_input == "InstaPrismfrac") {
    return(resolve_instaprism_fraction_path(paths, repo_root = repo_root))
  }

  paths$frac_path
}

read_frac_for_method <- function(paths, repo_root = find_repo_root()) {
  frac_path <- resolve_frac_input_path(paths, repo_root = repo_root)
  if (!file.exists(frac_path)) {
    stop(
      "Fraction input not found for config_id=", paths$config_id,
      " frac_input=", paths$config$frac_input[[1]],
      ". Resolved path: ", frac_path
    )
  }
  read_tsv_matrix(frac_path)
}



############# Sample restriction #############

restrict_to_test_samples <- function(bulk, paths, frac = NULL, use_test_samples = TRUE) {
  if (!isTRUE(use_test_samples)) {
    return(list(
      bulk = bulk,
      frac = frac,
      use_test_samples = FALSE,
      sample_split_path = NA_character_,
      n_samples = ncol(bulk),
      selected_samples = colnames(bulk),
      note = "running all samples"
    ))
  }

  sample_split_path <- file.path(paths$obj_dir, "self_reference", "sample_split.txt")
  if (!file.exists(sample_split_path)) {
    stop("use_test_samples=true but sample split file not found: ", sample_split_path)
  }

  sample_split <- read.delim(sample_split_path, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
  required_cols <- c("group", "sampleIDs")
  missing_cols <- setdiff(required_cols, colnames(sample_split))
  if (length(missing_cols) > 0) {
    stop("sample_split.txt missing columns: ", paste(missing_cols, collapse = ", "))
  }

  test_samples <- sample_split$sampleIDs[tolower(sample_split$group) == "test"]
  test_samples <- unique(test_samples[!is.na(test_samples) & nzchar(test_samples)])
  selected_samples <- intersect(colnames(bulk), test_samples)
  if (!is.null(frac)) {
    selected_samples <- intersect(selected_samples, rownames(frac))
  }

  if (length(selected_samples) == 0) {
    overlap_targets <- "bulk samples"
    if (!is.null(frac)) {
      overlap_targets <- "bulk and fraction samples"
    }
    stop(
      "use_test_samples=true but no test samples overlap with ", overlap_targets,
      ". sample_split_path=", sample_split_path,
      "; n_test_samples=", length(test_samples)
    )
  }

  bulk <- bulk[, selected_samples, drop = FALSE]
  if (!is.null(frac)) {
    frac <- frac[selected_samples, , drop = FALSE]
  }

  list(
    bulk = bulk,
    frac = frac,
    use_test_samples = TRUE,
    sample_split_path = sample_split_path,
    n_samples = length(selected_samples),
    selected_samples = selected_samples,
    note = paste0("restricted to ", length(selected_samples), " test samples")
  )
}



############# Output writing #############

write_tsv_matrix <- function(x, path, digits = 2) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  idx <- is.finite(x) & x != 0
  x[idx] <- round(x[idx], digits)

  if (grepl("\\.gz$", path)) {
    con <- gzfile(path, open = "wt")
    on.exit(close(con), add = TRUE)
    file_arg <- con
  } else {
    file_arg <- path
  }

  write.table(
    x,
    file = file_arg,
    sep = "\t",
    quote = FALSE,
    col.names = NA
  )
}

safe_cell_type_filename <- function(cell_type) {
  paste0(gsub("/", "_", cell_type), ".txt.gz")
}

message_existing_method_output <- function(paths, method) {
  method_dir <- file.path(paths$output_dir, method)
  if (dir.exists(method_dir) && length(list.files(method_dir, all.files = FALSE, no.. = TRUE)) > 0) {
    message(
      "  overwrite note: existing ", method,
      " outputs found for this config; files in ", method_dir,
      " may be overwritten"
    )
  }
}

write_ctse_array <- function(z_array, output_dir, method, digits = 2) {
  method_dir <- file.path(output_dir, method)
  dir.create(method_dir, recursive = TRUE, showWarnings = FALSE)
  cell_types <- dimnames(z_array)[[3]]
  if (is.null(cell_types)) {
    stop("Cannot export CTSE array without cell-type dimnames on dimension 3")
  }
  for (cell_type in cell_types) {
    out_path <- file.path(method_dir, safe_cell_type_filename(cell_type))
    write_tsv_matrix(z_array[, , cell_type, drop = FALSE][, , 1], out_path, digits = digits)
  }
}

write_fraction_input <- function(frac, obj_dir, name = "InstaPrismfrac", digits = 6) {
  frac_dir <- file.path(obj_dir, "frac_input")
  dir.create(frac_dir, recursive = TRUE, showWarnings = FALSE)
  idx <- is.finite(frac)
  frac[idx] <- round(frac[idx], digits)
  write.table(
    frac,
    file = file.path(frac_dir, paste0(name, ".txt")),
    sep = "\t",
    quote = FALSE,
    col.names = NA
  )
}

write_method_fraction <- function(frac, paths, method, name = "InstaPrismfrac", digits = 6) {
  frac_dir <- file.path(paths$output_dir, method)
  dir.create(frac_dir, recursive = TRUE, showWarnings = FALSE)
  idx <- is.finite(frac)
  frac[idx] <- round(frac[idx], digits)
  write.table(
    frac,
    file = file.path(frac_dir, paste0(name, ".txt")),
    sep = "\t",
    quote = FALSE,
    col.names = NA
  )
}



############# Runtime logging #############

message <- function(..., domain = NULL, appendLF = TRUE) {
  base::message(..., domain = domain, appendLF = appendLF)
  log_path <- getOption("ctse.message_log_path", "")
  if (nzchar(log_path)) {
    text <- paste0(..., collapse = "")
    cat(text, if (appendLF) "
" else "", file = log_path, append = TRUE, sep = "")
  }
}

runtime_log_columns <- function() {
  c(
    "run_id", "dataset", "method", "deconv_config", "bulk_input", "frac_input",
    "refType", "start_time", "n_core", "deconv_input", "status", "run_time"
  )
}

safe_log_component <- function(x) {
  gsub("[^A-Za-z0-9_.-]+", "_", x)
}

method_message_log_path <- function(paths, method, timestamp = format(Sys.time(), "%Y%m%d_%H%M%S")) {
  file.path(
    paths$obj_dir,
    "logs",
    "run_logs",
    paste0(safe_log_component(method), "_", safe_log_component(paths$config_id), "_", timestamp, ".log")
  )
}

make_run_id <- function(method, config_id, timestamp = format(Sys.time(), "%Y%m%d_%H%M%S")) {
  paste0(safe_log_component(method), "_", safe_log_component(config_id), "_", timestamp)
}

run_id_from_log_path <- function(log_path, method, config_id) {
  if (!is.null(log_path) && nzchar(log_path)) {
    stem <- sub("\\.log$", "", basename(log_path))
    expected_prefix <- paste0(safe_log_component(method), "_", safe_log_component(config_id), "_")
    if (startsWith(stem, expected_prefix)) {
      return(stem)
    }
  }
  make_run_id(method, config_id)
}

runtime_methods_without_frac_input <- function() {
  c("BLUE", "scTAPE", "InstaPrism", "InstaPrismUpdated")
}

runtime_methods_without_n_core <- function() {
  c("BLUE", "scTAPE", "CIBERSORTx", "ENIGMAL2", "ENIGMAtrace")
}

runtime_frac_input <- function(paths, method) {
  if (method %in% runtime_methods_without_frac_input()) {
    return("NA")
  }
  paths$config$frac_input[[1]]
}

runtime_n_core <- function(method, n_core) {
  if (method %in% runtime_methods_without_n_core() || length(n_core) == 0 || is.null(n_core) || is.na(n_core)) {
    return("NA")
  }
  as.character(n_core)
}


quote_command_arg <- function(x) {
  if (grepl("^[A-Za-z0-9_./:=,+-]+$", x)) {
    x
  } else {
    shQuote(x)
  }
}

format_runner_command <- function() {
  full_args <- commandArgs(FALSE)
  script <- NA_character_
  file_arg <- full_args[grepl("^--file=", full_args)]
  if (length(file_arg) > 0) {
    script <- sub("^--file=", "", file_arg[[1]])
  } else {
    file_idx <- match("--file", full_args)
    if (!is.na(file_idx) && length(full_args) >= file_idx + 1) {
      script <- full_args[[file_idx + 1]]
    }
  }
  if (is.na(script) || !nzchar(script)) {
    script <- "<script>"
  }
  args <- c("Rscript", script, commandArgs(TRUE))
  paste(vapply(args, quote_command_arg, character(1)), collapse = " ")
}

start_message_log <- function(log_path) {
  dir.create(dirname(log_path), recursive = TRUE, showWarnings = FALSE)
  old_log_path <- getOption("ctse.message_log_path", "")
  out_start <- sink.number(type = "output")
  out_con <- file(log_path, open = "at")
  options(ctse.message_log_path = log_path)
  sink(out_con, split = TRUE)
  message("message_log: ", log_path)
  message("command: ", format_runner_command())
  closed <- FALSE
  function() {
    if (closed) {
      return(invisible(TRUE))
    }
    while (sink.number(type = "output") > out_start) {
      sink(type = "output")
    }
    close(out_con)
    options(ctse.message_log_path = old_log_path)
    closed <<- TRUE
    invisible(TRUE)
  }
}

with_message_log <- function(log_path, expr) {
  cleanup_message_log <- start_message_log(log_path)
  on.exit(cleanup_message_log(), add = TRUE)
  force(expr)
}

format_runtime_timestamp <- function(x = Sys.time()) {
  format(x, "%Y-%m-%d %H:%M:%OS3 %Z")
}

format_elapsed_mins <- function(start_clock, end_clock = Sys.time()) {
  elapsed <- end_clock - start_clock
  units(elapsed) <- "mins"
  paste0(round(as.numeric(elapsed), 2), "mins")
}

read_runtime_log <- function(paths) {
  cols <- runtime_log_columns()
  if (!file.exists(paths$runtime_path)) {
    out <- as.data.frame(setNames(rep(list(character()), length(cols)), cols), stringsAsFactors = FALSE)
    return(out)
  }

  out <- read.delim(paths$runtime_path, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE, na.strings = character())
  missing_cols <- setdiff(cols, colnames(out))
  for (col in missing_cols) {
    out[[col]] <- ""
  }
  out[, cols, drop = FALSE]
}

write_runtime_log_table <- function(paths, runtime_log) {
  dir.create(dirname(paths$runtime_path), recursive = TRUE, showWarnings = FALSE)
  write.table(
    runtime_log[, runtime_log_columns(), drop = FALSE],
    file = paths$runtime_path,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
}

runtime_log_row <- function(paths, run_id, method, start_time, n_core, deconv_input, status, run_time = "") {
  data.frame(
    run_id = run_id,
    dataset = paths$dataset,
    method = method,
    deconv_config = paths$config_id,
    bulk_input = paths$config$bulk_input[[1]],
    frac_input = runtime_frac_input(paths, method),
    refType = paths$refType,
    start_time = start_time,
    n_core = runtime_n_core(method, n_core),
    deconv_input = deconv_input,
    status = status,
    run_time = run_time,
    stringsAsFactors = FALSE
  )
}

start_runtime_log <- function(paths, method, n_core = NA_integer_, deconv_input = "", log_path = getOption("ctse.message_log_path", "")) {
  start_clock <- Sys.time()
  start_time <- format_runtime_timestamp(start_clock)
  run_id <- run_id_from_log_path(log_path, method, paths$config_id)
  runtime_log <- read_runtime_log(paths)
  runtime_log <- rbind(runtime_log, runtime_log_row(paths, run_id, method, start_time, n_core, deconv_input, "running"))
  write_runtime_log_table(paths, runtime_log)

  list(
    run_id = run_id,
    method = method,
    start_time = start_time,
    start_clock = start_clock,
    n_core = n_core,
    deconv_input = deconv_input
  )
}

finish_runtime_log <- function(paths, runtime_entry, status = "completed", run_time = NULL) {
  if (is.null(run_time)) {
    run_time <- format_elapsed_mins(runtime_entry$start_clock)
  }

  runtime_log <- read_runtime_log(paths)
  match_idx <- which(runtime_log$run_id == runtime_entry$run_id)

  if (length(match_idx) == 0) {
    deconv_input <- if (!is.null(runtime_entry$deconv_input)) runtime_entry$deconv_input else ""
    runtime_log <- rbind(runtime_log, runtime_log_row(paths, runtime_entry$run_id, runtime_entry$method, runtime_entry$start_time, runtime_entry$n_core, deconv_input, status, run_time))
  } else {
    match_idx <- match_idx[[length(match_idx)]]
    runtime_log$status[match_idx] <- status
    runtime_log$run_time[match_idx] <- run_time
  }

  write_runtime_log_table(paths, runtime_log)
  invisible(run_time)
}

with_started_runtime_log <- function(paths, runtime_entry, expr) {
  tryCatch(
    {
      value <- force(expr)
      run_time <- finish_runtime_log(paths, runtime_entry, status = "completed")
      list(value = value, run_time = run_time)
    },
    interrupt = function(e) {
      finish_runtime_log(paths, runtime_entry, status = "killed")
      stop(e)
    },
    error = function(e) {
      finish_runtime_log(paths, runtime_entry, status = "failed")
      stop(e)
    }
  )
}

with_runtime_log <- function(paths, method, n_core = NA_integer_, expr, deconv_input = "", log_path = getOption("ctse.message_log_path", "")) {
  runtime_entry <- start_runtime_log(paths, method, n_core = n_core, deconv_input = deconv_input, log_path = log_path)
  with_started_runtime_log(paths, runtime_entry, force(expr))
}



############# Extra args parsing / validation #############

parse_extra_args <- function(extra_args) {
  if (is.null(extra_args) || is.na(extra_args) || !nzchar(extra_args)) {
    out <- list()
    attr(out, "user_keys") <- character()
    return(out)
  }
  pieces <- unlist(strsplit(extra_args, ";", fixed = TRUE))
  pieces <- trimws(pieces)
  pieces <- pieces[nzchar(pieces)]
  out <- list()
  for (piece in pieces) {
    kv <- strsplit(piece, "=", fixed = TRUE)[[1]]
    if (length(kv) != 2) {
      stop("Invalid extra_args piece: ", piece, ". Expected key=value")
    }
    key <- trimws(kv[[1]])
    value <- trimws(kv[[2]])
    value_lower <- tolower(value)
    if (value_lower %in% c("true", "false")) {
      value <- value_lower == "true"
    } else if (grepl("^[0-9]+$", value)) {
      value <- as.integer(value)
    } else if (grepl("^[0-9.]+$", value)) {
      value <- as.numeric(value)
    }
    out[[key]] <- value
  }
  attr(out, "user_keys") <- names(out)
  out
}

parse_method_default_value <- function(value, value_type) {
  value_type <- tolower(value_type)
  if (value_type == "null") {
    return(NULL)
  }
  if (value_type == "logical") {
    value_lower <- tolower(value)
    if (!value_lower %in% c("true", "false")) {
      stop("Invalid logical default value: ", value)
    }
    return(value_lower == "true")
  }
  if (value_type == "integer") {
    return(as.integer(value))
  }
  if (value_type == "numeric") {
    return(as.numeric(value))
  }
  if (value_type == "character") {
    return(value)
  }
  stop("Unsupported method default value_type: ", value_type)
}

read_method_default_extra_configs <- function(repo_root = find_repo_root()) {
  path <- file.path(repo_root, "DALE_Eval", "configs", "method_default_extra_configs.txt")
  if (!file.exists(path)) {
    stop("Method default extra config file not found: ", path)
  }
  out <- read.delim(path, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE, na.strings = character())
  required <- c("method", "extra_arg", "default_value", "value_type")
  missing <- setdiff(required, colnames(out))
  if (length(missing) > 0) {
    stop("method_default_extra_configs.txt missing columns: ", paste(missing, collapse = ", "))
  }
  out
}

apply_method_default_extra_args <- function(extra, method, repo_root = find_repo_root()) {
  user_keys <- attr(extra, "user_keys")
  if (is.null(user_keys)) {
    user_keys <- names(extra)
  }
  defaults <- read_method_default_extra_configs(repo_root = repo_root)
  defaults <- defaults[defaults$method == method, , drop = FALSE]
  merged <- list()
  for (i in seq_len(nrow(defaults))) {
    value <- parse_method_default_value(defaults$default_value[[i]], defaults$value_type[[i]])
    if (!is.null(value)) {
      merged[[defaults$extra_arg[[i]]]] <- value
    }
  }
  for (key in names(extra)) {
    merged[[key]] <- extra[[key]]
  }
  attr(merged, "user_keys") <- user_keys
  merged
}

validate_extra_args <- function(extra, allowed, method) {
  unknown <- setdiff(names(extra), allowed)
  if (length(unknown) > 0) {
    allowed_msg <- if (length(allowed) > 0) paste(allowed, collapse = ", ") else "none"
    stop(
      method, " received unsupported extra_args: ", paste(unknown, collapse = ", "),
      ". Allowed extra_args: ", allowed_msg
    )
  }
  invisible(TRUE)
}

validate_top_n_arg <- function(extra, use_limma_top_genes, method) {
  user_keys <- attr(extra, "user_keys")
  if (is.null(user_keys)) {
    user_keys <- names(extra)
  }
  if ("top_n" %in% user_keys && !isTRUE(use_limma_top_genes)) {
    stop(
      method, " received top_n, but use_limma_top_genes is false. ",
      "Set use_limma_top_genes=true to run on top marker genes."
    )
  }
  invisible(TRUE)
}

extra_arg <- function(extra, key, default = NULL) {
  if (!is.null(extra[[key]])) extra[[key]] else default
}



############# File lookup helpers #############

existing_file <- function(paths, label = "file") {
  hit <- paths[file.exists(paths)]
  if (length(hit) == 0) {
    stop(label, " not found. Checked: ", paste(paths, collapse = ", "))
  }
  hit[[1]]
}

ref_file <- function(ref_dir, basenames, label = "reference file") {
  existing_file(file.path(ref_dir, basenames), label)
}
