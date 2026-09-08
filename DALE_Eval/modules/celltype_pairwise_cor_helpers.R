# Helpers for within-result cell-type pairwise Pearson correlations.

pairwise_cor_separator <- function() {
  "--"
}

validate_pairwise_cell_types <- function(cell_types,
                                         separator = pairwise_cor_separator(),
                                         label = "cell types") {
  cell_types <- as.character(cell_types)
  if (length(cell_types) < 2L) {
    stop(label, " must contain at least two names")
  }
  if (anyNA(cell_types) || any(!nzchar(cell_types))) {
    stop(label, " contain missing or blank names")
  }
  if (any(grepl("[[:space:]]", cell_types))) {
    stop(label, " contain whitespace, which is not supported in pair names: ",
         paste(cell_types[grepl("[[:space:]]", cell_types)], collapse = ", "))
  }
  if (any(grepl(separator, cell_types, fixed = TRUE))) {
    stop(label, " contain the reserved pair separator '", separator, "': ",
         paste(cell_types[grepl(separator, cell_types, fixed = TRUE)], collapse = ", "))
  }
  duplicated_names <- unique(cell_types[duplicated(cell_types)])
  if (length(duplicated_names) > 0L) {
    stop(label, " contain duplicates: ", paste(duplicated_names, collapse = ", "))
  }
  invisible(cell_types)
}

make_celltype_pair_table <- function(cell_types,
                                     separator = pairwise_cor_separator()) {
  validate_pairwise_cell_types(cell_types, separator = separator)
  pair_matrix <- utils::combn(cell_types, 2L)
  data.frame(
    pair = paste(pair_matrix[1L, ], pair_matrix[2L, ], sep = separator),
    cell_type_a = pair_matrix[1L, ],
    cell_type_b = pair_matrix[2L, ],
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

parse_celltype_pair_names <- function(pair_names,
                                      separator = pairwise_cor_separator()) {
  pair_names <- as.character(pair_names)
  if (length(pair_names) == 0L) {
    return(data.frame(
      pair = character(),
      cell_type_a = character(),
      cell_type_b = character(),
      stringsAsFactors = FALSE
    ))
  }
  if (anyNA(pair_names) || any(!nzchar(pair_names))) {
    stop("Pair names contain missing or blank values")
  }
  pieces <- strsplit(pair_names, separator, fixed = TRUE)
  valid <- lengths(pieces) == 2L & vapply(
    pieces,
    function(x) all(nzchar(x)),
    logical(1)
  )
  if (any(!valid)) {
    stop(
      "Invalid cell-type pair names; expected exactly '<cell_type_a>",
      separator,
      "<cell_type_b>': ",
      paste(pair_names[!valid], collapse = ", ")
    )
  }
  if (anyDuplicated(pair_names)) {
    stop("Duplicated cell-type pair names: ",
         paste(unique(pair_names[duplicated(pair_names)]), collapse = ", "))
  }
  data.frame(
    pair = pair_names,
    cell_type_a = vapply(pieces, `[[`, character(1), 1L),
    cell_type_b = vapply(pieces, `[[`, character(1), 2L),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

celltype_pair_columns <- function(pair_names,
                                  cell_type,
                                  separator = pairwise_cor_separator()) {
  if (length(cell_type) != 1L || is.na(cell_type) || !nzchar(cell_type)) {
    stop("cell_type must be one non-empty name")
  }
  parsed <- parse_celltype_pair_names(pair_names, separator = separator)
  parsed$pair[parsed$cell_type_a == cell_type | parsed$cell_type_b == cell_type]
}

read_pairwise_ctse_matrix <- function(path, label = "CTSE matrix") {
  if (!file.exists(path)) {
    stop(label, " not found: ", path)
  }
  x <- read.delim(path, sep = "\t", check.names = FALSE, row.names = 1)
  x <- as.matrix(x)
  suppressWarnings(storage.mode(x) <- "double")
  if (is.null(rownames(x)) || is.null(colnames(x))) {
    stop(label, " must have gene row names and sample column names: ", path)
  }
  if (anyDuplicated(rownames(x)) || anyDuplicated(colnames(x))) {
    stop(label, " has duplicated gene or sample names: ", path)
  }
  if (any(is.infinite(x), na.rm = TRUE)) {
    stop(label, " contains infinite values: ", path)
  }
  x
}

read_pairwise_test_samples <- function(obj_dir) {
  path <- file.path(obj_dir, "self_reference", "sample_split.txt")
  if (!file.exists(path)) {
    stop("Testing-sample split not found: ", path)
  }
  split <- read.delim(
    path,
    sep = "\t",
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  required <- c("group", "sampleIDs")
  missing <- setdiff(required, colnames(split))
  if (length(missing) > 0L) {
    stop("sample_split.txt missing columns: ", paste(missing, collapse = ", "))
  }
  group <- tolower(trimws(as.character(split$group)))
  samples <- as.character(split$sampleIDs[group == "test"])
  samples <- unique(samples[!is.na(samples) & nzchar(samples)])
  if (length(samples) == 0L) {
    stop("No testing samples found in: ", path)
  }
  list(path = normalizePath(path), samples = samples)
}

read_pairwise_eval_config <- function(config_id, repo_root = find_repo_root()) {
  path <- file.path(repo_root, "DALE_Eval", "configs", "eval_configs.txt")
  if (!file.exists(path)) {
    stop("Evaluation config table not found: ", path)
  }
  configs <- read.delim(path, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
  required <- c("config_id", "expected_truth_type")
  missing <- setdiff(required, colnames(configs))
  if (length(missing) > 0L) {
    stop("eval_configs.txt missing columns: ", paste(missing, collapse = ", "))
  }
  row <- configs[configs$config_id == config_id, , drop = FALSE]
  if (nrow(row) != 1L) {
    supported <- paste(configs$config_id, collapse = ", ")
    stop(
      "Unsupported evaluation config_id=", config_id,
      ". Supported configs from eval_configs.txt: ", supported
    )
  }
  row
}

resolve_pairwise_method_cell_types <- function(paths,
                                               method_files,
                                               repo_root = find_repo_root()) {
  available <- names(method_files)
  if (paths$refType == "self") {
    selected <- sort(available)
    mapping <- NULL
  } else if (paths$refType == "indep") {
    mapping <- read_indep_ref_cell_type_mapping(paths, repo_root)
    if (nrow(mapping) == 0L) {
      stop(
        "No independent-reference mapping found for dataset=", paths$dataset,
        " and reference=", basename(paths$ref_dir)
      )
    }
    selected <- unique(as.character(mapping$target_cell_type))
    missing <- setdiff(selected, available)
    if (length(missing) > 0L) {
      stop(
        "Method output is missing mapped benchmark cell types: ",
        paste(missing, collapse = ", ")
      )
    }
  } else {
    stop("Unsupported refType: ", paths$refType)
  }
  validate_pairwise_cell_types(selected, label = "selected method cell types")
  list(
    selected = selected,
    extra = setdiff(available, selected),
    mapping = mapping
  )
}

compute_celltype_pairwise_pearson <- function(files,
                                              cell_types,
                                              test_samples,
                                              label = "CTSE") {
  validate_pairwise_cell_types(cell_types)
  missing_files <- setdiff(cell_types, names(files))
  if (length(missing_files) > 0L) {
    stop(label, " is missing cell-type files: ", paste(missing_files, collapse = ", "))
  }

  matrices <- setNames(vector("list", length(cell_types)), cell_types)
  for (cell_type in cell_types) {
    matrices[[cell_type]] <- read_pairwise_ctse_matrix(
      files[[cell_type]],
      paste0(label, " ", cell_type)
    )
  }

  common_genes <- sort(Reduce(intersect, lapply(matrices, rownames)))
  common_samples <- Reduce(intersect, lapply(matrices, colnames))
  common_samples <- test_samples[test_samples %in% common_samples]
  if (length(common_genes) == 0L) {
    stop(label, " has no genes shared across selected cell types")
  }
  if (length(common_samples) < 2L) {
    stop(label, " has fewer than two shared testing samples")
  }

  matrices <- lapply(
    matrices,
    function(x) x[common_genes, common_samples, drop = FALSE]
  )
  pairs <- make_celltype_pair_table(cell_types)
  values <- vector("list", nrow(pairs))
  for (i in seq_len(nrow(pairs))) {
    values[[i]] <- pair_cor_pearson_fast(
      matrices[[pairs$cell_type_a[[i]]]],
      matrices[[pairs$cell_type_b[[i]]]],
      margin = "row",
      pairwise_na = TRUE
    )
  }
  result <- do.call(cbind, values)
  rownames(result) <- common_genes
  colnames(result) <- pairs$pair

  list(
    matrix = result,
    pairs = pairs,
    n_genes = length(common_genes),
    n_samples = length(common_samples),
    samples = common_samples
  )
}

write_pairwise_cor_matrix <- function(x, path, digits = 6L) {
  if (length(digits) != 1L || is.na(digits) || digits < 0L) {
    stop("digits must be one non-negative integer")
  }
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  if (file.exists(path)) {
    warning("Overwriting existing pairwise-correlation file: ", path)
  }
  x <- as.matrix(x)
  storage.mode(x) <- "double"
  finite <- is.finite(x)
  x[finite] <- round(x[finite], digits = digits)
  write.table(
    x,
    file = path,
    sep = "\t",
    quote = FALSE,
    col.names = NA
  )
}

pairwise_truth_output_label <- function(truth_type) {
  labels <- c(
    meancpm = "truth_meancpm",
    sumcount_cpm = "truth_sumcount_cpm"
  )
  if (!truth_type %in% names(labels)) {
    stop(
      "Unsupported truth_type=", truth_type,
      ". Supported truth types: ", paste(names(labels), collapse = ", ")
    )
  }
  unname(labels[[truth_type]])
}
