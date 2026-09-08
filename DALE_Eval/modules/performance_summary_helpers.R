# Helpers for summarizing CTSE correlation performance files.

read_summary_numeric_matrix <- function(path, required = TRUE) {
  if (!file.exists(path)) {
    if (required) {
      stop("Missing matrix: ", path)
    }
    return(NULL)
  }

  x <- read.delim(path, sep = "\t", check.names = FALSE, row.names = 1)
  x <- as.matrix(x)
  storage.mode(x) <- "numeric"
  x
}

read_summary_csv_matrix <- function(path) {
  if (!file.exists(path)) {
    stop("Missing CSV matrix: ", path)
  }

  x <- read.delim(path, sep = ",", check.names = FALSE, row.names = 1)
  x <- as.matrix(x)
  storage.mode(x) <- "numeric"
  x
}

read_benchmark_dataset_info <- function(repo_root = find_repo_root()) {
  path <- file.path(repo_root, "DALE_Eval", "configs", "benchmark_dataset_info.txt")
  if (!file.exists(path)) {
    stop("benchmark_dataset_info.txt not found: ", path)
  }

  info <- read.delim(path, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
  required_cols <- c(
    "dataset", "cell_type_raw", "cell_type", "ct_abundance",
    "transcript_abundance", "sample_size", "tissue"
  )
  missing_cols <- setdiff(required_cols, colnames(info))
  if (length(missing_cols) > 0L) {
    stop("benchmark_dataset_info.txt missing columns: ", paste(missing_cols, collapse = ", "))
  }

  info
}

standardize_summary_cell_types <- function(dataset, cell_types_raw, dataset_info) {
  info <- dataset_info[dataset_info$dataset == dataset, , drop = FALSE]
  map <- setNames(info$cell_type, info$cell_type_raw)
  out <- ifelse(cell_types_raw %in% names(map), map[cell_types_raw], cell_types_raw)
  unname(out)
}

read_bulk_gene_names <- function(dataset, bulk_input, repo_root = find_repo_root()) {
  path <- file.path(repo_root, "Benchmarking_obj", dataset, "bulk_input", paste0(bulk_input, ".txt"))
  if (!file.exists(path)) {
    stop("Bulk input not found: ", path)
  }

  rownames(read.delim(path, sep = "\t", check.names = FALSE, row.names = 1))
}

read_cor_na_acceptable_mask <- function(dataset, truth_type, repo_root = find_repo_root()) {
  path <- file.path(
    repo_root,
    "Benchmarking_obj",
    dataset,
    "evaluation_metadata",
    "cor_na_acceptable_mask",
    paste0("cor_na_acceptable_mask_truth-", truth_type, ".txt")
  )
  read_summary_numeric_matrix(path, required = TRUE)
}

read_limma_top_genes_csv <- function(ref_dir) {
  path <- file.path(ref_dir, "limma_top_genes.csv")
  read_summary_csv_matrix(path)
}

rank_limma_genes_for_cell_type <- function(limma_stats, cell_type) {
  if (!cell_type %in% colnames(limma_stats)) {
    return(character(0))
  }

  stats <- limma_stats[, cell_type]
  keep <- is.finite(stats)
  if (!any(keep)) {
    return(character(0))
  }

  rownames(limma_stats)[keep][order(stats[keep], decreasing = TRUE)]
}

read_summary_indep_mapping <- function(dataset, indep_ref, repo_root = find_repo_root()) {
  path <- file.path(repo_root, "Indep_scReference", "cell_type_mapping.txt")
  if (!file.exists(path)) {
    stop("Independent-reference cell type mapping not found: ", path)
  }

  mapping <- read.delim(path, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
  required_cols <- c("dataset", "indep_ref", "target_cell_type", "indep_ref_cell_type")
  missing_cols <- setdiff(required_cols, colnames(mapping))
  if (length(missing_cols) > 0L) {
    stop("cell_type_mapping.txt missing columns: ", paste(missing_cols, collapse = ", "))
  }

  mapping[mapping$dataset == dataset & mapping$indep_ref == indep_ref, , drop = FALSE]
}

resolve_limma_cell_type <- function(cell_type_raw, ref_type, limma_stats, indep_mapping = NULL) {
  if (ref_type == "self") {
    return(cell_type_raw)
  }

  if (!is.null(indep_mapping) && cell_type_raw %in% indep_mapping$target_cell_type) {
    return(indep_mapping$indep_ref_cell_type[match(cell_type_raw, indep_mapping$target_cell_type)])
  }

  # Allow unmapped independent-reference output columns to use their original names.
  if (cell_type_raw %in% colnames(limma_stats)) {
    return(cell_type_raw)
  }

  NA_character_
}

build_top_gene_groups_by_cell_type <- function(dataset,
                                               config_row,
                                               cell_types_raw,
                                               top_n_values,
                                               repo_root = find_repo_root()) {
  ref_type <- config_row$refType[[1]]
  bulk_input <- config_row$bulk_input[[1]]
  bulk_genes <- read_bulk_gene_names(dataset, bulk_input, repo_root = repo_root)

  if (ref_type == "self") {
    ref_dir <- file.path(repo_root, "Benchmarking_obj", dataset, "self_reference")
    indep_mapping <- NULL
  } else if (ref_type == "indep") {
    ref_dir <- resolve_indep_ref_dir(dataset, repo_root = repo_root)
    indep_ref <- basename(ref_dir)
    indep_mapping <- read_summary_indep_mapping(dataset, indep_ref, repo_root = repo_root)
  } else {
    stop("Unsupported refType: ", ref_type)
  }

  limma_stats <- read_limma_top_genes_csv(ref_dir)
  groups_by_cell_type <- list()

  for (cell_type_raw in cell_types_raw) {
    limma_cell_type <- resolve_limma_cell_type(
      cell_type_raw = cell_type_raw,
      ref_type = ref_type,
      limma_stats = limma_stats,
      indep_mapping = indep_mapping
    )

    ranked_genes <- if (is.na(limma_cell_type)) {
      character(0)
    } else {
      rank_limma_genes_for_cell_type(limma_stats, limma_cell_type)
    }

    if (ref_type == "indep") {
      ranked_genes <- ranked_genes[ranked_genes %in% bulk_genes]
    }

    groups <- setNames(vector("list", length(top_n_values)), paste0("top_", top_n_values))
    for (i in seq_along(top_n_values)) {
      n <- top_n_values[[i]]
      groups[[i]] <- ranked_genes[seq_len(min(n, length(ranked_genes)))]
    }
    groups[["all_genes"]] <- NULL
    groups_by_cell_type[[cell_type_raw]] <- groups
  }

  groups_by_cell_type
}

summary_mean <- function(x) {
  if (length(x) == 0L || all(is.na(x))) {
    return(NA_real_)
  }
  mean(x, na.rm = TRUE)
}

summarize_cor_vector_by_groups <- function(cor_values,
                                           acceptable_mask_values,
                                           gene_groups) {
  out <- vector("list", length(gene_groups))
  names(out) <- names(gene_groups)

  for (group_name in names(gene_groups)) {
    genes <- gene_groups[[group_name]]
    genes <- intersect(genes, names(cor_values))

    if (length(genes) == 0L) {
      out[[group_name]] <- data.frame(
        group = group_name,
        avg_cor = NA_real_,
        avg_cor_with_NA_penalty = NA_real_,
        n_genes = 0L,
        n_constant_genes = 0L,
        stringsAsFactors = FALSE
      )
      next
    }

    values <- cor_values[genes]
    acceptable <- acceptable_mask_values[genes]
    acceptable[is.na(acceptable)] <- FALSE

    keep <- !acceptable
    values_keep <- values[keep]
    n_constant_genes <- sum(is.na(values_keep))
    values_penalized <- values_keep
    values_penalized[is.na(values_penalized)] <- 0

    out[[group_name]] <- data.frame(
      group = group_name,
      avg_cor = summary_mean(values_keep),
      avg_cor_with_NA_penalty = if (length(values_penalized) == 0L) NA_real_ else mean(values_penalized),
      n_genes = length(values_keep),
      n_constant_genes = n_constant_genes,
      stringsAsFactors = FALSE
    )
  }

  do.call(rbind, out)
}


summarize_baseline_cor_vector_by_groups <- function(cor_values,
                                                    acceptable_mask_values,
                                                    gene_groups,
                                                    value_col) {
  out <- vector("list", length(gene_groups))
  names(out) <- names(gene_groups)

  for (group_name in names(gene_groups)) {
    genes <- gene_groups[[group_name]]
    genes <- intersect(genes, names(cor_values))

    if (length(genes) == 0L) {
      out[[group_name]] <- data.frame(
        group = group_name,
        value = NA_real_,
        stringsAsFactors = FALSE
      )
      next
    }

    values <- cor_values[genes]
    acceptable <- acceptable_mask_values[genes]
    acceptable[is.na(acceptable)] <- FALSE
    values_keep <- values[!acceptable]

    out[[group_name]] <- data.frame(
      group = group_name,
      value = summary_mean(values_keep),
      stringsAsFactors = FALSE
    )
  }

  out <- do.call(rbind, out)
  colnames(out)[colnames(out) == "value"] <- value_col
  out
}

baseline_specs_for_config <- function(config_row, truth_type) {
  value_cols <- c(
    "avg_cor_bulk",
    "avg_cor_bulk_truthFrac_regressed",
    "avg_cor_bulk_InstaPrismFrac_regressed"
  )

  empty_specs <- data.frame(
    value_col = value_cols,
    baseline_file = NA_character_,
    stringsAsFactors = FALSE
  )

  config_id <- config_row$config_id[[1]]
  bulk_input <- config_row$bulk_input[[1]]
  bulk_normalization <- tolower(config_row$bulk_normalization[[1]])

  instaprism_config_id <- switch(
    config_id,
    config03 = "config03",
    config04 = "config04",
    config05 = "config05",
    config06 = "config06",
    config02 = "config02",
    config08 = "config08",
    config01 = "config01",
    config07 = "config07",
    config09 = "config01",
    NA_character_
  )

  if (is.na(instaprism_config_id)) {
    return(empty_specs)
  }

  if (truth_type == "meancpm") {
    baseline_id <- paste0("bulk-", bulk_input)
    truth_regress_frac <- "truth_cellfrac"
  } else if (truth_type == "sumcount_cpm" && bulk_normalization == "cpm") {
    baseline_id <- paste0("bulk-", bulk_input, "__norm-cpm")
    truth_regress_frac <- "truth_transcriptfrac"
  } else {
    return(empty_specs)
  }

  data.frame(
    value_col = value_cols,
    baseline_file = c(
      paste0("truth-", truth_type, "__", baseline_id, ".txt"),
      paste0("truth-", truth_type, "__", baseline_id, "__regress-", truth_regress_frac, ".txt"),
      paste0("truth-", truth_type, "__", baseline_id, "__regress-InstaPrismfrac_", instaprism_config_id, ".txt")
    ),
    stringsAsFactors = FALSE
  )
}

summarize_baseline_file <- function(dataset,
                                    baseline_path,
                                    value_col,
                                    mask,
                                    dataset_info,
                                    top_groups_by_cell_type,
                                    top_n_values) {
  cor_mat <- read_summary_numeric_matrix(baseline_path, required = TRUE)
  cell_types_raw <- colnames(cor_mat)
  cell_types <- standardize_summary_cell_types(dataset, cell_types_raw, dataset_info)

  rows <- list()

  for (i in seq_along(cell_types_raw)) {
    cell_type_raw <- cell_types_raw[[i]]
    cell_type <- cell_types[[i]]

    cor_values <- cor_mat[, cell_type_raw]
    names(cor_values) <- rownames(cor_mat)

    acceptable_mask_values <- rep(FALSE, length(cor_values))
    names(acceptable_mask_values) <- names(cor_values)
    if (cell_type_raw %in% colnames(mask)) {
      common_mask_genes <- intersect(names(cor_values), rownames(mask))
      acceptable_mask_values[common_mask_genes] <- mask[common_mask_genes, cell_type_raw] == 1
    }

    top_groups <- top_groups_by_cell_type[[cell_type_raw]]
    if (is.null(top_groups)) {
      empty_groups <- setNames(vector("list", length(top_n_values)), paste0("top_", top_n_values))
      top_groups <- empty_groups
    }
    gene_groups <- build_gene_groups_for_summary(rownames(cor_mat), top_groups)

    summary_df <- summarize_baseline_cor_vector_by_groups(
      cor_values = cor_values,
      acceptable_mask_values = acceptable_mask_values,
      gene_groups = gene_groups,
      value_col = value_col
    )
    summary_df$cell_type <- cell_type
    rows[[length(rows) + 1L]] <- summary_df[, c("cell_type", "group", value_col)]
  }

  do.call(rbind, rows)
}

summarize_baselines_for_config <- function(dataset,
                                           config_row,
                                           truth_type,
                                           mask,
                                           dataset_info,
                                           top_groups_by_cell_type,
                                           top_n_values,
                                           repo_root = find_repo_root()) {
  specs <- baseline_specs_for_config(config_row, truth_type)
  baseline_root <- file.path(
    repo_root,
    "Benchmarking_obj",
    dataset,
    "deconv_performance",
    "bulk_baseline",
    "spearman_cor"
  )

  baseline_summary <- NULL

  for (i in seq_len(nrow(specs))) {
    value_col <- specs$value_col[[i]]
    baseline_file <- specs$baseline_file[[i]]

    if (is.na(baseline_file)) {
      next
    }

    baseline_path <- file.path(baseline_root, baseline_file)
    if (!file.exists(baseline_path)) {
      next
    }

    table_one <- summarize_baseline_file(
      dataset = dataset,
      baseline_path = baseline_path,
      value_col = value_col,
      mask = mask,
      dataset_info = dataset_info,
      top_groups_by_cell_type = top_groups_by_cell_type,
      top_n_values = top_n_values
    )

    if (is.null(baseline_summary)) {
      baseline_summary <- table_one
    } else {
      baseline_summary <- merge(
        baseline_summary,
        table_one,
        by = c("cell_type", "group"),
        all = TRUE,
        sort = FALSE
      )
    }
  }

  for (value_col in specs$value_col) {
    if (is.null(baseline_summary) || !value_col %in% colnames(baseline_summary)) {
      if (is.null(baseline_summary)) {
        baseline_summary <- data.frame(cell_type = character(0), group = character(0), stringsAsFactors = FALSE)
      }
      baseline_summary[[value_col]] <- numeric(0)
    }
  }

  baseline_summary[, c("cell_type", "group", specs$value_col), drop = FALSE]
}

build_gene_groups_for_summary <- function(cor_genes, top_groups) {
  groups <- top_groups
  groups[["all_genes"]] <- cor_genes
  groups
}

summarize_spearman_file <- function(dataset,
                                    method,
                                    cor_path,
                                    mask,
                                    config_id,
                                    config_row,
                                    truth_type,
                                    dataset_info,
                                    top_groups_by_cell_type,
                                    top_n_values) {
  cor_mat <- read_summary_numeric_matrix(cor_path, required = TRUE)
  cell_types_raw <- colnames(cor_mat)
  cell_types <- standardize_summary_cell_types(dataset, cell_types_raw, dataset_info)

  rows <- list()

  for (i in seq_along(cell_types_raw)) {
    cell_type_raw <- cell_types_raw[[i]]
    cell_type <- cell_types[[i]]

    cor_values <- cor_mat[, cell_type_raw]
    names(cor_values) <- rownames(cor_mat)

    acceptable_mask_values <- rep(FALSE, length(cor_values))
    names(acceptable_mask_values) <- names(cor_values)
    if (cell_type_raw %in% colnames(mask)) {
      common_mask_genes <- intersect(names(cor_values), rownames(mask))
      acceptable_mask_values[common_mask_genes] <- mask[common_mask_genes, cell_type_raw] == 1
    }

    top_groups <- top_groups_by_cell_type[[cell_type_raw]]
    if (is.null(top_groups)) {
      empty_groups <- setNames(vector("list", length(top_n_values)), paste0("top_", top_n_values))
      top_groups <- empty_groups
    }
    gene_groups <- build_gene_groups_for_summary(rownames(cor_mat), top_groups)

    summary_df <- summarize_cor_vector_by_groups(
      cor_values = cor_values,
      acceptable_mask_values = acceptable_mask_values,
      gene_groups = gene_groups
    )

    summary_df$method <- method
    summary_df$cell_type <- cell_type
    summary_df$refType <- config_row$refType[[1]]
    summary_df$bulk_input <- config_row$bulk_input[[1]]
    summary_df$bulk_normalization <- config_row$bulk_normalization[[1]]
    summary_df$frac_input <- config_row$frac_input[[1]]
    summary_df$config_id <- config_id
    summary_df$truth_type <- truth_type

    rows[[length(rows) + 1L]] <- summary_df[, c(
      "method",
      "cell_type",
      "group",
      "avg_cor",
      "avg_cor_with_NA_penalty",
      "n_genes",
      "n_constant_genes",
      "refType",
      "bulk_input",
      "bulk_normalization",
      "frac_input",
      "config_id",
      "truth_type"
    )]
  }

  do.call(rbind, rows)
}

find_performance_dirs_for_config <- function(dataset,
                                             config_row,
                                             truth_type_include = NULL,
                                             repo_root = find_repo_root()) {
  slug <- config_slug(config_row)
  perf_root <- file.path(repo_root, "Benchmarking_obj", dataset, "deconv_performance")
  if (!dir.exists(perf_root)) {
    return(data.frame())
  }

  dirs <- list.dirs(perf_root, full.names = TRUE, recursive = FALSE)
  dirs <- dirs[grepl(paste0("^", slug, "__truth-"), basename(dirs))]
  if (length(dirs) == 0L) {
    return(data.frame())
  }

  truth_type <- sub("^.*__truth-", "", basename(dirs))
  out <- data.frame(
    performance_dir = dirs,
    truth_type = truth_type,
    stringsAsFactors = FALSE
  )

  if (!is.null(truth_type_include)) {
    out <- out[out$truth_type %in% truth_type_include, , drop = FALSE]
  }

  out[order(out$truth_type), , drop = FALSE]
}
