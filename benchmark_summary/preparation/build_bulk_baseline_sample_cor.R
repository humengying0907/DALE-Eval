# Within-sample Spearman correlation for bulk and InstaPrism-regressed bulk.
#
# Run from the repository root, scripts/, or benchmark_summary/preparation/.
# First run build_bulk_baseline_cor.R to prepare the regressed bulk matrices.
# Edit the settings below to choose datasets, a configuration, and gene scopes.
# Each output is a sample-by-cell-type matrix; correlations are across genes.
# Only samples marked test in self_reference/sample_split.txt are evaluated.
# Samples with truth fractions at or below the evaluation threshold retain NA.
# Running this script replaces the selected baseline tables and their provenance.

## settings ----

dataset_include <- NULL
config_id <- "config01"
gene_scopes <- c("all_genes", "hallmark_genes")
digits <- 6L
hallmark_gmt <- "other_source_data/h.all.v7.5.1.symbols.gmt"


## paths and configuration ----

repo_root_candidates <- unique(normalizePath(
  c(getwd(), file.path(getwd(), ".."), file.path(getwd(), "../..")),
  mustWork = FALSE
))
repo_root <- repo_root_candidates[file.exists(file.path(
  repo_root_candidates, "DALE_Eval", "configs", "deconv_configs.txt"
))][1L]
if (is.na(repo_root)) {
  stop("Run from the repository root, scripts/, or benchmark_summary/preparation/")
}

source(file.path(repo_root, "DALE_Eval", "modules", "config_helpers.R"))
source(file.path(repo_root, "DALE_Eval", "modules", "runner_helpers.R"))
source(file.path(repo_root, "DALE_Eval", "modules", "evalu.R"))

eval_configs <- read.delim(
  file.path(repo_root, "DALE_Eval", "configs", "eval_configs.txt"),
  stringsAsFactors = FALSE,
  check.names = FALSE
)
eval_contract <- eval_configs[eval_configs$config_id == config_id, , drop = FALSE]
if (nrow(eval_contract) != 1L) {
  stop("Expected one evaluation contract for config_id=", config_id)
}
truth_type <- eval_contract$expected_truth_type[[1L]]
filter_frac_name <- eval_contract$filter_frac[[1L]]
min_frac <- eval_contract$min_frac[[1L]]
baseline_ids <- c(
  bulk = eval_contract$baseline_bulk_id[[1L]],
  bulk_regressed = eval_contract$baseline_InstaPrismFrac_regressed_id[[1L]]
)
if (anyNA(baseline_ids) || any(!nzchar(baseline_ids))) {
  stop("Bulk baselines are not defined for config_id=", config_id)
}
if (length(gene_scopes) == 0L || anyDuplicated(gene_scopes) ||
    !all(gene_scopes %in% c("all_genes", "hallmark_genes"))) {
  stop("gene_scopes must select all_genes and/or hallmark_genes without duplicates")
}

hallmark_genes <- if ("hallmark_genes" %in% gene_scopes) {
  read_gmt_gene_union(file.path(repo_root, hallmark_gmt))
} else {
  NULL
}

dataset_info <- read.delim(
  file.path(repo_root, "DALE_Eval", "configs", "benchmark_dataset_info.txt"),
  stringsAsFactors = FALSE
)
datasets <- unique(dataset_info$dataset)
if (!is.null(dataset_include)) {
  if (!all(dataset_include %in% datasets)) {
    stop("dataset_include contains an unknown benchmark dataset")
  }
  datasets <- intersect(datasets, dataset_include)
}


## calculate and export ----

for (dataset_name in datasets) {
  message("Dataset: ", dataset_name, "; config: ", config_id)
  paths <- resolve_deconv_paths(dataset_name, config_id, repo_root = repo_root)
  truth_files <- ctse_cell_type_files(file.path(paths$obj_dir, "ctse_truth", truth_type))
  filter_frac <- as.matrix(read.delim(
    file.path(paths$obj_dir, "frac_input", paste0(filter_frac_name, ".txt")),
    row.names = 1,
    check.names = FALSE
  ))
  cell_types <- sort(intersect(names(truth_files), colnames(filter_frac)))
  if (length(cell_types) == 0L) {
    stop("No overlapping truth/fraction cell types for ", dataset_name)
  }

  output_root <- file.path(paths$obj_dir, "deconv_performance", "bulk_baseline")
  manifest_rows <- list()

  for (baseline_name in names(baseline_ids)) {
    baseline_id <- unname(baseline_ids[[baseline_name]])
    source_file <- if (baseline_name == "bulk") {
      file.path("Benchmarking_obj", dataset_name, "bulk_input",
                paste0(paths$config$bulk_input[[1L]], ".txt"))
    } else {
      file.path("Benchmarking_obj", dataset_name, "regressed_bulk",
                paste0(baseline_id, ".txt"))
    }
    if (baseline_name == "bulk") {
      bulk_prep <- prepare_bulk_for_deconv(read_bulk(paths), paths, "Bulk baseline")
      estimate <- bulk_prep$bulk_expr
    } else {
      # Regression was fitted by build_bulk_baseline_cor.R on all shared samples,
      # without an intercept, with the original gene means added to residuals.
      estimate <- read_ctse_matrix(file.path(repo_root, source_file), "Regressed bulk")
    }

    # Apply the same declared test split to both raw and regressed baselines.
    # Subsetting an existing regressed matrix does not refit its regression.
    estimate <- restrict_to_test_samples(
      bulk = estimate, paths = paths, frac = filter_frac, use_test_samples = TRUE
    )$bulk

    sample_vectors <- setNames(lapply(gene_scopes, function(x) list()), gene_scopes)
    for (cell_type in cell_types) {
      truth <- read_ctse_matrix(truth_files[[cell_type]], paste0("Truth ", cell_type))
      samples <- intersect(colnames(truth), colnames(estimate))
      samples <- intersect(samples, rownames(filter_frac))
      fractions <- filter_frac[samples, cell_type]
      eligible <- is.finite(fractions) & fractions > min_frac

      for (gene_scope in gene_scopes) {
        genes <- if (gene_scope == "hallmark_genes") hallmark_genes else NULL
        values <- compute_ctse_sample_cor(truth, estimate, genes = genes, samples = samples)
        values[!eligible] <- NA_real_
        sample_vectors[[gene_scope]][[cell_type]] <- values
      }
    }

    for (gene_scope in gene_scopes) {
      sample_cor <- merge_sample_metric_vectors(sample_vectors[[gene_scope]])
      metric_file <- file.path(
        paste0("sample_cor_", gene_scope),
        paste0("truth-", truth_type, "__", baseline_id, ".txt")
      )
      write_metric_matrix(sample_cor, file.path(output_root, metric_file), digits = digits)
      manifest_rows[[length(manifest_rows) + 1L]] <- data.frame(
        dataset = dataset_name,
        config_id = config_id,
        truth_type = truth_type,
        baseline_id = baseline_id,
        gene_scope = gene_scope,
        filter_frac = filter_frac_name,
        min_frac = min_frac,
        result_source = "current_inputs",
        source_file = source_file,
        metric_file = metric_file,
        n_samples = nrow(sample_cor),
        n_cell_types = ncol(sample_cor),
        digits = digits,
        sample_policy = "test_samples",
        stringsAsFactors = FALSE
      )
      message("  wrote ", metric_file, ": ", nrow(sample_cor),
              " samples x ", ncol(sample_cor), " cell types")
    }
  }

  # Replace provenance only for the tables generated in this run.
  manifest <- do.call(rbind, manifest_rows)
  manifest_path <- file.path(output_root, "sample_cor_manifest.tsv")
  if (file.exists(manifest_path)) {
    previous_manifest <- read.delim(
      manifest_path, stringsAsFactors = FALSE, check.names = FALSE
    )
    if (!"sample_policy" %in% colnames(previous_manifest)) {
      previous_manifest$sample_policy <- rep(NA_character_, nrow(previous_manifest))
    }
    previous_manifest <- previous_manifest[
      !previous_manifest$metric_file %in% manifest$metric_file, , drop = FALSE
    ]
    manifest <- rbind(previous_manifest, manifest)
  }
  write.table(manifest, manifest_path, sep = "\t", quote = FALSE,
              row.names = FALSE, na = "NA")
}
