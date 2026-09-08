# Build config-specific high-fraction sample masks for CTSE correlation tests.
#
# Notebook-style workflow:
#   1. Edit dataset_include or thresholds below if needed.
#   2. Source or run this script from repo root or benchmark_summary/preparation/.
#   3. Inspect the masks under each dataset's evaluation_metadata/sample_mask/.
#
# Two masks are built when their source files are available:
#
#   config01:
#     test sample AND config01 InstaPrismfrac > 0.1
#
#   config02:
#     test sample AND config02 InstaPrismfrac > 0.1
#
# Output files:
#   Benchmarking_obj/<dataset>/evaluation_metadata/sample_mask/
#     InstaPrismfrac_config01_gt0.1.txt
#     InstaPrismfrac_config02_gt0.1.txt
#
# Each output is a sample-by-target-cell-type 0/1 matrix containing test samples
# present in the config-specific InstaPrism fraction matrix. A value of 1 means
# that the mapped InstaPrism fraction is strictly greater than 0.1. The later
# evaluator will apply its normal truth-fraction filter separately.
#
# InstaPrism independent-reference cell types are mapped to dataset target cell
# types using Indep_scReference/cell_type_mapping.txt. Unmapped cell types are
# excluded.
#
# A mask column is retained even when fewer than min_n_sample samples pass.
# The later evaluator should use min_n_sample=10 and return NA correlations
# for that cell type when the retained sample count is below this threshold.

## ---------------------------------------------------------------------------
## Settings
## ---------------------------------------------------------------------------

repo_root_candidates <- unique(normalizePath(
  c(".", "..", "../.."),
  mustWork = FALSE
))
repo_root_matches <- repo_root_candidates[
  file.exists(file.path(repo_root_candidates, "Benchmarking_obj"))
]
if (length(repo_root_matches) == 0L) {
  stop("Run from repo root or benchmark_summary/preparation/")
}
repo_root <- normalizePath(repo_root_matches[[1]])

dataset_include <- NULL
instaprism_min_frac <- 0.1
min_n_sample <- 10L

config_ids <- c("config01", "config02")

## ---------------------------------------------------------------------------
## Setup
## ---------------------------------------------------------------------------

source(file.path(repo_root, "DALE_Eval", "modules", "config_helpers.R"))

read_numeric_matrix <- function(path) {
  if (!file.exists(path)) {
    stop("Missing matrix: ", path)
  }

  x <- read.delim(path, sep = "\t", check.names = FALSE, row.names = 1)
  x <- as.matrix(x)
  storage.mode(x) <- "numeric"
  x
}

read_test_samples <- function(path) {
  if (!file.exists(path)) {
    stop("Missing sample split: ", path)
  }

  split <- read.delim(path, sep = "\t", check.names = FALSE)
  required_cols <- c("group", "sampleIDs")
  missing_cols <- setdiff(required_cols, colnames(split))
  if (length(missing_cols) > 0L) {
    stop(
      "sample_split.txt missing columns: ",
      paste(missing_cols, collapse = ", ")
    )
  }

  samples <- split$sampleIDs[split$group == "test"]
  unique(samples[!is.na(samples) & samples != ""])
}

read_cell_type_mapping <- function(path) {
  if (!file.exists(path)) {
    stop("Missing cell-type mapping: ", path)
  }

  mapping <- read.delim(
    path,
    sep = "\t",
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  required_cols <- c(
    "dataset",
    "indep_ref",
    "target_cell_type",
    "indep_ref_cell_type"
  )
  missing_cols <- setdiff(required_cols, colnames(mapping))
  if (length(missing_cols) > 0L) {
    stop(
      "cell_type_mapping.txt missing columns: ",
      paste(missing_cols, collapse = ", ")
    )
  }

  mapping
}

write_sample_mask <- function(mask, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  if (file.exists(path)) {
    warning("Overwriting existing sample mask: ", path)
  }

  write.table(
    mask,
    file = path,
    sep = "\t",
    quote = FALSE,
    col.names = NA
  )
}

format_cutoff <- function(value) {
  format(value, scientific = FALSE, trim = TRUE)
}

config_rows <- setNames(
  lapply(config_ids, function(config_id) {
    config_row <- read_deconv_config(config_id, repo_root = repo_root)
    if (config_row$refType[[1]] != "indep") {
      stop(config_id, " must use refType=indep")
    }
    config_row
  }),
  config_ids
)

mapping_all <- read_cell_type_mapping(file.path(
  repo_root,
  "Indep_scReference",
  "cell_type_mapping.txt"
))

benchmark_root <- file.path(repo_root, "Benchmarking_obj")
datasets <- sort(list.dirs(
  benchmark_root,
  full.names = FALSE,
  recursive = FALSE
))
if (!is.null(dataset_include)) {
  datasets <- intersect(datasets, dataset_include)
}
if (length(datasets) == 0L) {
  stop("No datasets selected")
}

## ---------------------------------------------------------------------------
## Build Masks
## ---------------------------------------------------------------------------

for (dataset_name in datasets) {
  dataset_root <- file.path(benchmark_root, dataset_name)
  sample_split_path <- file.path(
    dataset_root,
    "self_reference",
    "sample_split.txt"
  )
  if (!file.exists(sample_split_path)) {
    message("Dataset: ", dataset_name, " [skip: missing sample_split.txt]")
    next
  }

  message("Dataset: ", dataset_name)
  test_samples <- read_test_samples(sample_split_path)

  indep_ref_dir <- tryCatch(
    resolve_indep_ref_dir(dataset_name, repo_root = repo_root),
    error = function(error) error
  )
  if (inherits(indep_ref_dir, "error")) {
    message("  skip: ", conditionMessage(indep_ref_dir))
    next
  }
  indep_ref <- basename(indep_ref_dir)

  mapping <- mapping_all[
    mapping_all$dataset == dataset_name &
      mapping_all$indep_ref == indep_ref,
    c("target_cell_type", "indep_ref_cell_type"),
    drop = FALSE
  ]
  mapping <- unique(mapping)
  if (nrow(mapping) == 0L) {
    message("  skip: no target/independent-reference cell-type mappings")
    next
  }
  if (anyDuplicated(mapping$target_cell_type)) {
    stop(
      dataset_name,
      ": multiple independent-reference cell types map to one target cell type"
    )
  }

  for (config_id in config_ids) {
    config_row <- config_rows[[config_id]]
    slug <- config_slug(config_row)

    instaprism_frac_path <- file.path(
      dataset_root,
      "deconv_res",
      slug,
      "InstaPrism",
      "InstaPrismfrac.txt"
    )

    if (!file.exists(instaprism_frac_path)) {
      message("  ", config_id, " skip: missing ", instaprism_frac_path)
      next
    }

    instaprism_frac <- read_numeric_matrix(instaprism_frac_path)

    mapping_one <- mapping[
      mapping$indep_ref_cell_type %in% colnames(instaprism_frac),
      ,
      drop = FALSE
    ]
    if (nrow(mapping_one) == 0L) {
      message(
        "  ", config_id,
        " skip: no mapped cell types overlap the InstaPrism fraction matrix"
      )
      next
    }

    common_samples <- test_samples[
      test_samples %in% rownames(instaprism_frac)
    ]
    if (length(common_samples) == 0L) {
      message("  ", config_id, " skip: no overlapping test samples")
      next
    }

    mask <- matrix(
      0L,
      nrow = length(common_samples),
      ncol = nrow(mapping_one),
      dimnames = list(common_samples, mapping_one$target_cell_type)
    )

    for (cell_type_index in seq_len(nrow(mapping_one))) {
      target_cell_type <- mapping_one$target_cell_type[[cell_type_index]]
      indep_ref_cell_type <- mapping_one$indep_ref_cell_type[[cell_type_index]]

      instaprism_values <- instaprism_frac[
        common_samples,
        indep_ref_cell_type
      ]
      retained <- is.finite(instaprism_values) &
        instaprism_values > instaprism_min_frac
      mask[, target_cell_type] <- as.integer(retained)

      message(
        "    ", config_id, " / ", target_cell_type,
        ": retained=", sum(retained),
        "/", length(common_samples),
        if (sum(retained) < min_n_sample) " [below min_n_sample]" else ""
      )
    }

    output_path <- file.path(
      dataset_root,
      "evaluation_metadata",
      "sample_mask",
      paste0(
        "InstaPrismfrac_",
        config_id,
        "_gt",
        format_cutoff(instaprism_min_frac),
        ".txt"
      )
    )
    write_sample_mask(mask, output_path)
    message(
      "  wrote ", output_path,
      " (", nrow(mask), " samples x ", ncol(mask), " cell types)"
    )
  }
}

message("Done.")
