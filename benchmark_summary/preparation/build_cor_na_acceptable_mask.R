# Build test-sample truth/filter-derived acceptable-NA masks for correlation metrics.
#
# Notebook-style workflow:
#   1. Edit the settings block.
#   2. Source or run this script from repo root or benchmark_summary/preparation/.
#   3. Inspect messages and output files.
#
# Output:
#   Benchmarking_obj/<dataset>/evaluation_metadata/cor_na_acceptable_mask/
#     cor_na_acceptable_mask_truth-<truth_type>.txt
#
# Matrix values:
#   1 = correlation NA is acceptable and should remain NA
#       (too few evaluable samples or truth CTSE is constant)
#   0 = correlation NA is not explained by truth/filter constraints

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
truth_types <- c("meancpm", "sumcount_cpm")

filter_frac_by_truth <- c(
  meancpm = "truth_cellfrac",
  sumcount_cpm = "truth_transcriptfrac"
)

min_frac <- 0.001
min_n_sample <- 10L
output_folder_name <- "cor_na_acceptable_mask"

## ---------------------------------------------------------------------------
## Small Helpers
## ---------------------------------------------------------------------------

read_numeric_matrix <- function(path) {
  x <- read.delim(path, sep = "\t", check.names = FALSE, row.names = 1)
  x <- as.matrix(x)
  storage.mode(x) <- "numeric"
  x[is.na(x)] <- 0
  x
}

read_test_samples <- function(path) {
  if (!file.exists(path)) {
    stop("Missing sample split: ", path)
  }

  split <- read.delim(path, sep = "	", check.names = FALSE)
  required_cols <- c("group", "sampleIDs")
  missing_cols <- setdiff(required_cols, colnames(split))
  if (length(missing_cols) > 0L) {
    stop(
      "sample_split.txt missing columns: ",
      paste(missing_cols, collapse = ", ")
    )
  }

  is_test <- !is.na(split$group) & trimws(as.character(split$group)) == "test"
  samples <- trimws(as.character(split$sampleIDs[is_test]))
  samples <- unique(samples[!is.na(samples) & nzchar(samples)])

  if (length(samples) == 0L) {
    stop("No test samples found in sample split: ", path)
  }

  samples
}

safe_cell_type_name <- function(path) {
  sub("\\.txt\\.gz$", "", basename(path))
}

truth_is_constant <- function(x) {
  x <- x[is.finite(x)]
  length(unique(x)) <= 1L
}

write_mask_txt <- function(mask, path) {
  write.table(
    mask,
    file = path,
    sep = "\t",
    quote = FALSE,
    col.names = NA
  )
}

## ---------------------------------------------------------------------------
## Discover Datasets
## ---------------------------------------------------------------------------

benchmark_root <- file.path(repo_root, "Benchmarking_obj")
datasets <- list.dirs(benchmark_root, full.names = FALSE, recursive = FALSE)
datasets <- datasets[file.exists(file.path(benchmark_root, datasets, "ctse_truth"))]

if (!is.null(dataset_include)) {
  datasets <- intersect(datasets, dataset_include)
}

if (length(datasets) == 0) {
  stop("No datasets selected")
}

## ---------------------------------------------------------------------------
## Build Masks
## ---------------------------------------------------------------------------

for (dataset_name in datasets) {
  message("Dataset: ", dataset_name)

  dataset_root <- file.path(benchmark_root, dataset_name)
  frac_root <- file.path(dataset_root, "frac_input")
  out_dir <- file.path(dataset_root, "evaluation_metadata", output_folder_name)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  test_samples <- read_test_samples(
    file.path(dataset_root, "self_reference", "sample_split.txt")
  )

  for (truth_type in truth_types) {
    filter_frac <- filter_frac_by_truth[[truth_type]]
    if (is.null(filter_frac) || is.na(filter_frac)) {
      message("  skip ", truth_type, ": no default filter fraction configured")
      next
    }

    truth_dir <- file.path(dataset_root, "ctse_truth", truth_type)
    frac_path <- file.path(frac_root, paste0(filter_frac, ".txt"))

    if (!dir.exists(truth_dir)) {
      message("  skip ", truth_type, ": missing truth dir ", truth_dir)
      next
    }
    if (!file.exists(frac_path)) {
      message("  skip ", truth_type, ": missing fraction file ", frac_path)
      next
    }

    truth_paths <- sort(list.files(truth_dir, pattern = "\\.txt\\.gz$", full.names = TRUE))
    if (length(truth_paths) == 0) {
      message("  skip ", truth_type, ": no truth .txt.gz files")
      next
    }

    frac <- read_numeric_matrix(frac_path)
    truth_files <- setNames(truth_paths, vapply(truth_paths, safe_cell_type_name, character(1)))
    cell_types <- intersect(names(truth_files), colnames(frac))

    if (length(cell_types) == 0) {
      message("  skip ", truth_type, ": no overlap between truth cell types and ", filter_frac)
      next
    }

    truth_by_cell_type <- list()
    all_genes <- character(0)

    for (cell_type in cell_types) {
      truth <- read_numeric_matrix(truth_files[[cell_type]])
      truth_by_cell_type[[cell_type]] <- truth
      all_genes <- union(all_genes, rownames(truth))
    }

    all_genes <- sort(all_genes)
    mask <- matrix(
      1L,
      nrow = length(all_genes),
      ncol = length(cell_types),
      dimnames = list(all_genes, cell_types)
    )

    for (cell_type in cell_types) {
      truth <- truth_by_cell_type[[cell_type]]
      filter_values <- frac[, cell_type]
      filter_samples <- rownames(frac)[is.finite(filter_values) & filter_values > min_frac]
      eval_samples <- Reduce(
        intersect,
        list(test_samples, colnames(truth), filter_samples)
      )

      if (length(eval_samples) < min_n_sample) {
        message(
          "  ", truth_type, " / ", cell_type,
          ": all genes acceptable NA; eval_samples=", length(eval_samples)
        )
        next
      }

      for (gene in rownames(truth)) {
        truth_values <- truth[gene, eval_samples]
        if (!truth_is_constant(truth_values)) {
          mask[gene, cell_type] <- 0L
        }
      }

      message(
        "  ", truth_type, " / ", cell_type,
        ": genes=", nrow(truth),
        ", eval_samples=", length(eval_samples),
        ", acceptable_NA=", sum(mask[rownames(truth), cell_type] == 1L)
      )
    }

    out_path <- file.path(out_dir, paste0("cor_na_acceptable_mask_truth-", truth_type, ".txt"))
    write_mask_txt(mask, out_path)

    message(
      "  wrote ", out_path,
      " (", nrow(mask), " genes x ", ncol(mask), " cell types)"
    )
  }
}

message("Done.")
