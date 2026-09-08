# Build unified dataset/cell-type annotation for performance summaries.
#
# Notebook-style workflow:
#   1. Edit the settings block if the dataset/tissue map changes.
#   2. Source or run this script from repo root or benchmark_summary/preparation/.
#   3. Inspect DALE_Eval/configs/benchmark_dataset_info.txt.
#
# Output:
#   DALE_Eval/configs/benchmark_dataset_info.txt
#
# Columns:
#   dataset, cell_type_raw, cell_type, ct_abundance, transcript_abundance,
#   sample_size, n_eval_samples_meancpm, n_eval_samples_sumcount_cpm,
#   n_truth_genes, n_genes_indep_ref, n_genes_truth_indep_ref_overlap,
#   n_eval_genes_meancpm, n_eval_genes_sumcount_cpm, tissue

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

output_path <- file.path(
  repo_root,
  "DALE_Eval",
  "configs",
  "benchmark_dataset_info.txt"
)

# Canonical fraction threshold used by active CTSE correlation evaluation.
eval_min_frac <- 0.001

dataset_tissue <- data.frame(
  dataset = c(
    "BRCA_Bassez2021",
    "CRC_Pelka2021",
    "LUAD_Kim2020",
    "PBMC_1k1k",
    "PBMC_Perez2022",
    "PBMC_refined_1k1k",
    "PBMC_refined_Perez2022",
    "ROSMAP_AD430_Mathys2023",
    "ROSMAP_AD92_Xiong2023"
  ),
  tissue = c(
    "Tumor",
    "Tumor",
    "Tumor",
    "Blood",
    "Blood",
    "Blood_refined",
    "Blood_refined",
    "Brain",
    "Brain"
  ),
  stringsAsFactors = FALSE
)

tumor_lineage_alias <- c(
  Malignant = "Malignant",
  Cancer_cell = "Malignant",
  Myeloid_cell = "Myeloid",
  Myeloid = "Myeloid",
  T_cell = "T_cell",
  Fibroblast = "Fibroblasts",
  NK_cell = "NK",
  NK = "NK",
  Endothelial = "Endothelial",
  Endothelial_cell = "Endothelial",
  B_cell = "B_cell",
  B = "B_cell",
  Plasma_cell = "Plasma_B",
  Plasma = "Plasma_B",
  pDC = "DC",
  Dendritic = "DC",
  Dendritic_cell = "DC",
  Mast_cell = "Mast",
  Mast = "Mast",
  Epithelial_cell = "Epithelial"
)

## ---------------------------------------------------------------------------
## Small Helpers
## ---------------------------------------------------------------------------

read_numeric_matrix <- function(path, required = TRUE) {
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

validate_gene_names <- function(genes, label) {
  if (length(genes) == 0L) {
    stop(label, ": no genes found")
  }

  invalid <- is.na(genes) | !nzchar(trimws(genes))
  if (any(invalid)) {
    stop(label, ": blank or missing gene identifiers")
  }

  duplicated_genes <- unique(genes[duplicated(genes)])
  if (length(duplicated_genes) > 0L) {
    stop(
      label,
      ": duplicated gene identifiers: ",
      paste(head(duplicated_genes, 5L), collapse = ", ")
    )
  }

  genes
}

read_matrix_gene_names <- function(path, sep) {
  if (!file.exists(path)) {
    stop("Missing gene matrix: ", path)
  }

  con <- if (grepl("\\.gz$", path)) {
    gzfile(path, open = "rt")
  } else {
    file(path, open = "rt")
  }
  header <- tryCatch(
    readLines(con, n = 1L, warn = FALSE),
    finally = close(con)
  )
  if (length(header) != 1L) {
    stop("Missing matrix header: ", path)
  }

  n_fields <- length(strsplit(header, sep, fixed = TRUE)[[1]])
  if (n_fields < 2L) {
    stop("Expected a row-name column and at least one data column: ", path)
  }

  matrix_stub <- read.delim(
    path,
    sep = sep,
    check.names = FALSE,
    row.names = 1,
    colClasses = c("character", rep("NULL", n_fields - 1L))
  )

  validate_gene_names(rownames(matrix_stub), path)
}

read_truth_gene_universe <- function(dataset_root, dataset_name) {
  truth_root <- file.path(dataset_root, "ctse_truth")
  truth_paths <- sort(list.files(
    truth_root,
    pattern = "\\.txt\\.gz$",
    recursive = TRUE,
    full.names = TRUE
  ))
  if (length(truth_paths) == 0L) {
    stop(dataset_name, ": no CTSE truth matrices found under ", truth_root)
  }

  truth_genes <- read_matrix_gene_names(truth_paths[[1]], sep = "\t")
  if (length(truth_paths) > 1L) {
    for (truth_path in truth_paths[-1]) {
      genes <- read_matrix_gene_names(truth_path, sep = "\t")
      if (!setequal(genes, truth_genes)) {
        stop(
          dataset_name,
          ": truth gene set differs in ",
          truth_path
        )
      }
    }
  }

  truth_genes
}

read_indep_ref_assignments <- function(path) {
  if (!file.exists(path)) {
    stop("Missing independent-reference assignments: ", path)
  }

  assignments <- read.delim(
    path,
    sep = "\t",
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  required_cols <- c("dataset", "indep_ref")
  missing_cols <- setdiff(required_cols, colnames(assignments))
  if (length(missing_cols) > 0L) {
    stop(
      "benchmark_ref_assignment.txt missing columns: ",
      paste(missing_cols, collapse = ", ")
    )
  }
  if (anyDuplicated(assignments$dataset)) {
    stop("benchmark_ref_assignment.txt has duplicated dataset assignments")
  }

  assignments
}

count_evaluable_genes <- function(dataset_root,
                                  dataset_name,
                                  truth_type,
                                  truth_genes,
                                  cell_types) {
  truth_dir <- file.path(dataset_root, "ctse_truth", truth_type)
  if (!dir.exists(truth_dir)) {
    return(setNames(rep(NA_integer_, length(cell_types)), cell_types))
  }

  truth_paths <- list.files(
    truth_dir,
    pattern = "\\.txt\\.gz$",
    full.names = TRUE
  )
  if (length(truth_paths) == 0L) {
    return(setNames(rep(NA_integer_, length(cell_types)), cell_types))
  }

  mask_path <- file.path(
    dataset_root,
    "evaluation_metadata",
    "cor_na_acceptable_mask",
    paste0("cor_na_acceptable_mask_truth-", truth_type, ".txt")
  )
  mask <- read_numeric_matrix(mask_path)

  if (!setequal(rownames(mask), truth_genes)) {
    stop(dataset_name, " / ", truth_type, ": mask genes differ from truth genes")
  }
  if (!setequal(colnames(mask), cell_types)) {
    stop(dataset_name, " / ", truth_type, ": mask cell types differ from truth fractions")
  }
  valid_mask <- is.finite(mask) & mask %in% c(0, 1)
  if (!all(valid_mask)) {
    stop(dataset_name, " / ", truth_type, ": mask contains values other than 0 or 1")
  }

  counts <- colSums(mask[, cell_types, drop = FALSE] == 0)
  setNames(as.integer(counts), cell_types)
}

read_matrix_sample_names <- function(path, required = FALSE) {
  if (!file.exists(path)) {
    if (required) {
      stop("Missing matrix: ", path)
    }
    return(NULL)
  }

  matrix_header <- read.delim(
    path,
    sep = "	",
    check.names = FALSE,
    row.names = 1,
    nrows = 1
  )
  colnames(matrix_header)
}

count_truth_eval_samples <- function(truth_path,
                                     frac,
                                     cell_type,
                                     test_samples,
                                     min_frac = eval_min_frac) {
  if (is.null(frac) || !cell_type %in% colnames(frac)) {
    return(NA_integer_)
  }

  truth_samples <- read_matrix_sample_names(truth_path, required = FALSE)
  if (is.null(truth_samples)) {
    return(NA_integer_)
  }

  filter_values <- frac[, cell_type]
  filter_samples <- rownames(frac)[
    is.finite(filter_values) & filter_values > min_frac
  ]

  eval_samples <- Reduce(
    intersect,
    list(test_samples, truth_samples, filter_samples)
  )

  as.integer(length(eval_samples))
}

read_test_samples <- function(path) {
  if (!file.exists(path)) {
    stop("Missing sample split: ", path)
  }

  split <- read.delim(path, sep = "\t", check.names = FALSE)
  required_cols <- c("group", "sampleIDs")
  missing_cols <- setdiff(required_cols, colnames(split))
  if (length(missing_cols) > 0L) {
    stop("sample_split.txt missing columns: ", paste(missing_cols, collapse = ", "))
  }

  samples <- split$sampleIDs[split$group == "test"]
  samples[!is.na(samples) & samples != ""]
}

standardize_cell_type <- function(cell_type_raw, tissue) {
  if (identical(tissue, "Tumor") && cell_type_raw %in% names(tumor_lineage_alias)) {
    return(unname(tumor_lineage_alias[[cell_type_raw]]))
  }
  cell_type_raw
}

format_numeric <- function(x) {
  ifelse(is.na(x), NA, signif(x, digits = 10))
}

## ---------------------------------------------------------------------------
## Build Table
## ---------------------------------------------------------------------------

benchmark_root <- file.path(repo_root, "Benchmarking_obj")
ref_assignments <- read_indep_ref_assignments(file.path(
  repo_root,
  "DALE_Eval",
  "configs",
  "benchmark_ref_assignment.txt"
))
rows <- list()

for (i in seq_len(nrow(dataset_tissue))) {
  dataset_name <- dataset_tissue$dataset[[i]]
  tissue <- dataset_tissue$tissue[[i]]
  dataset_root <- file.path(benchmark_root, dataset_name)

  message("Dataset: ", dataset_name)

  cellfrac <- read_numeric_matrix(
    file.path(dataset_root, "frac_input", "truth_cellfrac.txt")
  )
  transcriptfrac <- read_numeric_matrix(
    file.path(dataset_root, "frac_input", "truth_transcriptfrac.txt"),
    required = FALSE
  )
  test_samples <- read_test_samples(
    file.path(dataset_root, "self_reference", "sample_split.txt")
  )

  truth_genes <- read_truth_gene_universe(dataset_root, dataset_name)
  n_truth_genes <- as.integer(length(truth_genes))

  ref_assignment <- ref_assignments[
    ref_assignments$dataset == dataset_name,
    ,
    drop = FALSE
  ]
  if (nrow(ref_assignment) != 1L) {
    stop(dataset_name, ": expected exactly one independent-reference assignment")
  }
  indep_ref <- ref_assignment$indep_ref[[1]]
  indep_ref_genes <- read_matrix_gene_names(
    file.path(
      repo_root,
      "Indep_scReference",
      indep_ref,
      "rowMeans_sig.csv"
    ),
    sep = ","
  )
  n_genes_indep_ref <- as.integer(length(indep_ref_genes))
  n_genes_truth_indep_ref_overlap <- as.integer(length(intersect(
    truth_genes,
    indep_ref_genes
  )))

  eval_genes_meancpm <- count_evaluable_genes(
    dataset_root = dataset_root,
    dataset_name = dataset_name,
    truth_type = "meancpm",
    truth_genes = truth_genes,
    cell_types = colnames(cellfrac)
  )
  eval_genes_sumcount_cpm <- count_evaluable_genes(
    dataset_root = dataset_root,
    dataset_name = dataset_name,
    truth_type = "sumcount_cpm",
    truth_genes = truth_genes,
    cell_types = colnames(cellfrac)
  )

  cellfrac_samples <- intersect(test_samples, rownames(cellfrac))
  if (length(cellfrac_samples) == 0L) {
    stop(dataset_name, ": no test samples overlap truth_cellfrac.txt")
  }

  for (cell_type_raw in colnames(cellfrac)) {
    cell_type <- standardize_cell_type(cell_type_raw, tissue)
    ct_abundance <- mean(cellfrac[cellfrac_samples, cell_type_raw], na.rm = TRUE)

    transcript_abundance <- NA_real_
    if (!is.null(transcriptfrac) && cell_type_raw %in% colnames(transcriptfrac)) {
      transcript_samples <- intersect(cellfrac_samples, rownames(transcriptfrac))
      if (length(transcript_samples) > 0L) {
        transcript_abundance <- mean(
          transcriptfrac[transcript_samples, cell_type_raw],
          na.rm = TRUE
        )
      }
    }

    n_eval_samples_meancpm <- count_truth_eval_samples(
      truth_path = file.path(
        dataset_root,
        "ctse_truth",
        "meancpm",
        paste0(cell_type_raw, ".txt.gz")
      ),
      frac = cellfrac,
      cell_type = cell_type_raw,
      test_samples = test_samples
    )
    n_eval_samples_sumcount_cpm <- count_truth_eval_samples(
      truth_path = file.path(
        dataset_root,
        "ctse_truth",
        "sumcount_cpm",
        paste0(cell_type_raw, ".txt.gz")
      ),
      frac = transcriptfrac,
      cell_type = cell_type_raw,
      test_samples = test_samples
    )

    rows[[length(rows) + 1L]] <- data.frame(
      dataset = dataset_name,
      cell_type_raw = cell_type_raw,
      cell_type = cell_type,
      ct_abundance = format_numeric(ct_abundance),
      transcript_abundance = format_numeric(transcript_abundance),
      sample_size = length(cellfrac_samples),
      n_eval_samples_meancpm = n_eval_samples_meancpm,
      n_eval_samples_sumcount_cpm = n_eval_samples_sumcount_cpm,
      n_truth_genes = n_truth_genes,
      n_genes_indep_ref = n_genes_indep_ref,
      n_genes_truth_indep_ref_overlap = n_genes_truth_indep_ref_overlap,
      n_eval_genes_meancpm = eval_genes_meancpm[[cell_type_raw]],
      n_eval_genes_sumcount_cpm = eval_genes_sumcount_cpm[[cell_type_raw]],
      tissue = tissue,
      stringsAsFactors = FALSE
    )
  }
}

dataset_info <- do.call(rbind, rows)

row_keys <- paste(dataset_info$dataset, dataset_info$cell_type_raw, sep = "::")
if (anyDuplicated(row_keys)) {
  stop("Dataset information has duplicated dataset/cell_type_raw rows")
}

dataset_level_gene_cols <- c(
  "n_truth_genes",
  "n_genes_indep_ref",
  "n_genes_truth_indep_ref_overlap"
)
for (dataset_name in unique(dataset_info$dataset)) {
  dataset_rows <- dataset_info$dataset == dataset_name
  for (column in dataset_level_gene_cols) {
    if (length(unique(dataset_info[[column]][dataset_rows])) != 1L) {
      stop(dataset_name, ": ", column, " is not constant across cell types")
    }
  }
}

invalid_overlap <- dataset_info$n_genes_truth_indep_ref_overlap > pmin(
  dataset_info$n_truth_genes,
  dataset_info$n_genes_indep_ref
)
if (any(invalid_overlap)) {
  stop("Truth/independent-reference overlap exceeds a source gene count")
}

for (column in c("n_eval_genes_meancpm", "n_eval_genes_sumcount_cpm")) {
  invalid <- !is.na(dataset_info[[column]]) & (
    dataset_info[[column]] < 0L |
      dataset_info[[column]] > dataset_info$n_truth_genes
  )
  if (any(invalid)) {
    stop(column, " is outside [0, n_truth_genes]")
  }
}

write.table(
  dataset_info,
  file = output_path,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  na = "NA"
)

message("Wrote: ", output_path)
