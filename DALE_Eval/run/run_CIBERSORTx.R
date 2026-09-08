#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(argparse)
})

args0 <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args0, value = TRUE)
run_dir <- if (length(file_arg) > 0) dirname(normalizePath(sub("^--file=", "", file_arg[[1]]))) else getwd()
source(file.path(dirname(run_dir), "modules", "config_helpers.R"))
source(file.path(dirname(run_dir), "modules", "runner_helpers.R"))
source(file.path(dirname(run_dir), "modules", "marker_helpers.R"))
source(file.path(dirname(run_dir), "modules", "mapping_helpers.R"))

write_cibersortx_matrix <- function(x, path, row_label) {
  out <- data.frame(x, check.names = FALSE)
  out <- cbind(data.frame(row_label = rownames(x), check.names = FALSE), out)
  colnames(out)[[1]] <- row_label
  write.table(out, file = path, sep = "\t", row.names = FALSE, quote = FALSE)
}

cibersortx_gene_repair_candidates <- function(gene) {
  unique(c(
    sub("-([0-9]+)$", ".\\1", gene),
    sub("\\.", "-", gene),
    gsub("-", ".", gene),
    gsub("\\.", "-", gene),
    gene
  ))
}

read_cibersortx_cell_type_output <- function(path, genes, samples, label) {
  x <- read.delim(path, sep = "\t", check.names = FALSE, row.names = 1)
  x <- as.matrix(x)
  storage.mode(x) <- "numeric"
  x[is.na(x)] <- 0

  missing_samples <- setdiff(samples, colnames(x))
  if (length(missing_samples) > 0) {
    stop(label, " output missing samples: ", paste(head(missing_samples, 10), collapse = ", "))
  }

  output_genes <- rownames(x)
  repaired_from <- character(0)
  repaired_to <- character(0)
  dropped_output_genes <- character(0)
  gene_lookup <- setNames(genes, genes)

  for (i in seq_along(output_genes)) {
    gene <- output_genes[[i]]
    if (gene %in% gene_lookup) {
      next
    }

    candidates <- cibersortx_gene_repair_candidates(gene)
    matched <- candidates[candidates %in% gene_lookup]
    if (length(matched) > 0) {
      output_genes[[i]] <- matched[[1]]
      repaired_from <- c(repaired_from, gene)
      repaired_to <- c(repaired_to, matched[[1]])
    } else {
      dropped_output_genes <- c(dropped_output_genes, gene)
    }
  }
  rownames(x) <- output_genes

  keep <- rownames(x) %in% genes
  x <- x[keep, , drop = FALSE]
  if (length(dropped_output_genes) > 0) {
    message(label, " gene-name repair dropped unmatched output genes: ", length(dropped_output_genes))
    message("  drop examples: ", paste(head(dropped_output_genes, 3), collapse = "; "))
  }

  duplicated_genes <- unique(rownames(x)[duplicated(rownames(x))])
  if (length(duplicated_genes) > 0) {
    message(label, " gene-name repair dropped duplicated repaired output rows: ", length(duplicated_genes))
    message("  duplicate target examples: ", paste(head(duplicated_genes, 3), collapse = "; "))
    x <- x[!duplicated(rownames(x)), , drop = FALSE]
  }

  if (length(repaired_from) > 0) {
    examples <- paste0(head(repaired_from, 3), " -> ", head(repaired_to, 3))
    message(label, " gene-name repair renamed output genes: ", length(repaired_from))
    message("  repair examples: ", paste(examples, collapse = "; "))
  }

  matched_genes <- intersect(genes, rownames(x))
  list(
    matrix = x[matched_genes, samples, drop = FALSE],
    repaired_count = length(repaired_from),
    dropped_output_count = length(dropped_output_genes),
    matched_genes = matched_genes
  )
}

find_cibersortx_cell_type_file <- function(files, cell_type) {
  matches <- files[grepl(paste0(cell_type, "_Window"), basename(files), fixed = TRUE)]
  if (length(matches) != 1) {
    stop(
      "Expected exactly one CIBERSORTx HiRes output for cell type ", cell_type,
      "; found ", length(matches), ". Matches: ", paste(basename(matches), collapse = ", ")
    )
  }
  matches[[1]]
}

run_cibersortx_hires <- function(temp_dir, singularity_container_path, username, token) {
  args <- c(
    "exec", "-c",
    "-B", shQuote(paste0(temp_dir, ":/src/data")),
    "-B", shQuote(paste0(temp_dir, ":/src/outdir")),
    shQuote(singularity_container_path),
    "/src/CIBERSORTxHiRes",
    "--mixture", "bulk.txt",
    "--sigmatrix", "cbsx_sig.txt",
    "--cibresults", "frac.txt",
    "--label", "HiRes",
    "--threads", "1",
    "--heatmap", "FALSE",
    "--username", shQuote(username),
    "--token", shQuote(token)
  )
  out <- system2("singularity", args = args, stdout = TRUE, stderr = TRUE)
  exit_code <- attr(out, "status")
  if (length(out) > 0) {
    message(paste(out, collapse = "\n"))
  }
  if (is.null(exit_code)) 0L else exit_code
}

parser <- ArgumentParser(add_help = TRUE)
parser$add_argument("--dataset", type = "character", required = TRUE)
parser$add_argument("--config_id", type = "character", required = TRUE)
parser$add_argument("--n_core", type = "integer", required = FALSE, default = NULL)
parser$add_argument("--extra_args", type = "character", required = FALSE, default = "")

args <- parser$parse_args()
repo_root <- find_repo_root(run_dir)
extra <- apply_method_default_extra_args(parse_extra_args(args$extra_args), "CIBERSORTx", repo_root = repo_root)
validate_extra_args(
  extra,
  c("username", "token", "singularity_container_path", "map_cell_types", "use_test_samples", "use_limma_top_genes", "top_n", "keep_temp"),
  "CIBERSORTx"
)

dataset <- args$dataset
config_id <- args$config_id
n_core <- args$n_core
map_cell_types <- extra_arg(extra, "map_cell_types", TRUE)
use_test_samples <- extra_arg(extra, "use_test_samples", TRUE)
use_limma_top_genes <- extra_arg(extra, "use_limma_top_genes", FALSE)
top_n <- extra_arg(extra, "top_n", 100)
keep_temp <- extra_arg(extra, "keep_temp", FALSE)
validate_top_n_arg(extra, use_limma_top_genes, "CIBERSORTx")

default_container_path <- file.path(repo_root, "DALE_Eval", "external_modules", "CIBERSORTx", "hires.sif")
singularity_container_path <- extra_arg(extra, "singularity_container_path", default_container_path)
if (!grepl("^/", singularity_container_path)) {
  singularity_container_path <- file.path(repo_root, singularity_container_path)
}
username <- extra_arg(extra, "username")
token <- extra_arg(extra, "token")
if (is.null(username) || !nzchar(username)) {
  stop("CIBERSORTx requires username in extra_args")
}
if (is.null(token) || !nzchar(token)) {
  stop("CIBERSORTx requires token in extra_args")
}
if (!file.exists(singularity_container_path)) {
  stop("CIBERSORTx singularity container not found: ", singularity_container_path)
}

paths <- resolve_deconv_paths(dataset, config_id, repo_root = repo_root)
bulk_input <- paths$config$bulk_input[[1]]
bulk_scale <- paths$config$bulk_scale[[1]]
bulk_normalization <- paths$config$bulk_normalization[[1]]
frac_input <- paths$config$frac_input[[1]]
refType <- paths$refType
message_log_path <- method_message_log_path(paths, "CIBERSORTx")
cleanup_message_log <- start_message_log(message_log_path)
on.exit(cleanup_message_log(), add = TRUE)

bulk_expr <- read_bulk(paths)
bulk_prep <- prepare_bulk_for_deconv(bulk_expr, paths, method = "CIBERSORTx")
bulk_expr <- bulk_prep$bulk_expr
validate_non_log_bulk_input(bulk_expr, bulk_input, config_id, method = "CIBERSORTx")
frac <- read_frac_for_method(paths, repo_root = repo_root)
sample_restriction <- restrict_to_test_samples(
  bulk_expr,
  paths,
  frac = frac,
  use_test_samples = use_test_samples
)
bulk_expr <- sample_restriction$bulk
frac <- sample_restriction$frac
frac <- validate_bulk_frac(bulk_expr, frac)
bulk_input_dim <- format_deconv_dim(bulk_expr)

cbsx_sig <- read_reference_matrix(paths$ref_dir, "cbsx_sig.txt", "CIBERSORTx signature")
aligned <- align_fraction_and_reference_cell_types(
  frac,
  colnames(cbsx_sig),
  paths,
  repo_root
)
frac <- aligned$frac
cbsx_sig <- cbsx_sig[, aligned$reference_cell_types, drop = FALSE]

if (!identical(colnames(cbsx_sig), colnames(frac))) {
  stop("CIBERSORTx signature and fraction cell types are not aligned after subsetting")
}

signature_genes <- intersect(rownames(cbsx_sig), rownames(bulk_expr))
if (length(signature_genes) == 0) {
  stop("No CIBERSORTx signature genes overlap bulk input")
}
cbsx_sig <- cbsx_sig[signature_genes, , drop = FALSE]

limma_gene_count <- NA_integer_
limma_bulk_gene_count <- NA_integer_
if (use_limma_top_genes) {
  limma_genes <- read_limma_top_gene_union(paths, top_n = top_n)
  limma_gene_count <- length(limma_genes)
  limma_bulk_genes <- intersect(rownames(bulk_expr), limma_genes)
  limma_bulk_gene_count <- length(limma_bulk_genes)
  if (limma_bulk_gene_count == 0) {
    stop("No genes available for CIBERSORTx after applying limma top gene filter")
  }

  bulk_genes <- union(rownames(cbsx_sig), limma_bulk_genes)
  bulk_genes <- intersect(rownames(bulk_expr), bulk_genes)
  bulk_expr <- bulk_expr[bulk_genes, , drop = FALSE]
}

deconv_input_dim <- format_deconv_dim(bulk_expr)

message("Running CIBERSORTx")
message("  dataset: ", dataset)
message("  config_id: ", config_id)
message("  bulk_input: ", bulk_input, " (", bulk_input_dim, ")")
message("  bulk_scale: ", bulk_scale)
message("  bulk_normalization: ", bulk_normalization)
message("  bulk_preparation: ", bulk_prep$action)
message("  refType: ", refType)
message("  frac_input: ", frac_input, " (", nrow(frac), " samples x ", ncol(frac), " cell types)")
message("  map_cell_types: ", map_cell_types)
message("  use_test_samples: ", use_test_samples)
message("  use_limma_top_genes: ", use_limma_top_genes)
message("  keep_temp: ", keep_temp)
if (use_limma_top_genes) {
  message("  top_n: ", top_n)
  message("  limma_top_gene_union: ", limma_gene_count, " genes before intersecting bulk")
  message("  limma_top_gene_bulk_overlap: ", limma_bulk_gene_count, " genes")
  message("  note: limma mode keeps all CIBERSORTx signature genes in addition to limma top genes")
}
message("  signature_input: cbsx_sig.txt (", nrow(cbsx_sig), " genes x ", ncol(cbsx_sig), " cell types)")
message("  singularity_container_path: ", singularity_container_path)
if (!is.null(n_core) && !is.na(n_core)) {
  message("  n_core argument ignored: ", n_core)
}
message("  note: CIBERSORTx HiRes is forced to 1 thread because n_core > 1 is known to fail")

message("Starting deconvolution")
message("  deconv_input: ", deconv_input_dim)

with_runtime_log(paths, "CIBERSORTx", n_core = NA_integer_, expr = {
  if (keep_temp) {
    temp_dir <- file.path(
      paths$obj_dir,
      "logs",
      paste0("CIBERSORTx_temp_", config_id, "_", format(Sys.time(), "%Y%m%d_%H%M%S"))
    )
  } else {
    temp_dir <- tempfile("CIBERSORTx_")
  }
  dir.create(temp_dir, recursive = TRUE, showWarnings = FALSE)
  if (keep_temp) {
    message("  temp_dir: ", temp_dir)
  } else {
    on.exit(unlink(temp_dir, recursive = TRUE, force = TRUE), add = TRUE)
  }

  write_cibersortx_matrix(bulk_expr, file.path(temp_dir, "bulk.txt"), "GeneSymbol")
  write_cibersortx_matrix(frac, file.path(temp_dir, "frac.txt"), "Mixture")
  write_cibersortx_matrix(cbsx_sig, file.path(temp_dir, "cbsx_sig.txt"), "GeneSymbol")

  exit_code <- run_cibersortx_hires(temp_dir, singularity_container_path, username, token)
  if (!identical(exit_code, 0L)) {
    stop("CIBERSORTx singularity job failed with exit code ", exit_code)
  }

  output_files <- list.files(temp_dir, full.names = TRUE)
  cell_types <- colnames(cbsx_sig)
  cell_type_outputs <- list()

  for (cell_type in cell_types) {
    cell_type_path <- find_cibersortx_cell_type_file(output_files, cell_type)
    cell_type_outputs[[cell_type]] <- read_cibersortx_cell_type_output(
      cell_type_path,
      rownames(bulk_expr),
      colnames(bulk_expr),
      paste0("CIBERSORTx ", cell_type)
    )
  }

  common_genes <- Reduce(intersect, lapply(cell_type_outputs, function(x) x$matched_genes))
  common_genes <- rownames(bulk_expr)[rownames(bulk_expr) %in% common_genes]
  if (length(common_genes) == 0) {
    stop("No CIBERSORTx output genes matched the bulk input after gene-name repair")
  }
  dropped_bulk_genes <- setdiff(rownames(bulk_expr), common_genes)
  if (length(dropped_bulk_genes) > 0) {
    message("CIBERSORTx gene-name repair/export keeps ", length(common_genes), " of ", nrow(bulk_expr), " bulk genes")
    message("  dropped bulk gene examples: ", paste(head(dropped_bulk_genes, 3), collapse = "; "))
  }

  Z <- array(
    NA_real_,
    dim = c(length(common_genes), ncol(bulk_expr), length(cell_types)),
    dimnames = list(common_genes, colnames(bulk_expr), cell_types)
  )
  for (cell_type in cell_types) {
    Z[, , cell_type] <- cell_type_outputs[[cell_type]]$matrix[common_genes, colnames(bulk_expr), drop = FALSE]
  }

  Z_export <- Z
  if (map_cell_types) {
    Z_export <- apply_indep_ref_cell_type_mapping(Z_export, paths, repo_root = repo_root)
  }
  message_existing_method_output(paths, "CIBERSORTx")
  write_ctse_array(Z_export, paths$output_dir, "CIBERSORTx")
}, deconv_input = deconv_input_dim)

message("Output saved")
message("  CIBERSORTx ctse: ", file.path(paths$output_dir, "CIBERSORTx", "<cell_type>.txt.gz"))
if (map_cell_types && refType == "indep") {
  message("  ctse naming note: mapped cell types use benchmark names; unmapped reference cell types keep reference names")
} else {
  message("  ctse naming note: exported cell types use reference cell type names")
}
message("done")
