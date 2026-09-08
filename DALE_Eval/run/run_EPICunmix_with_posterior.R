#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(argparse)
  library(EPICunmix)
})

args0 <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args0, value = TRUE)
run_dir <- if (length(file_arg) > 0) dirname(normalizePath(sub("^--file=", "", file_arg[[1]]))) else getwd()
source(file.path(dirname(run_dir), "modules", "config_helpers.R"))
source(file.path(dirname(run_dir), "modules", "runner_helpers.R"))
source(file.path(dirname(run_dir), "modules", "marker_helpers.R"))
source(file.path(dirname(run_dir), "modules", "mapping_helpers.R"))

posterior_gene_names <- function(posterior) {
  if (!is.null(posterior$A) && length(dim(posterior$A)) >= 1) {
    genes <- dimnames(posterior$A)[[1]]
    if (!is.null(genes)) return(genes)
  }
  if (!is.null(posterior$mu) && length(dim(posterior$mu)) >= 1) {
    genes <- dimnames(posterior$mu)[[1]]
    if (!is.null(genes)) return(genes)
  }
  stop("bMIND posterior must include gene dimnames on posterior$A or posterior$mu")
}

posterior_cell_type_names <- function(posterior) {
  if (!is.null(posterior$A) && length(dim(posterior$A)) >= 2) {
    cell_types <- dimnames(posterior$A)[[2]]
    if (!is.null(cell_types)) return(cell_types)
  }
  if (!is.null(posterior$mu) && length(dim(posterior$mu)) >= 2) {
    cell_types <- dimnames(posterior$mu)[[2]]
    if (!is.null(cell_types)) return(cell_types)
  }
  stop("bMIND posterior must include cell-type dimnames on posterior$A or posterior$mu")
}

subset_posterior_genes <- function(posterior, all_genes, genes_keep) {
  subset_gene_dim <- function(x) {
    if (is.list(x)) {
      return(lapply(x, subset_gene_dim))
    }

    dx <- dim(x)
    if (is.null(dx)) {
      if (!is.null(names(x)) && all(genes_keep %in% names(x))) {
        return(x[genes_keep])
      }
      return(x)
    }

    dn <- dimnames(x)
    gene_dim <- integer(0)
    if (!is.null(dn)) {
      gene_dim <- which(vapply(dn, function(z) {
        !is.null(z) && identical(z, all_genes)
      }, logical(1)))
    }

    if (length(gene_dim) == 0 && dx[[1]] == length(all_genes)) {
      gene_dim <- 1
    }

    if (length(gene_dim) == 0) {
      return(x)
    }

    idx <- rep(list(TRUE), length(dx))
    idx[[gene_dim[[1]]]] <- genes_keep
    do.call(`[`, c(list(x), idx, list(drop = FALSE)))
  }

  subset_gene_dim(posterior)
}

split_gene_chunks <- function(genes, chunk_size) {
  if (chunk_size < 1) {
    stop("chunk_size must be a positive integer")
  }
  split(genes, ceiling(seq_along(genes) / chunk_size))
}

write_chunk_status <- function(status, chunk_root) {
  write.table(
    status,
    file = file.path(chunk_root, "chunk_status.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
}

merge_chunk_outputs <- function(chunk_output_dirs, final_temp_dir) {
  if (length(chunk_output_dirs) == 0) {
    stop("No successful EPICunmix chunk output directories to merge")
  }

  cell_files_by_chunk <- lapply(chunk_output_dirs, function(d) {
    sort(list.files(d, pattern = "\\.txt\\.gz$", full.names = FALSE))
  })
  cell_files <- cell_files_by_chunk[[1]]
  if (length(cell_files) == 0) {
    stop("No CTSE files found in first chunk output directory: ", chunk_output_dirs[[1]])
  }
  for (i in seq_along(cell_files_by_chunk)) {
    if (!identical(cell_files_by_chunk[[i]], cell_files)) {
      stop("Chunk output cell-type files differ for chunk output directory: ", chunk_output_dirs[[i]])
    }
  }

  dir.create(final_temp_dir, recursive = TRUE, showWarnings = FALSE)
  for (cell_file in cell_files) {
    mats <- lapply(chunk_output_dirs, function(method_dir) {
      path <- file.path(method_dir, cell_file)
      read.delim(path, sep = "\t", check.names = FALSE, row.names = 1)
    })

    col_ref <- colnames(mats[[1]])
    for (mat in mats) {
      if (!identical(colnames(mat), col_ref)) {
        stop("Sample columns differ across chunks for ", cell_file)
      }
    }

    merged <- do.call(rbind, mats)
    if (anyDuplicated(rownames(merged))) {
      stop("Duplicated genes after merging chunks for ", cell_file)
    }

    out_path <- file.path(final_temp_dir, cell_file)
    con <- gzfile(out_path, "wt")
    tryCatch(
      write.table(merged, file = con, sep = "\t", quote = FALSE, col.names = NA),
      finally = close(con)
    )
  }
  invisible(cell_files)
}

replace_final_method_dir <- function(final_temp_dir, final_method_dir, backup_dir) {
  dir.create(dirname(final_method_dir), recursive = TRUE, showWarnings = FALSE)
  if (dir.exists(final_method_dir)) {
    if (dir.exists(backup_dir)) {
      unlink(backup_dir, recursive = TRUE, force = TRUE)
    }
    if (!file.rename(final_method_dir, backup_dir)) {
      stop("Could not move existing EPICunmix output aside: ", final_method_dir)
    }
  }

  ok <- file.rename(final_temp_dir, final_method_dir)
  if (!ok) {
    if (dir.exists(backup_dir) && !dir.exists(final_method_dir)) {
      file.rename(backup_dir, final_method_dir)
    }
    stop("Could not move merged EPICunmix output into final directory: ", final_method_dir)
  }
}

parser <- ArgumentParser(add_help = TRUE)
parser$add_argument("--dataset", type = "character", required = TRUE)
parser$add_argument("--config_id", type = "character", required = TRUE)
parser$add_argument("--n_core", type = "integer", required = FALSE, default = 15)
parser$add_argument("--extra_args", type = "character", required = FALSE, default = "")

args <- parser$parse_args()
repo_root <- find_repo_root(run_dir)
extra <- apply_method_default_extra_args(parse_extra_args(args$extra_args), "EPICunmix_with_posterior", repo_root = repo_root)
validate_extra_args(
  extra,
  c("map_cell_types", "use_test_samples", "use_limma_top_genes", "top_n", "chunk_size"),
  "EPICunmix_with_posterior"
)

dataset <- args$dataset
config_id <- args$config_id
n_core <- args$n_core
map_cell_types <- extra_arg(extra, "map_cell_types", TRUE)
use_test_samples <- extra_arg(extra, "use_test_samples", TRUE)
use_limma_top_genes <- extra_arg(extra, "use_limma_top_genes", FALSE)
top_n <- extra_arg(extra, "top_n", 100)
chunk_size <- extra_arg(extra, "chunk_size", 1000)
validate_top_n_arg(extra, use_limma_top_genes, "EPICunmix_with_posterior")
if (is.na(chunk_size) || chunk_size < 1) {
  stop("chunk_size must be a positive integer")
}

paths <- resolve_deconv_paths(dataset, config_id, repo_root = repo_root)
bulk_input <- paths$config$bulk_input[[1]]
bulk_scale <- paths$config$bulk_scale[[1]]
bulk_normalization <- paths$config$bulk_normalization[[1]]
frac_input <- paths$config$frac_input[[1]]
refType <- paths$refType
message_log_path <- method_message_log_path(paths, "EPICunmix")
cleanup_message_log <- start_message_log(message_log_path)
on.exit(cleanup_message_log(), add = TRUE)

bulk_expr <- read_bulk(paths)
bulk_prep <- prepare_bulk_for_deconv(bulk_expr, paths, method = "EPICunmix_with_posterior")
bulk_expr <- bulk_prep$bulk_expr
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

posterior_path <- file.path(paths$output_dir, "bMIND", "bMIND_posterior.RDS")
if (!file.exists(posterior_path)) {
  stop(
    "bMIND posterior not found for this config: ", posterior_path,
    ". Run run_bMIND.R for dataset=", dataset,
    " config_id=", config_id,
    " with extra_args='export_posterior=true' first."
  )
}
posterior <- readRDS(posterior_path)
posterior_genes <- posterior_gene_names(posterior)
posterior_cell_types <- posterior_cell_type_names(posterior)

aligned <- align_fraction_and_reference_cell_types(
  frac,
  posterior_cell_types,
  paths,
  repo_root
)
frac <- aligned$frac

log_bulk <- prepare_log_bulk_input(bulk_expr, bulk_input, config_id, method = "EPICunmix_with_posterior")
bulk_expr_logged <- log_bulk$bulk_expr

limma_gene_count <- NA_integer_
if (use_limma_top_genes) {
  selected_genes <- read_limma_top_gene_union(paths, top_n = top_n)
  limma_gene_count <- length(selected_genes)
} else {
  selected_genes <- NULL
}

genes <- intersect(rownames(bulk_expr_logged), posterior_genes)
if (!is.null(selected_genes)) {
  genes <- intersect(genes, selected_genes)
}
if (length(genes) == 0) {
  stop("No genes available for EPICunmix after intersecting bulk, posterior, and configured gene filters")
}

bulk_expr_logged <- bulk_expr_logged[genes, , drop = FALSE]
deconv_input_dim <- format_deconv_dim(bulk_expr_logged)
gene_chunks <- split_gene_chunks(genes, chunk_size)

if (!identical(colnames(frac), posterior_cell_types)) {
  stop("EPICunmix fraction cell types are not aligned to posterior cell types after subsetting")
}
if (!is.null(posterior$A) && !is.null(dimnames(posterior$A)[[3]])) {
  posterior_samples <- dimnames(posterior$A)[[3]]
  if (!identical(posterior_samples, colnames(bulk_expr_logged))) {
    stop("EPICunmix posterior samples are not aligned to bulk samples. Re-run bMIND posterior with matching use_test_samples/sample settings.")
  }
}

run_stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
gene_set_label <- if (use_limma_top_genes) paste0("top", top_n) else "all"
chunk_root <- file.path(
  paths$obj_dir,
  "logs",
  "EPICunmix_chunks",
  paste0(config_id, "_", gene_set_label, "_chunk", chunk_size, "_", run_stamp)
)
outputs_dir <- file.path(chunk_root, "outputs")
final_temp_dir <- file.path(chunk_root, "merged_EPICunmix")
final_method_dir <- file.path(paths$output_dir, "EPICunmix")
backup_dir <- file.path(chunk_root, "existing_EPICunmix_before_chunked_run")
dir.create(outputs_dir, recursive = TRUE, showWarnings = FALSE)

status <- data.frame(
  chunk = sprintf("chunk_%03d", seq_along(gene_chunks)),
  gene_count = vapply(gene_chunks, length, integer(1)),
  status = rep("pending", length(gene_chunks)),
  message = rep("", length(gene_chunks)),
  output_dir = file.path(outputs_dir, sprintf("chunk_%03d", seq_along(gene_chunks))),
  start_time = rep("", length(gene_chunks)),
  end_time = rep("", length(gene_chunks)),
  run_time = rep("", length(gene_chunks)),
  stringsAsFactors = FALSE
)
write_chunk_status(status, chunk_root)

message("Running EPICunmix with bMIND posterior")
message("  dataset: ", dataset)
message("  config_id: ", config_id)
message("  bulk_input: ", bulk_input, " (", bulk_input_dim, ")")
message("  bulk_scale: ", bulk_scale)
message("  bulk_normalization: ", bulk_normalization)
message("  bulk_preparation: ", bulk_prep$action)
message("  refType: ", refType)
message("  frac_input: ", frac_input, " (", nrow(frac), " samples x ", ncol(frac), " cell types)")
message("  use_test_samples: ", use_test_samples)
message("  n_core: ", n_core)
message("  use_limma_top_genes: ", use_limma_top_genes)
if (use_limma_top_genes) {
  message("  top_n: ", top_n)
  message("  limma_top_gene_union: ", limma_gene_count, " genes")
}
message("  chunk_size: ", chunk_size)
message("  n_chunks: ", length(gene_chunks))
message("  posterior_dim: ", length(posterior_genes), " genes x ", length(posterior_cell_types), " cell types")
message("  posterior: ", posterior_path)
message("  chunk_root: ", chunk_root)
message("  note: ", log_bulk$note)

message("Starting deconvolution")
message("  deconv_input: ", deconv_input_dim)

with_runtime_log(paths, "EPICunmix", n_core = n_core, deconv_input = deconv_input_dim, expr = {
  successful_output_dirs <- character(0)

  for (i in seq_along(gene_chunks)) {
    chunk_name <- status$chunk[[i]]
    chunk_genes <- gene_chunks[[i]]
    chunk_output_dir <- status$output_dir[[i]]
    chunk_start <- Sys.time()
    status$status[[i]] <- "running"
    status$start_time[[i]] <- format_runtime_timestamp(chunk_start)
    write_chunk_status(status, chunk_root)

    message("  chunk ", i, "/", length(gene_chunks), ": ", chunk_name, " (", length(chunk_genes), " genes)")
    chunk_ok <- tryCatch(
      {
        bulk_chunk <- bulk_expr_logged[chunk_genes, , drop = FALSE]
        posterior_chunk <- subset_posterior_genes(posterior, posterior_genes, chunk_genes)

        if (!is.null(posterior_chunk$A) && !identical(dimnames(posterior_chunk$A)[[1]], chunk_genes)) {
          stop("EPICunmix posterior gene names are not aligned for ", chunk_name)
        }

        epic_unmix <- EPICunmix::run_epic_unmix(
          bulk_chunk,
          frac,
          posterior_chunk,
          outf = FALSE,
          ncore = n_core
        )
        Z <- aperm(epic_unmix$A, c(1, 3, 2))
        dimnames(Z)[[2]] <- colnames(bulk_chunk)
        dimnames(Z)[[3]] <- colnames(frac)

        Z_export <- Z
        if (map_cell_types) {
          Z_export <- apply_indep_ref_cell_type_mapping(Z_export, paths, repo_root = repo_root)
        }
        write_ctse_array(Z_export, outputs_dir, chunk_name)

        rm(bulk_chunk, posterior_chunk, epic_unmix, Z, Z_export)
        TRUE
      },
      error = function(e) {
        status$message[[i]] <<- conditionMessage(e)
        FALSE
      }
    )

    chunk_end <- Sys.time()
    status$end_time[[i]] <- format_runtime_timestamp(chunk_end)
    status$run_time[[i]] <- format_elapsed_mins(chunk_start, chunk_end)
    if (chunk_ok) {
      status$status[[i]] <- "completed"
      successful_output_dirs <- c(successful_output_dirs, chunk_output_dir)
    } else {
      status$status[[i]] <- "failed"
      message("  FAILED ", chunk_name, ": ", status$message[[i]])
    }
    write_chunk_status(status, chunk_root)
    gc()
  }

  failed <- status$status != "completed"
  if (any(failed)) {
    stop(
      "EPICunmix failed for ", sum(failed), " of ", nrow(status),
      " chunk(s). Kept chunk diagnostics under: ", chunk_root
    )
  }

  message("Merging ", length(successful_output_dirs), " EPICunmix chunk(s)")
  merge_chunk_outputs(successful_output_dirs, final_temp_dir)
  message_existing_method_output(paths, "EPICunmix")
  replace_final_method_dir(final_temp_dir, final_method_dir, backup_dir)
})

message("Output saved")
message("  EPICunmix ctse: ", file.path(final_method_dir, "<cell_type>.txt.gz"))
if (map_cell_types && refType == "indep") {
  message("  ctse naming note: mapped cell types use benchmark names; unmapped reference cell types keep reference names")
} else {
  message("  ctse naming note: exported cell types use reference cell type names")
}
message("  cleaned_chunk_root: ", chunk_root)
unlink(chunk_root, recursive = TRUE, force = TRUE)
message("done")
