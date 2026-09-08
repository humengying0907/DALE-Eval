limma_top_gene_union <- function(limma_stats, top_n = 100) {
  if (!is.matrix(limma_stats)) {
    limma_stats <- as.matrix(limma_stats)
  }
  storage.mode(limma_stats) <- "numeric"
  if (is.null(rownames(limma_stats))) {
    stop("limma_top_genes matrix must have gene row names")
  }
  if (is.null(colnames(limma_stats))) {
    stop("limma_top_genes matrix must have cell-type columns")
  }
  if (!is.numeric(top_n) || length(top_n) != 1 || is.na(top_n) || top_n < 1) {
    stop("top_n must be a positive number")
  }

  top_n <- as.integer(top_n)
  genes <- unlist(
    lapply(seq_len(ncol(limma_stats)), function(j) {
      stats <- limma_stats[, j]
      keep <- is.finite(stats)
      if (!any(keep)) {
        return(character(0))
      }
      ranked <- order(stats[keep], decreasing = TRUE)
      rownames(limma_stats)[keep][ranked][seq_len(min(top_n, length(ranked)))]
    }),
    use.names = FALSE
  )
  unique(genes)
}


read_limma_top_gene_union <- function(paths, top_n = 100) {
  limma_stats <- read_reference_matrix(paths$ref_dir, c("limma_top_genes.csv", "limma_top_genes.txt"), "limma top genes")
  limma_top_gene_union(limma_stats, top_n = top_n)
}


read_gene_list_file <- function(path) {
  if (is.null(path) || is.na(path) || !nzchar(path)) {
    stop("gene_list_path must be a non-empty path")
  }
  if (!file.exists(path)) {
    stop("gene_list_path not found: ", path)
  }

  lines <- readLines(path, warn = FALSE)
  lines <- trimws(lines)
  lines <- lines[nzchar(lines)]
  lines <- lines[!startsWith(lines, "#")]
  genes <- vapply(strsplit(lines, "[[:space:]]+"), `[`, character(1), 1)
  genes <- unique(genes[nzchar(genes)])
  if (length(genes) == 0) {
    stop("gene_list_path contains no genes: ", path)
  }
  genes
}
