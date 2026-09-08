pair_cor_pearson_fast <- function(a, b,
                                  margin = c("row", "column"),
                                  min_both = NULL,
                                  pairwise_na = TRUE,
                                  chunk_size = 1000,
                                  check_names = TRUE) {
  margin <- match.arg(margin)
  
  if (check_names) {
    stopifnot(
      isTRUE(all.equal(colnames(a), colnames(b))),
      isTRUE(all.equal(rownames(a), rownames(b)))
    )
  }
  
  stopifnot(identical(dim(a), dim(b)))
  
  a <- as.matrix(a)
  b <- as.matrix(b)
  
  storage.mode(a) <- "double"
  storage.mode(b) <- "double"
  
  finish_cor <- function(n, ssx, ssy, sxy) {
    out <- rep(NA_real_, length(n))
    
    good <- n >= 2 &
      is.finite(ssx) & is.finite(ssy) & is.finite(sxy) &
      ssx > 0 & ssy > 0

    if (any(good)) {
      good_idx <- which(good)
      candidate <- sxy[good_idx] /
        sqrt(ssx[good_idx] * ssy[good_idx])

      # A centered calculation can exceed the theoretical boundary only by
      # negligible floating-point error. Keep larger violations undefined.
      in_range <- is.finite(candidate) & abs(candidate) <= 1 + 1e-12
      out_idx <- good_idx[in_range]
      out[out_idx] <- pmax(-1, pmin(1, candidate[in_range]))
    }

    out
  }
  
  calc_cols <- function(x, y) {
    if (pairwise_na || !is.null(min_both)) {
      ok <- !is.na(x) & !is.na(y)
      
      if (!is.null(min_both)) {
        ok <- ok & x >= min_both & y >= min_both
      }
      
      xx <- x
      yy <- y
      xx[!ok] <- 0
      yy[!ok] <- 0
      
      n <- colSums(ok)

      xx <- sweep(xx, 2L, colSums(xx) / n, FUN = "-")
      yy <- sweep(yy, 2L, colSums(yy) / n, FUN = "-")
      xx[!ok] <- 0
      yy[!ok] <- 0
    } else {
      n <- rep.int(nrow(x), ncol(x))

      xx <- sweep(x, 2L, colMeans(x), FUN = "-")
      yy <- sweep(y, 2L, colMeans(y), FUN = "-")
    }

    finish_cor(
      n = n,
      ssx = colSums(xx * xx),
      ssy = colSums(yy * yy),
      sxy = colSums(xx * yy)
    )
  }
  
  calc_rows <- function(x, y) {
    if (pairwise_na || !is.null(min_both)) {
      ok <- !is.na(x) & !is.na(y)
      
      if (!is.null(min_both)) {
        ok <- ok & x >= min_both & y >= min_both
      }
      
      xx <- x
      yy <- y
      xx[!ok] <- 0
      yy[!ok] <- 0
      
      n <- rowSums(ok)

      xx <- sweep(xx, 1L, rowSums(xx) / n, FUN = "-")
      yy <- sweep(yy, 1L, rowSums(yy) / n, FUN = "-")
      xx[!ok] <- 0
      yy[!ok] <- 0
    } else {
      n <- rep.int(ncol(x), nrow(x))

      xx <- sweep(x, 1L, rowMeans(x), FUN = "-")
      yy <- sweep(y, 1L, rowMeans(y), FUN = "-")
    }

    finish_cor(
      n = n,
      ssx = rowSums(xx * xx),
      ssy = rowSums(yy * yy),
      sxy = rowSums(xx * yy)
    )
  }
  
  n_out <- if (margin == "column") ncol(a) else nrow(a)
  out <- rep(NA_real_, n_out)
  
  if (is.null(chunk_size)) {
    chunk_size <- n_out
  }
  
  starts <- seq.int(1L, n_out, by = chunk_size)
  
  for (lo in starts) {
    hi <- min(lo + chunk_size - 1L, n_out)
    ii <- lo:hi
    
    out[ii] <- if (margin == "column") {
      calc_cols(a[, ii, drop = FALSE], b[, ii, drop = FALSE])
    } else {
      calc_rows(a[ii, , drop = FALSE], b[ii, , drop = FALSE])
    }
  }
  
  names(out) <- if (margin == "column") colnames(a) else rownames(a)
  out
}

pair_cor_fast <- function(a, b,
                          margin   = c("row", "column"),
                          method   = c("pearson", "spearman"),
                          min_both = NULL) {
  margin <- match.arg(margin)
  method <- match.arg(tolower(method), c("pearson", "spearman"))
  
  # Always check names
  stopifnot(
    isTRUE(all.equal(colnames(a), colnames(b))),
    isTRUE(all.equal(rownames(a), rownames(b)))
  )
  
  stopifnot(identical(dim(a), dim(b)))
  
  # Keep sparse Matrix objects sparse when possible.
  is_matrix_like <- function(x) {
    is.matrix(x) || inherits(x, "Matrix")
  }
  
  if (!is_matrix_like(a)) a <- as.matrix(a)
  if (!is_matrix_like(b)) b <- as.matrix(b)
  
  nr <- nrow(a)
  nc <- ncol(a)
  
  n_out <- if (margin == "column") nc else nr
  vec_len <- if (margin == "column") nr else nc
  
  if (n_out == 0L) {
    return(numeric(0))
  }
  
  # Internal adaptive chunk size.
  # Increase target_cells if you have lots of RAM.
  # Decrease it if memory spikes.
  target_cells <- 5e6
  chunk_size <- max(1L, min(n_out, floor(target_cells / max(1L, vec_len))))
  
  dense_double <- function(x) {
    x <- as.matrix(x)
    storage.mode(x) <- "double"
    x
  }
  
  finish_cor <- function(n, sx, sy, sxx, syy, sxy) {
    vx <- n * sxx - sx * sx
    vy <- n * syy - sy * sy
    
    out <- rep(NA_real_, length(n))
    
    good <- n >= 2L &
      is.finite(vx) & is.finite(vy) &
      vx > 0 & vy > 0
    
    out[good] <- (n[good] * sxy[good] - sx[good] * sy[good]) /
      sqrt(vx[good] * vy[good])
    
    out[!is.finite(out)] <- NA_real_
    out
  }
  
  pearson_cols <- function(x, y, keep = NULL) {
    if (!is.null(keep)) {
      x <- dense_double(x)
      y <- dense_double(y)
      
      x[!keep] <- 0
      y[!keep] <- 0
      
      n <- colSums(keep)
    } else {
      n <- rep.int(nrow(x), ncol(x))
    }
    
    finish_cor(
      n   = as.numeric(n),
      sx  = as.numeric(colSums(x)),
      sy  = as.numeric(colSums(y)),
      sxx = as.numeric(colSums(x * x)),
      syy = as.numeric(colSums(y * y)),
      sxy = as.numeric(colSums(x * y))
    )
  }
  
  rank_cols_average <- function(x) {
    if (requireNamespace("matrixStats", quietly = TRUE)) {
      r <- matrixStats::colRanks(
        x,
        ties.method   = "average",
        preserveShape = TRUE
      )
      
      storage.mode(r) <- "double"
      
      # Keep NA positions as NA.
      r[is.na(x)] <- NA_real_
      
      r
    } else {
      r <- matrix(NA_real_, nrow = nrow(x), ncol = ncol(x))
      
      for (j in seq_len(ncol(x))) {
        r[, j] <- rank(
          x[, j],
          ties.method = "average",
          na.last = "keep"
        )
      }
      
      r
    }
  }
  
  make_keep <- function(x, y) {
    keep <- !is.na(x) & !is.na(y)
    
    if (!is.null(min_both)) {
      keep <- keep & x >= min_both & y >= min_both
    }
    
    keep
  }
  
  cor_cols <- function(x, y) {
    if (method == "pearson") {
      pair_cor_pearson_fast(
        x,
        y,
        margin = "column",
        min_both = min_both,
        pairwise_na = TRUE,
        chunk_size = NULL,
        check_names = FALSE
      )
    } else {
      # Spearman = Pearson correlation of ranks.
      # Always pairwise-complete: remove bad paired positions before ranking.
      x <- dense_double(x)
      y <- dense_double(y)
      
      keep <- make_keep(x, y)
      
      x[!keep] <- NA_real_
      y[!keep] <- NA_real_
      
      rx <- rank_cols_average(x)
      ry <- rank_cols_average(y)
      
      if (all(keep)) {
        pearson_cols(rx, ry, keep = NULL)
      } else {
        pearson_cols(rx, ry, keep = keep)
      }
    }
  }
  
  out <- rep(NA_real_, n_out)
  
  starts <- seq.int(1L, n_out, by = chunk_size)
  
  # Always show progress
  pb <- utils::txtProgressBar(
    min = 0,
    max = length(starts),
    style = 3
  )
  on.exit(close(pb), add = TRUE)
  
  for (kk in seq_along(starts)) {
    lo <- starts[kk]
    hi <- min(lo + chunk_size - 1L, n_out)
    ii <- lo:hi
    
    if (margin == "column") {
      x <- a[, ii, drop = FALSE]
      y <- b[, ii, drop = FALSE]
    } else {
      # Convert row-wise problem into column-wise problem.
      # Each original row becomes one column after transpose.
      x <- t(a[ii, , drop = FALSE])
      y <- t(b[ii, , drop = FALSE])
    }
    
    out[ii] <- cor_cols(x, y)
    
    utils::setTxtProgressBar(pb, kk)
  }
  
  names(out) <- if (margin == "column") colnames(a) else rownames(a)
  
  out
}



############# CTSE performance helpers #############

read_ctse_matrix <- function(path, label = "CTSE matrix") {
  if (!file.exists(path)) {
    stop(label, " not found: ", path)
  }
  x <- read.delim(path, sep = "\t", check.names = FALSE, row.names = 1)
  x <- as.matrix(x)
  storage.mode(x) <- "numeric"
  x[is.na(x)] <- 0
  x
}

write_metric_matrix <- function(x, path, digits = 2) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  if (file.exists(path)) {
    warning("Overwriting existing metric file: ", path)
  }

  x <- as.matrix(x)
  storage.mode(x) <- "numeric"
  idx <- is.finite(x)
  x[idx] <- round(x[idx], digits)

  write.table(
    x,
    file = path,
    sep = "\t",
    quote = FALSE,
    col.names = NA
  )
}

ctse_cell_type_files <- function(dir) {
  if (!dir.exists(dir)) {
    stop("CTSE directory not found: ", dir)
  }
  files <- list.files(dir, pattern = "\\.txt\\.gz$", full.names = TRUE)
  names(files) <- sub("\\.txt\\.gz$", "", basename(files))
  files[order(names(files))]
}

apply_ctse_transform <- function(x, transform = c("none", "log2p1")) {
  transform <- match.arg(transform)
  if (transform == "none") {
    return(x)
  }
  if (any(x < 0, na.rm = TRUE)) {
    stop("Cannot apply log2p1 transform to matrix with negative values")
  }
  log2(x + 1)
}

compute_ctse_metric_rows <- function(truth,
                                     estimate,
                                     metric = c("spearman_cor", "pearson_cor"),
                                     samples = NULL,
                                     min_n_sample = 10) {
  metric <- match.arg(metric)

  genes <- intersect(rownames(truth), rownames(estimate))
  common_samples <- intersect(colnames(truth), colnames(estimate))
  if (!is.null(samples)) {
    common_samples <- intersect(common_samples, samples)
  }

  if (length(genes) == 0) {
    stop("No common genes between truth and estimate")
  }

  truth <- truth[genes, common_samples, drop = FALSE]
  estimate <- estimate[genes, common_samples, drop = FALSE]

  if (length(common_samples) < min_n_sample) {
    out <- rep(NA_real_, length(genes))
    names(out) <- genes
    return(out)
  }

  if (metric == "spearman_cor") {
    return(pair_cor_fast(truth, estimate, margin = "row", method = "spearman"))
  }
  if (metric == "pearson_cor") {
    return(pair_cor_fast(truth, estimate, margin = "row", method = "pearson"))
  }

  stop("Unsupported metric: ", metric)
}

merge_metric_vectors <- function(metric_vectors) {
  if (length(metric_vectors) == 0) {
    stop("No metric vectors to merge")
  }

  genes <- sort(unique(unlist(lapply(metric_vectors, names), use.names = FALSE)))
  out <- matrix(
    NA_real_,
    nrow = length(genes),
    ncol = length(metric_vectors),
    dimnames = list(genes, names(metric_vectors))
  )

  for (cell_type in names(metric_vectors)) {
    values <- metric_vectors[[cell_type]]
    out[names(values), cell_type] <- values
  }

  out
}

read_gmt_gene_union <- function(path) {
  if (!file.exists(path)) {
    stop("GMT file not found: ", path)
  }

  lines <- readLines(path, warn = FALSE)
  lines <- lines[nzchar(trimws(lines))]
  if (length(lines) == 0) {
    stop("GMT file contains no gene sets: ", path)
  }

  fields <- strsplit(lines, "\t", fixed = TRUE)
  malformed <- which(lengths(fields) < 3L)
  if (length(malformed) > 0) {
    stop(
      "GMT lines must contain a set name, description, and at least one gene; ",
      "malformed line(s): ", paste(malformed, collapse = ", ")
    )
  }

  genes <- unique(unlist(lapply(fields, function(x) x[-c(1L, 2L)]), use.names = FALSE))
  genes <- trimws(genes)
  genes[nzchar(genes)]
}

compute_ctse_sample_cor <- function(truth, estimate, genes = NULL, samples = NULL) {
  common_genes <- intersect(rownames(truth), rownames(estimate))
  if (!is.null(genes)) {
    common_genes <- intersect(common_genes, genes)
  }

  common_samples <- intersect(colnames(truth), colnames(estimate))
  if (!is.null(samples)) {
    common_samples <- intersect(common_samples, samples)
  }

  if (length(common_samples) == 0L) {
    return(stats::setNames(numeric(), character()))
  }

  # Correlation is mathematically undefined with fewer than two paired genes.
  if (length(common_genes) < 2L) {
    out <- rep(NA_real_, length(common_samples))
    names(out) <- common_samples
    return(out)
  }

  truth <- truth[common_genes, common_samples, drop = FALSE]
  estimate <- estimate[common_genes, common_samples, drop = FALSE]
  pair_cor_fast(truth, estimate, margin = "column", method = "spearman")
}

merge_sample_metric_vectors <- function(metric_vectors) {
  if (length(metric_vectors) == 0) {
    stop("No sample metric vectors to merge")
  }

  samples <- unique(unlist(lapply(metric_vectors, names), use.names = FALSE))
  out <- matrix(
    NA_real_,
    nrow = length(samples),
    ncol = length(metric_vectors),
    dimnames = list(samples, names(metric_vectors))
  )

  for (cell_type in names(metric_vectors)) {
    values <- metric_vectors[[cell_type]]
    out[names(values), cell_type] <- values
  }

  out
}

