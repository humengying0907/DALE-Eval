############################################################
# General pseudobulk and bulk-noise helpers
#
# This module is intended for benchmark-wide reuse. It contains
# sparse-safe matrix helpers, CPM/log transforms, CTSE/pseudobulk
# construction helpers, and multiplicative bulk noise simulation.
############################################################

is_sparse_matrix <- function(x) {
  inherits(x, "sparseMatrix")
}

require_matrix_package <- function() {
  if (!requireNamespace("Matrix", quietly = TRUE)) {
    stop("The Matrix package is required for sparse Matrix input")
  }
  invisible(TRUE)
}

matrix_col_sums <- function(x, na.rm = TRUE) {
  if (is_sparse_matrix(x)) {
    require_matrix_package()
    return(Matrix::colSums(x, na.rm = na.rm))
  }
  colSums(x, na.rm = na.rm)
}

matrix_row_sums <- function(x, na.rm = TRUE) {
  if (is_sparse_matrix(x)) {
    require_matrix_package()
    return(Matrix::rowSums(x, na.rm = na.rm))
  }
  rowSums(x, na.rm = na.rm)
}

matrix_row_means <- function(x, na.rm = TRUE) {
  if (is_sparse_matrix(x)) {
    require_matrix_package()
    return(Matrix::rowMeans(x, na.rm = na.rm))
  }
  rowMeans(x, na.rm = na.rm)
}

renorm_cpm <- function(mat, scale = 1e6) {
  if (is_sparse_matrix(mat)) {
    require_matrix_package()
    denom <- matrix_col_sums(mat, na.rm = TRUE)
    factors <- numeric(ncol(mat))
    good <- is.finite(denom) & denom > 0
    factors[good] <- scale / denom[good]
    out <- mat %*% Matrix::Diagonal(x = factors)
    dimnames(out) <- dimnames(mat)
    return(out)
  }

  mat <- as.matrix(mat)
  denom <- colSums(mat, na.rm = TRUE)
  out <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  good <- is.finite(denom) & denom > 0
  if (any(good)) {
    out[, good] <- sweep(mat[, good, drop = FALSE], 2, denom[good], "/") * scale
  }
  out
}

counts_to_cpm <- function(counts, scale = 1e6) {
  renorm_cpm(counts, scale = scale)
}

log2_cpm_plus1 <- function(mat, input_scale = c("counts", "cpm")) {
  input_scale <- match.arg(input_scale)
  if (input_scale == "counts") {
    mat <- counts_to_cpm(mat)
  }
  log2(as.matrix(mat) + 1)
}

check_counts_meta <- function(counts, meta, sample_col = "sample_id", celltype_col = "cell_type") {
  if (ncol(counts) != nrow(meta)) {
    stop("ncol(counts) must equal nrow(meta)")
  }
  if (!all(colnames(counts) == rownames(meta))) {
    stop("colnames(counts) must match rownames(meta) in the same order")
  }
  missing_cols <- setdiff(c(sample_col, celltype_col), colnames(meta))
  if (length(missing_cols) > 0) {
    stop("meta is missing columns: ", paste(missing_cols, collapse = ", "))
  }
  invisible(TRUE)
}

make_cellfrac <- function(meta, sample_col = "sample_id", celltype_col = "cell_type") {
  samples <- unique(as.character(meta[[sample_col]]))
  cell_types <- unique(as.character(meta[[celltype_col]]))
  tab <- table(
    factor(as.character(meta[[sample_col]]), levels = samples),
    factor(as.character(meta[[celltype_col]]), levels = cell_types)
  )
  frac <- sweep(tab, 1, rowSums(tab), "/")
  frac[!is.finite(frac)] <- 0
  as.matrix(frac)
}

make_ctse_sumcount <- function(counts, meta, sample_col = "sample_id", celltype_col = "cell_type") {
  samples <- unique(as.character(meta[[sample_col]]))
  cell_types <- unique(as.character(meta[[celltype_col]]))

  out <- setNames(vector("list", length(cell_types)), cell_types)
  for (ct in cell_types) {
    mat <- matrix(0, nrow = nrow(counts), ncol = length(samples), dimnames = list(rownames(counts), samples))
    for (s in samples) {
      idx <- which(as.character(meta[[sample_col]]) == s & as.character(meta[[celltype_col]]) == ct)
      if (length(idx) > 0) {
        mat[, s] <- matrix_row_sums(counts[, idx, drop = FALSE])
      }
    }
    out[[ct]] <- mat
  }
  out
}

make_ctse_meancpm <- function(counts, meta, sample_col = "sample_id", celltype_col = "cell_type") {
  cell_cpm <- counts_to_cpm(counts)
  samples <- unique(as.character(meta[[sample_col]]))
  cell_types <- unique(as.character(meta[[celltype_col]]))

  out <- setNames(vector("list", length(cell_types)), cell_types)
  for (ct in cell_types) {
    mat <- matrix(0, nrow = nrow(counts), ncol = length(samples), dimnames = list(rownames(counts), samples))
    for (s in samples) {
      idx <- which(as.character(meta[[sample_col]]) == s & as.character(meta[[celltype_col]]) == ct)
      if (length(idx) > 0) {
        mat[, s] <- matrix_row_means(cell_cpm[, idx, drop = FALSE])
      }
    }
    out[[ct]] <- mat
  }
  out
}

mix_ctse_by_fraction <- function(ctse_by_ct, cellfrac) {
  genes <- rownames(ctse_by_ct[[1]])
  samples <- colnames(ctse_by_ct[[1]])
  bulk <- matrix(0, nrow = length(genes), ncol = length(samples), dimnames = list(genes, samples))

  for (ct in names(ctse_by_ct)) {
    weights <- if (ct %in% colnames(cellfrac)) cellfrac[samples, ct] else rep(0, length(samples))
    bulk <- bulk + sweep(ctse_by_ct[[ct]], 2, weights, "*")
  }
  bulk
}

sum_ctse_counts <- function(ctse_sumcount) {
  Reduce("+", ctse_sumcount)
}

sample_multinomial_counts <- function(prob_mat, library_sizes) {
  prob_mat <- as.matrix(prob_mat)
  if (length(library_sizes) == 1) {
    library_sizes <- rep(library_sizes, ncol(prob_mat))
  }
  stopifnot(length(library_sizes) == ncol(prob_mat))

  out <- matrix(0, nrow = nrow(prob_mat), ncol = ncol(prob_mat), dimnames = dimnames(prob_mat))
  for (j in seq_len(ncol(prob_mat))) {
    p <- prob_mat[, j]
    p[!is.finite(p) | p < 0] <- 0
    if (sum(p) == 0) {
      p <- rep(1 / length(p), length(p))
    }
    p <- p / sum(p)
    out[, j] <- as.vector(stats::rmultinom(1, size = round(library_sizes[j]), prob = p))
  }
  out
}

# Explicit multiplicative bulk-noise helpers.
# epsilon[g,i] ~ Normal(mu_shift_log2[g], sigma^2), multiplier[g,i] = 2^epsilon[g,i].
#
# simulate_cpm_direct_multiplicative_noise():
#   CPM -> multiply by 2^epsilon -> optional CPM renormalization.
#
# simulate_cpm_probability_multinomial_noise():
#   CPM -> probabilities -> multiply probabilities by 2^epsilon
#   -> renormalize probabilities -> multinomial sampling.
#
# simulate_count_probability_multinomial_noise():
#   counts -> probabilities -> multiply probabilities by 2^epsilon
#   -> renormalize probabilities -> multinomial sampling.

make_log2_noise_components <- function(clean_mat,
                                       mu_shift_log2 = 0,
                                       sigma = 0.5) {
  sigma <- as.numeric(sigma)
  if (length(sigma) != 1 || !is.finite(sigma) || sigma < 0) {
    stop("sigma must be one finite non-negative number")
  }
  clean_mat <- as.matrix(clean_mat)

  mu_shift_input <- mu_shift_log2
  if (!is.null(names(mu_shift_input)) && all(rownames(clean_mat) %in% names(mu_shift_input))) {
    mu_shift_log2 <- as.numeric(mu_shift_input[rownames(clean_mat)])
  } else {
    mu_shift_log2 <- as.numeric(mu_shift_input)
  }

  if (length(mu_shift_log2) == 1) {
    mu_shift_log2 <- rep(mu_shift_log2, nrow(clean_mat))
  }
  if (length(mu_shift_log2) != nrow(clean_mat) || any(!is.finite(mu_shift_log2))) {
    stop("mu_shift_log2 must be one finite number or one finite value per gene")
  }
  names(mu_shift_log2) <- rownames(clean_mat)

  sigma_log2 <- rep(sigma, nrow(clean_mat))
  names(sigma_log2) <- rownames(clean_mat)

  epsilon_log2 <- matrix(
    stats::rnorm(
      length(clean_mat),
      mean = rep(mu_shift_log2, times = ncol(clean_mat)),
      sd = sigma
    ),
    nrow = nrow(clean_mat),
    ncol = ncol(clean_mat),
    dimnames = dimnames(clean_mat)
  )

  list(
    mu_shift_log2 = mu_shift_log2,
    sigma_log2 = sigma_log2,
    sigma = sigma,
    epsilon_log2 = epsilon_log2,
    multiplier = 2^epsilon_log2
  )
}

simulate_cpm_direct_multiplicative_noise <- function(clean_cpm,
                                                     mu_shift_log2 = 0,
                                                     sigma = 0.5,
                                                     renormalize_cpm = TRUE) {
  clean_cpm <- as.matrix(clean_cpm)
  clean_cpm[!is.finite(clean_cpm) | clean_cpm < 0] <- 0

  noise <- make_log2_noise_components(clean_cpm, mu_shift_log2 = mu_shift_log2, sigma = sigma)
  bulk_raw <- clean_cpm * noise$multiplier
  bulk <- if (isTRUE(renormalize_cpm)) renorm_cpm(bulk_raw) else bulk_raw

  list(
    bulk = bulk,
    noise = c(
      noise,
      list(
        mode = "cpm_direct",
        bulk_raw = bulk_raw,
        renormalize_cpm = renormalize_cpm,
        output_scale = if (isTRUE(renormalize_cpm)) "cpm" else "cpm_raw",
        colsum_before_renorm = matrix_col_sums(bulk_raw),
        colsum_after_renorm = matrix_col_sums(bulk)
      )
    )
  )
}

simulate_cpm_probability_multinomial_noise <- function(clean_cpm,
                                                       mu_shift_log2 = 0,
                                                       sigma = 0.5,
                                                       library_sizes = NULL,
                                                       output_scale = c("cpm", "counts")) {
  output_scale <- match.arg(output_scale)
  clean_cpm <- as.matrix(clean_cpm)
  clean_cpm[!is.finite(clean_cpm) | clean_cpm < 0] <- 0

  if (is.null(library_sizes)) {
    library_sizes <- matrix_col_sums(clean_cpm)
  }
  if (length(library_sizes) == 1) {
    library_sizes <- rep(library_sizes, ncol(clean_cpm))
  }
  stopifnot(length(library_sizes) == ncol(clean_cpm))
  names(library_sizes) <- colnames(clean_cpm)

  noise <- make_log2_noise_components(clean_cpm, mu_shift_log2 = mu_shift_log2, sigma = sigma)
  probability <- as.matrix(renorm_cpm(clean_cpm, scale = 1))
  noisy_probability_raw <- probability * noise$multiplier
  noisy_probability <- as.matrix(renorm_cpm(noisy_probability_raw, scale = 1))
  expected_count <- sweep(noisy_probability, 2, library_sizes, "*")
  sampled_count <- sample_multinomial_counts(noisy_probability, library_sizes)
  bulk <- if (output_scale == "counts") sampled_count else counts_to_cpm(sampled_count)

  list(
    bulk = bulk,
    noise = c(
      noise,
      list(
        mode = "cpm_probability_multinomial",
        probability = probability,
        noisy_probability_raw = noisy_probability_raw,
        noisy_probability = noisy_probability,
        library_sizes = library_sizes,
        expected_count = expected_count,
        sampled_count = sampled_count,
        output_scale = output_scale,
        realized_library_sizes = matrix_col_sums(sampled_count)
      )
    )
  )
}

simulate_count_probability_multinomial_noise <- function(clean_counts,
                                                         mu_shift_log2 = 0,
                                                         sigma = 0.5,
                                                         library_sizes = NULL,
                                                         count_sampler = c("multinomial", "none"),
                                                         output_scale = c("counts", "cpm")) {
  count_sampler <- match.arg(count_sampler)
  output_scale <- match.arg(output_scale)
  clean_counts <- as.matrix(clean_counts)
  clean_counts[!is.finite(clean_counts) | clean_counts < 0] <- 0

  if (is.null(library_sizes)) {
    library_sizes <- matrix_col_sums(clean_counts)
  }
  if (length(library_sizes) == 1) {
    library_sizes <- rep(library_sizes, ncol(clean_counts))
  }
  stopifnot(length(library_sizes) == ncol(clean_counts))
  names(library_sizes) <- colnames(clean_counts)

  noise <- make_log2_noise_components(clean_counts, mu_shift_log2 = mu_shift_log2, sigma = sigma)
  probability <- as.matrix(renorm_cpm(clean_counts, scale = 1))
  noisy_probability_raw <- probability * noise$multiplier
  noisy_probability <- as.matrix(renorm_cpm(noisy_probability_raw, scale = 1))
  expected_count <- sweep(noisy_probability, 2, library_sizes, "*")
  bulk_counts <- if (count_sampler == "multinomial") {
    sample_multinomial_counts(noisy_probability, library_sizes)
  } else {
    expected_count
  }
  bulk <- if (output_scale == "cpm") counts_to_cpm(bulk_counts) else bulk_counts

  list(
    bulk = bulk,
    noise = c(
      noise,
      list(
        mode = "count_probability_multinomial",
        probability = probability,
        noisy_probability_raw = noisy_probability_raw,
        noisy_probability = noisy_probability,
        library_sizes = library_sizes,
        expected_count = expected_count,
        sampled_count = bulk_counts,
        count_sampler = count_sampler,
        output_scale = output_scale,
        realized_library_sizes = matrix_col_sums(bulk_counts)
      )
    )
  )
}

# Backward-compatible alias for direct CPM perturbation.
simulate_cpm_multiplicative_noise <- function(clean_cpm,
                                              mu_shift_log2 = 0,
                                              sigma = 0.5,
                                              renormalize_cpm = TRUE) {
  simulate_cpm_direct_multiplicative_noise(
    clean_cpm = clean_cpm,
    mu_shift_log2 = mu_shift_log2,
    sigma = sigma,
    renormalize_cpm = renormalize_cpm
  )
}


# Sparse gene-level shift and quick QC helpers for noisy pseudobulk simulations.

make_sparse_expression_mu_shift <- function(pseudobulk_cpm,
                                            candidate_expression_threshold_cpm = 1,
                                            positive_fraction = 0.15,
                                            positive_mu_shift_mean = 0.4,
                                            positive_mu_shift_sd = 0.15,
                                            expressed_background_mu_shift = 0,
                                            low_mu_shift_mean = -0.1,
                                            low_mu_shift_sd = 0.05,
                                            seed = NULL) {
  if (!is.null(seed)) {
    set.seed(seed)
  }

  pseudobulk_cpm <- as.matrix(pseudobulk_cpm)
  params <- c(
    candidate_expression_threshold_cpm = candidate_expression_threshold_cpm,
    positive_fraction = positive_fraction,
    positive_mu_shift_mean = positive_mu_shift_mean,
    positive_mu_shift_sd = positive_mu_shift_sd,
    expressed_background_mu_shift = expressed_background_mu_shift,
    low_mu_shift_mean = low_mu_shift_mean,
    low_mu_shift_sd = low_mu_shift_sd
  )
  if (any(!is.finite(params))) {
    stop("All sparse mu-shift parameters must be finite")
  }
  if (candidate_expression_threshold_cpm < 0) {
    stop("candidate_expression_threshold_cpm must be non-negative")
  }
  if (positive_fraction < 0 || positive_fraction > 1) {
    stop("positive_fraction must be between 0 and 1")
  }
  if (positive_mu_shift_sd < 0 || low_mu_shift_sd < 0) {
    stop("positive_mu_shift_sd and low_mu_shift_sd must be non-negative")
  }

  genes <- rownames(pseudobulk_cpm)
  mean_cpm <- rowMeans(pseudobulk_cpm, na.rm = TRUE)
  expressed <- is.finite(mean_cpm) & mean_cpm >= candidate_expression_threshold_cpm
  low_expression <- !expressed
  expressed_genes <- genes[expressed]

  mu_shift_log2 <- rep(expressed_background_mu_shift, length(genes))
  names(mu_shift_log2) <- genes

  if (any(low_expression)) {
    mu_shift_log2[low_expression] <- stats::rnorm(
      n = sum(low_expression),
      mean = low_mu_shift_mean,
      sd = low_mu_shift_sd
    )
  }

  n_expressed <- length(expressed_genes)
  if (n_expressed == 0) {
    return(mu_shift_log2)
  }

  n_positive <- min(round(n_expressed * positive_fraction), n_expressed)
  positive_genes <- if (n_positive > 0) sample(expressed_genes, n_positive) else character()

  if (length(positive_genes) > 0) {
    mu_shift_log2[positive_genes] <- stats::rnorm(
      n = length(positive_genes),
      mean = positive_mu_shift_mean,
      sd = positive_mu_shift_sd
    )
  }

  mu_shift_log2
}

plot_pseudobulk_noise_qc_grid <- function(clean_bulk,
                                          simulated_obj,
                                          clean_scale = c("cpm", "counts"),
                                          label = "noisy pseudobulk vs clean pseudobulk",
                                          bins = 80,
                                          trim_quantiles = c(0.01, 0.99),
                                          point_size = 0.35,
                                          expressed_threshold_cpm = 1) {
  if (!requireNamespace("ggplot2", quietly = TRUE) ||
      !requireNamespace("gridExtra", quietly = TRUE)) {
    stop("Packages ggplot2 and gridExtra are required")
  }

  clean_scale <- match.arg(clean_scale)
  noisy_bulk <- if (is.list(simulated_obj) && "bulk" %in% names(simulated_obj)) {
    simulated_obj$bulk
  } else {
    simulated_obj
  }

  clean_cpm <- if (clean_scale == "counts") counts_to_cpm(clean_bulk) else as.matrix(clean_bulk)
  noisy_cpm <- as.matrix(noisy_bulk)
  if (clean_scale == "counts") {
    noisy_cpm <- counts_to_cpm(noisy_cpm)
  }

  common_residual_genes <- intersect(rownames(clean_cpm), rownames(noisy_cpm))
  common_residual_samples <- intersect(colnames(clean_cpm), colnames(noisy_cpm))
  clean_log <- log2(clean_cpm[common_residual_genes, common_residual_samples, drop = FALSE] + 1)
  noisy_log <- log2(noisy_cpm[common_residual_genes, common_residual_samples, drop = FALSE] + 1)
  residual <- noisy_log - clean_log
  residual_vec <- as.vector(residual)
  residual_vec <- residual_vec[is.finite(residual_vec)]
  mu_g <- rowMeans(residual, na.rm = TRUE)

  residual_df <- data.frame(
    residual = residual_vec,
    stringsAsFactors = FALSE
  )
  residual_df <- residual_df[is.finite(residual_df$residual), , drop = FALSE]
  residual_xlim <- stats::quantile(residual_df$residual, probs = trim_quantiles, na.rm = TRUE)

  p_residual <- ggplot2::ggplot(residual_df, ggplot2::aes(x = residual)) +
    ggplot2::geom_histogram(bins = bins, fill = "gray60", color = "white") +
    ggplot2::geom_vline(xintercept = 0, color = "red", linetype = "dashed") +
    ggplot2::coord_cartesian(xlim = residual_xlim) +
    ggplot2::theme_classic() +
    ggplot2::labs(
      x = "log2(noisy CPM + 1) - log2(clean CPM + 1)",
      y = "Gene-sample count",
      title = "Residual"
    )

  clean_mean_cpm <- rowMeans(clean_cpm, na.rm = TRUE)
  mu_df <- data.frame(
    gene = names(mu_g),
    mu_g = mu_g,
    clean_mean_cpm = clean_mean_cpm[names(mu_g)],
    stringsAsFactors = FALSE
  )
  mu_df <- mu_df[is.finite(mu_df$mu_g), , drop = FALSE]
  mu_xlim <- stats::quantile(mu_df$mu_g, probs = trim_quantiles, na.rm = TRUE)

  p_mu <- ggplot2::ggplot(mu_df, ggplot2::aes(x = mu_g)) +
    ggplot2::geom_histogram(bins = bins, fill = "gray60", color = "white") +
    ggplot2::geom_vline(xintercept = 0, color = "red", linetype = "dashed") +
    ggplot2::coord_cartesian(xlim = mu_xlim) +
    ggplot2::theme_classic() +
    ggplot2::labs(
      x = "Gene-wise mean residual (mu_g)",
      y = "Gene count",
      title = "mu_g: all genes"
    )

  mu_expressed_df <- mu_df[
    is.finite(mu_df$clean_mean_cpm) &
      mu_df$clean_mean_cpm >= expressed_threshold_cpm,
    ,
    drop = FALSE
  ]
  if (nrow(mu_expressed_df) == 0) {
    mu_expressed_df <- data.frame(mu_g = numeric())
  }
  mu_expressed_xlim <- if (nrow(mu_expressed_df) > 0) {
    stats::quantile(mu_expressed_df$mu_g, probs = trim_quantiles, na.rm = TRUE)
  } else {
    mu_xlim
  }

  p_mu_expressed <- ggplot2::ggplot(mu_expressed_df, ggplot2::aes(x = mu_g)) +
    ggplot2::geom_histogram(bins = bins, fill = "gray60", color = "white") +
    ggplot2::geom_vline(xintercept = 0, color = "red", linetype = "dashed") +
    ggplot2::coord_cartesian(xlim = mu_expressed_xlim) +
    ggplot2::theme_classic() +
    ggplot2::labs(
      x = "Gene-wise mean residual (mu_g)",
      y = "Gene count",
      title = paste0("mu_g: clean mean CPM >= ", expressed_threshold_cpm)
    )

  common_genes <- intersect(rownames(clean_cpm), rownames(noisy_cpm))
  mean_df <- data.frame(
    clean_mean = rowMeans(clean_cpm[common_genes, , drop = FALSE], na.rm = TRUE),
    noisy_mean = rowMeans(noisy_cpm[common_genes, , drop = FALSE], na.rm = TRUE),
    stringsAsFactors = FALSE
  )
  mean_df <- mean_df[is.finite(mean_df$clean_mean) & is.finite(mean_df$noisy_mean), , drop = FALSE]

  p_mean <- ggplot2::ggplot(mean_df, ggplot2::aes(x = log2(clean_mean + 1), y = log2(noisy_mean + 1)))
  if (requireNamespace("ggpointdensity", quietly = TRUE) &&
      requireNamespace("viridis", quietly = TRUE)) {
    p_mean <- p_mean +
      ggpointdensity::geom_pointdensity(size = point_size) +
      viridis::scale_color_viridis(option = "magma") +
      ggplot2::labs(color = "Density")
  } else {
    p_mean <- p_mean +
      ggplot2::geom_point(size = point_size, alpha = 0.35, color = "gray30")
  }
  p_mean <- p_mean +
    ggplot2::geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
    ggplot2::theme_classic() +
    ggplot2::labs(
      x = "Clean mean log2(CPM + 1)",
      y = "Noisy mean log2(CPM + 1)",
      title = "Mean comparison"
    )

  gridExtra::grid.arrange(p_residual, p_mu, p_mu_expressed, p_mean, ncol = 4, top = label)
}


make_gene_mean_variance_df <- function(pseudobulk,
                                       real_bulk_cpm_by_protocol,
                                       pseudobulk_scale = c("cpm", "counts"),
                                       protocols = names(real_bulk_cpm_by_protocol),
                                       genes = NULL) {
  pseudobulk_scale <- match.arg(pseudobulk_scale)
  pseudobulk_log <- log2_cpm_plus1(pseudobulk, input_scale = pseudobulk_scale)

  out <- lapply(protocols, function(protocol) {
    if (!protocol %in% names(real_bulk_cpm_by_protocol)) {
      stop("real_bulk_cpm_by_protocol is missing protocol: ", protocol)
    }

    real_log <- log2_cpm_plus1(real_bulk_cpm_by_protocol[[protocol]], input_scale = "cpm")
    common_genes <- intersect(rownames(pseudobulk_log), rownames(real_log))
    if (!is.null(genes)) {
      common_genes <- intersect(common_genes, genes)
    }
    common_samples <- intersect(colnames(pseudobulk_log), colnames(real_log))

    pb <- pseudobulk_log[common_genes, common_samples, drop = FALSE]
    rb <- real_log[common_genes, common_samples, drop = FALSE]

    data.frame(
      protocol = protocol,
      gene = common_genes,
      pseudobulk_mean = rowMeans(pb, na.rm = TRUE),
      real_bulk_mean = rowMeans(rb, na.rm = TRUE),
      pseudobulk_variance = apply(pb, 1, stats::var, na.rm = TRUE),
      real_bulk_variance = apply(rb, 1, stats::var, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  })

  out <- do.call(rbind, out)
  out[
    is.finite(out$pseudobulk_mean) &
      is.finite(out$real_bulk_mean) &
      is.finite(out$pseudobulk_variance) &
      is.finite(out$real_bulk_variance),
    ,
    drop = FALSE
  ]
}

plot_gene_mean_variance_diagnostics <- function(mean_variance_df,
                                                title = "Gene-wise mean/variance diagnostics",
                                                point_size = 0.35,
                                                trim_quantiles = c(0.01, 0.99)) {
  if (!requireNamespace("ggplot2", quietly = TRUE) ||
      !requireNamespace("gridExtra", quietly = TRUE)) {
    stop("Packages ggplot2 and gridExtra are required")
  }

  df <- mean_variance_df
  mean_xlim <- stats::quantile(df$pseudobulk_mean, trim_quantiles, na.rm = TRUE)
  mean_ylim <- stats::quantile(df$real_bulk_mean, trim_quantiles, na.rm = TRUE)
  var_xlim <- stats::quantile(log2(df$pseudobulk_variance + 1), trim_quantiles, na.rm = TRUE)
  var_ylim <- stats::quantile(log2(df$real_bulk_variance + 1), trim_quantiles, na.rm = TRUE)

  point_layer <- function() {
    if (requireNamespace("ggpointdensity", quietly = TRUE) &&
        requireNamespace("viridis", quietly = TRUE)) {
      list(
        ggpointdensity::geom_pointdensity(size = point_size),
        viridis::scale_color_viridis(option = "magma"),
        ggplot2::labs(color = "Density")
      )
    } else {
      ggplot2::geom_point(size = point_size, alpha = 0.35, color = "gray30")
    }
  }

  p_mean <- ggplot2::ggplot(df, ggplot2::aes(x = pseudobulk_mean, y = real_bulk_mean)) +
    point_layer() +
    ggplot2::geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
    ggplot2::coord_cartesian(xlim = mean_xlim, ylim = mean_ylim) +
    ggplot2::facet_wrap(~ protocol) +
    ggplot2::theme_classic() +
    ggplot2::labs(
      x = "Pseudobulk gene mean: mean log2(CPM + 1)",
      y = "Real bulk gene mean: mean log2(CPM + 1)",
      title = "Gene-wise mean"
    )

  p_variance <- ggplot2::ggplot(
    df,
    ggplot2::aes(x = log2(pseudobulk_variance + 1), y = log2(real_bulk_variance + 1))
  ) +
    point_layer() +
    ggplot2::geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
    ggplot2::coord_cartesian(xlim = var_xlim, ylim = var_ylim) +
    ggplot2::facet_wrap(~ protocol) +
    ggplot2::theme_classic() +
    ggplot2::labs(
      x = "Pseudobulk gene variance: log2(var + 1)",
      y = "Real bulk gene variance: log2(var + 1)",
      title = "Gene-wise variance"
    )

  mv_pb <- data.frame(
    protocol = df$protocol,
    source = "pseudobulk",
    mean = df$pseudobulk_mean,
    variance = df$pseudobulk_variance,
    stringsAsFactors = FALSE
  )
  mv_real <- data.frame(
    protocol = df$protocol,
    source = "real_bulk",
    mean = df$real_bulk_mean,
    variance = df$real_bulk_variance,
    stringsAsFactors = FALSE
  )
  mv_df <- rbind(mv_pb, mv_real)
  mv_df <- mv_df[is.finite(mv_df$mean) & is.finite(mv_df$variance), , drop = FALSE]

  p_mean_variance <- ggplot2::ggplot(
    mv_df,
    ggplot2::aes(x = mean, y = log2(variance + 1), color = source)
  ) +
    ggplot2::geom_point(size = point_size, alpha = 0.35) +
    ggplot2::facet_wrap(~ protocol) +
    ggplot2::theme_classic() +
    ggplot2::labs(
      x = "Gene mean: mean log2(CPM + 1)",
      y = "Gene variance: log2(var + 1)",
      color = "Source",
      title = "Mean-variance structure"
    )

  gridExtra::grid.arrange(p_mean, p_variance, p_mean_variance, ncol = 3, top = title)
}

