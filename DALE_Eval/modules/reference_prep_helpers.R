cpm_normalization <- function(scExpr, target_sum = 1e6) {
  col_sums <- Matrix::colSums(scExpr)
  # Preserve legacy behavior from archived_folder/scripts/helpers.R.
  scaling_factors <- 1e6 / col_sums
  cpm <- scExpr %*% Matrix::Diagonal(x = scaling_factors)
  colnames(cpm) <- colnames(scExpr)
  cpm
}

get_mu_by_group <- function(Expr, cell_type_labels) {
  cell_type_labels <- as.vector(cell_type_labels)
  stopifnot(ncol(Expr) == length(cell_type_labels))

  group <- list()
  for (i in unique(cell_type_labels)) {
    group[[i]] <- which(cell_type_labels %in% i)
  }
  do.call(cbind, lapply(group, function(x) Matrix::rowMeans(Expr[, x, drop = FALSE])))
}

get_top_variable_genes <- function(expression_data, top_n = 1000) {
  log_transformed <- log2(expression_data + 1)
  gene_variances <- apply(log_transformed, 1, stats::var)
  names(sort(gene_variances, decreasing = TRUE))[seq_len(min(top_n, length(gene_variances)))]
}

compute_limma_statistics <- function(pseudobulk_expr, cell_type_labels) {
  stopifnot(ncol(pseudobulk_expr) == length(cell_type_labels))

  if (max(pseudobulk_expr) > 100) {
    warning("please make sure pseudobulk_expr is already in log scale (dismiss if this is for TCA Z-inferred)")
  }

  annotation <- factor(cell_type_labels)
  design <- model.matrix(~ 0 + annotation)
  colnames(design) <- unlist(lapply(strsplit(colnames(design), "annotation"), function(x) x[2]))

  cont_matrix <- matrix((-1 / ncol(design)), nrow = ncol(design), ncol = ncol(design))
  colnames(cont_matrix) <- colnames(design)
  diag(cont_matrix) <- (ncol(design) - 1) / ncol(design)

  fit <- limma::lmFit(pseudobulk_expr, design)
  fit2 <- limma::contrasts.fit(fit, cont_matrix)
  fit2 <- limma::eBayes(fit2, trend = TRUE)
  fit2[["coefficients"]]
}

ct_recode <- function(cell_types, mapping_list) {
  matched <- sapply(cell_types, function(ct) {
    matched_name <- names(mapping_list)[sapply(mapping_list, function(x) ct %in% x)]
    if (length(matched_name) == 0) {
      NA
    } else {
      matched_name
    }
  })
  unname(unlist(matched))
}
