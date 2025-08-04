
cpm_normalization = function(scExpr,target_sum = 1e6){
  col_sums <- Matrix::colSums(scExpr)
  scaling_factors <- 1e6 / col_sums
  cpm <- scExpr %*% Matrix::Diagonal(x = scaling_factors)
  colnames(cpm) = colnames(scExpr)
  return(cpm)
}

get_mu_by_group = function(Expr,cell_type_labels){
  cell_type_labels = as.vector(cell_type_labels)
  stopifnot(ncol(Expr)==length(cell_type_labels))
  group = list()
  for(i in unique(cell_type_labels)){
    group[[i]] <- which(cell_type_labels %in% i)
  }
  C = do.call(cbind, lapply(group,function(x) Matrix::rowMeans(Expr[,x,drop=F])))
  C
}

get_reconstructed_X_from_weighted_Z = function(Z,F){
  stopifnot(all(colnames(F) %in% dimnames(Z)[[3]]))
  F = F[,dimnames(Z)[[3]]]
  
  k = dim(Z)[3]
  X_reconstructed <- matrix(0,nrow = dim(Z)[1],ncol = dim(Z)[2])
  
  for (i in 1:k) {
    Z_slice <- Z[,,i]  # Extract the i-th slice of Z
    F_slice <- F[,i] 
    X_reconstructed <- X_reconstructed + sweep(Z_slice,2,F_slice,'*')
  }
  return(X_reconstructed)
}

get_top_variable_genes <- function(expression_data, top_n = 1000) {
  
  log_transformed <- log2(expression_data + 1)
  gene_variances <- apply(log_transformed, 1, stats::var)
  top_genes <- names(sort(gene_variances, decreasing = TRUE))[1:top_n]
  
  return(top_genes)
}

indep_colnames_rename = function(mat,indep_mapping){
  # keep only matched cell types, and rename current colnames by source cell-type name
  
  stopifnot(all(indep_mapping$maxCorName %in% colnames(mat)))
  
  mat = mat[,indep_mapping$maxCorName,drop = F]
  colnames(mat) = indep_mapping$cell_type
  return(mat)
}


anno_by_list <- function(x, named_list) {
  vapply(x, function(item) {
    match <- names(named_list)[sapply(named_list, function(vals) item %in% vals)]
    if (length(match) > 0) match else item
  }, FUN.VALUE = character(1))
}

find_gene_specificity_group = function(gene,statistics_vec,group_type = 'top_n',
                                       n_vec = c(10,30,100,300,1000,3000,10000),
                                       cutoff_vec = c(0,1,2,3,4,5,6)){
  
  if (!(gene %in% names(statistics_vec))) {
    return("others")
  }
  statistics_vec = statistics_vec[order(statistics_vec,decreasing = T)]
  
  if(group_type == 'top_n'){
    gene_rank <- which(names(statistics_vec) == gene)
    group <- min(n_vec[n_vec >= gene_rank])
    return(paste0("top_", group))
  }else if (group_type == 'cutoff') {
    gene_value <- statistics_vec[gene]
    group <- max(cutoff_vec[cutoff_vec <= gene_value], na.rm = TRUE)
    return(paste0("cutoff_", group))
  } else {
    stop("Invalid group_type. Choose 'top_n' or 'cutoff'.")
  }
}


get_marker_list <- function(cellType_stats_matrix,select_method = 'top_n', n = 10, cutoff = 1){
  
  if(select_method == 'top_n'){
    if(n > nrow(cellType_stats_matrix)){
      warning('n exceeds the number of available rows in cellType_stats_matrix, will use available rows only!')
      n = nrow(cellType_stats_matrix)
    }
    top_list = lapply(colnames(cellType_stats_matrix), function(colname) {
      col <- cellType_stats_matrix[, colname]
      rownames(cellType_stats_matrix)[order(col, decreasing = TRUE)[1:n]]
    })
  }else if(select_method == 'cutoff'){
    
    top_list = lapply(colnames(cellType_stats_matrix), function(colname) {
      col <- setNames(cellType_stats_matrix[, colname],rownames(cellType_stats_matrix))
      selected_rows <- rownames(cellType_stats_matrix)[col >= cutoff]
      selected_rows <- selected_rows[order(col[selected_rows], decreasing = TRUE)]
      return(selected_rows)
    })
  }else{
    stop('invalid select_method provided. Available options are "top_n" and "cutoff"')
  }
  names(top_list) <- colnames(cellType_stats_matrix)    
  
  return(top_list)
}


compute_limma_statistics = function(pseudobulk_expr,cell_type_labels){
  stopifnot(ncol(pseudobulk_expr)==length(cell_type_labels))

  if(max(pseudobulk_expr)>100){
    warning('please make sure pseudobulk_expr is already in log scale (dismiss if this is for TCA Z-inferred)')
  }
  
  annotation=factor(cell_type_labels)
  design <- model.matrix(~0+annotation)
  colnames(design) <- unlist(lapply(strsplit(colnames(design),"annotation"), function(x) x[2]))
  cont.matrix <- matrix((-1/ncol(design)),nrow=ncol(design),ncol=ncol(design))
  colnames(cont.matrix) <- colnames(design)
  diag(cont.matrix) <- (ncol(design)-1)/ncol(design)
  
  fit <- limma::lmFit(pseudobulk_expr, design) 
  fit2 <- limma::contrasts.fit(fit, cont.matrix)
  fit2 <- limma::eBayes(fit2, trend=TRUE)
  
  limma_top_genes = fit2[["coefficients"]]
  return(limma_top_genes)
}

read.gmt<-function(gmt.file){
  g<-GSA::GSA.read.gmt(gmt.file)
  genelist<-g$genesets
  names(genelist)<-g$geneset.names
  rm(g)
  return(genelist)
}

linear_regress_out <- function(y,mod){
  
  y = as.matrix(y)
  
  # mod = model.matrix(~0+truthFrac %>% as.matrix())
  Hat = solve(crossprod(mod))%*%t(mod)
  
  # Multiplies the gene expression matrix (y) by the transposed hat matrix to compute the coefficients (or "weights") for each gene-cell type combination
  beta = y %*% t(Hat) 
  
  # Subtracts the explained portion from the original bulk expression to get the residuals, 
  # which are the parts of the expression not attributable to cell type fractions
  residuals = y - beta %*% t(mod) 
  
  # Add back the mean of each gene
  yr = residuals + rowMeans(y)  
  return(yr)
}

poisson_regress_out <- function(y, mod){
  # mod = model.matrix(~0+truthFrac %>% as.matrix())
  if (ncol(y) != nrow(mod)) {
    stop("The number of samples (columns of y) must match the number of rows in mod.")
  }
  # Initialize matrices for coefficients, residuals, and adjusted expression
  beta <- matrix(0, nrow = nrow(y), ncol = ncol(mod))  # Coefficients
  residuals <- matrix(0, nrow = nrow(y), ncol = ncol(y))  # Residuals
  yr <- matrix(0, nrow = nrow(y), ncol = ncol(y))  # Adjusted expression
  
  # Loop over each gene
  for (i in 1:nrow(y)) {
    # Fit Poisson regression for the current gene
    fit <- glm(y[i, ] ~ mod - 1, family = poisson(link = "log"))  # "-1" for no intercept
    beta[i, ] <- coef(fit) 
    fitted_values <- exp(mod %*% coef(fit))  
    residuals[i, ] <- y[i, ] - fitted_values
    yr[i, ] <- residuals[i, ] + mean(y[i, ])  
  }
  dimnames(yr) = dimnames(y)
  return(yr)
}

define_gene_specificity_groups <- function(gene_vec,
                                           statistics_vec, 
                                           n_vec = c(100, 300, 1000, 3000, 10000)) {
  ordered_genes <- names(sort(statistics_vec, decreasing = TRUE))
  group_map <- rep("all_genes", length(ordered_genes))
  names(group_map) <- ordered_genes
  
  for (n in n_vec) {
    group_map[ordered_genes[1:n][group_map[ordered_genes[1:n]] == "all_genes"]] <- paste0("top_", n)
  }
  result <- group_map[gene_vec]
  result[is.na(result)] <- "all_genes"
  return(result)
}

GeneSigTest_adapted = function(Bulk,frac,p_threshold = 0.05){
  writeLines("Using Linear Regression To Infer Probability of Expression...")
  p_tab <- NULL
  df <- (nrow(frac) - ncol(frac))
  for (i in 1:nrow(Bulk)) {
    exp <- Bulk[i, ]
    rlm.o <- lm(exp ~ as.matrix(frac) - 1)
    p <- pt(summary(rlm.o)$coef[, 3], df, lower.tail = FALSE)
    p_tab <- rbind(p_tab, p)
  }
  p_tab[is.nan(p_tab)] <- 1
  pvalue.m <- p_tab
  for (i in 1:ncol(pvalue.m)) pvalue.m[, i] <- p.adjust(pvalue.m[,i], method = "BH")
  colnames(pvalue.m) <- colnames(frac)
  rownames(pvalue.m) <- rownames(Bulk)
  
  score <- pvalue.m
  call <- matrix(0, nrow = nrow(score), ncol = ncol(score))
  rownames(call) <- rownames(score)
  colnames(call) <- colnames(score)
  for (ct in 1:ncol(score)) {
    gene <- rownames(score)[score[, ct] < p_threshold]
    call[gene, ct] <- 1
  }
  gene.call <- rowSums(call)
  
  res <- list(call = call, pval = score)
  return(res)
}

