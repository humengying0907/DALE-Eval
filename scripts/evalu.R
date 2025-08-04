
########### functions for correlation calculation ############
pair_cor = function(a,b,margin = 'row',cor_method = 'pearson'){
  # a: gene * sample
  # b: gene * sample
  
  # make sure a,b has the same rows and columns
  stopifnot(all.equal(colnames(a),colnames(b)))
  stopifnot(all.equal(rownames(a),rownames(b)))
  
  corr = c()
  
  if(margin == 'row'){
    for(i in 1:nrow(a)){
      corr[i] = cor(a[i,],b[i,],method = cor_method,use = 'pairwise.complete.obs')
    }
  }else if(margin == 'column'){
    for(i in 1:ncol(a)){
      corr[i] = cor(a[,i],b[,i],method = cor_method,use = 'pairwise.complete.obs')
    }
  }
  
  return(corr)
}

get_per_gene_cor = function(Z_inferred,
                            Z_truth,
                            truthFrac,
                            margin = 'row',
                            cor_method = 'spearman',
                            min_frac = 0.001,
                            min_nSample = 10){
  
  # Ensure the cell types in Z_inferred match those in Z_truth
  stopifnot(all(dimnames(Z_inferred)[[3]] %in% dimnames(Z_truth)[[3]]))
  stopifnot(all(colnames(truthFrac)%in% dimnames(Z_truth)[[3]]))
  
  # for matched cell_types only
  ct = intersect(dimnames(Z_inferred)[[3]],dimnames(Z_truth)[[3]])
  ct = ct[order(ct)]
  Z_inferred = Z_inferred[,,ct]
  Z_truth = Z_truth[,,ct]
  
  # keep genes with no missing data in Z_inferred
  g_with_no_na <- apply(Z_inferred, 1, function(x) !anyNA(x))
  Z_inferred <- Z_inferred[g_with_no_na,,]
  
  # Subset Z_truth and Z_inferred to common genes
  cms <- intersect(dimnames(Z_truth)[[1]], dimnames(Z_inferred)[[1]])
  Z_truth <- Z_truth[cms,,]
  Z_inferred <- Z_inferred[cms,,]
  
  cell_types <- dimnames(Z_inferred)[[3]]
  truthFrac <- truthFrac[,cell_types]
  
  # Initialize result matrix
  if (margin == 'row') {
    result_matrix <- matrix(NA, nrow = length(cms), ncol = length(cell_types), dimnames = list(cms, cell_types))
  } else if (margin == 'column') {
    result_matrix <- matrix(NA, nrow = dim(Z_inferred)[2], ncol = length(cell_types), dimnames = list(dimnames(Z_inferred)[[2]], cell_types))
  }
  
  # Compute correlations for each cell type
  for (ct in cell_types) {
    # only consider samples with frac > min_frac for cor calculation 
    id = which(truthFrac[,ct]>min_frac)
    
    if(margin == 'row'){
      if(length(id)<min_nSample){
        warning(paste('less than',min_nSample,'samples available for correlation calculation in',ct,'\n will return NA for gene correlation'))
        correlations = rep(NA,length(cms))
      }else{
        correlations <- pair_cor(Z_truth[,id,ct], Z_inferred[,id,ct], margin = margin, cor_method = cor_method)
      }
      result_matrix[,ct] = correlations
      
    }else if(margin == 'column'){
      correlations <- pair_cor(Z_truth[,id,ct], Z_inferred[,id,ct], margin = margin, cor_method = cor_method)
      result_matrix[id,ct] = correlations
    }
  }
  
  return(as.data.frame(result_matrix))
}

########### functions for RMSE calculation #############
compute_nRMSE <- function(true, pred, normalization = 'by_mean') {
  rmse <- function(true, pred) {
    sqrt(mean((true - pred)^2, na.rm = TRUE))
  }
  
  if (normalization == 'by_mean') {
    mean_val <- mean(true, na.rm = TRUE)
    if (is.na(mean_val) || mean_val == 0) return(NA)
    return(rmse(true, pred) / mean_val)
    
  } else if (normalization == 'by_range') {
    range_val <- max(true, na.rm = TRUE) - min(true, na.rm = TRUE)
    if (is.na(range_val) || range_val == 0) return(NA)
    return(rmse(true, pred) / range_val)
    
  } else {
    stop('Please provide a valid normalization method: "by_mean" or "by_range"')
  }
}


pair_nRMSE = function(a, b, margin = 'row',normalization = 'by_mean'){
  # a: gene * sample (true)
  # b: gene * sample (pred)
  
  # make sure a,b has the same rows and columns
  stopifnot(all.equal(colnames(a),colnames(b)))
  stopifnot(all.equal(rownames(a),rownames(b)))
  
  RMSE = c()
  
  if(margin == 'row'){
    for(i in 1:nrow(a)){
      RMSE[i] = compute_nRMSE(a[i,],b[i,],normalization)
    }
  }else if(margin == 'column'){
    for(i in 1:ncol(a)){
      RMSE[i] = compute_nRMSE(a[,i],b[,i],normalization)
    }
  }
  return(RMSE)
}

get_per_gene_nRMSE = function(Z_inferred,
                              Z_truth,
                              margin = 'row',
                              normalization = 'by_mean'){
  # by default, we will only consider RMSE in log scale
  if(max(Z_truth)>20){
    warning('Z_truth seems not in log scale, please note we recommend nRMSE calculation for log scale only; dismiss this if this is for RMSE calculation on shifted values')
  }
  
  # Ensure the cell types in Z_inferred match those in Z_truth
  stopifnot(all(dimnames(Z_inferred)[[3]] %in% dimnames(Z_truth)[[3]]))
  
  # for matched cell_types only
  ct = intersect(dimnames(Z_inferred)[[3]],dimnames(Z_truth)[[3]])
  ct = ct[order(ct)]
  Z_inferred = Z_inferred[,,ct]
  Z_truth = Z_truth[,,ct]
  
  # keep genes with no missing data in Z_inferred
  g_with_no_na <- apply(Z_inferred, 1, function(x) !anyNA(x))
  Z_inferred <- Z_inferred[g_with_no_na,,]
  
  # Subset Z_truth and Z_inferred to common genes
  cms <- intersect(dimnames(Z_truth)[[1]], dimnames(Z_inferred)[[1]])
  Z_truth <- Z_truth[cms,,]
  Z_inferred <- Z_inferred[cms,,]
  
  # Initialize result matrix
  cell_types <- dimnames(Z_inferred)[[3]]
  if (margin == 'row') {
    result_matrix <- matrix(NA, nrow = length(cms), ncol = length(cell_types), dimnames = list(cms, cell_types))
  } else if (margin == 'column') {
    result_matrix <- matrix(NA, nrow = dim(Z_inferred)[2], ncol = length(cell_types), dimnames = list(dimnames(Z_inferred)[[2]], cell_types))
  }
  
  for (ct in cell_types){
    result_matrix[,ct] = pair_nRMSE(Z_truth[,,ct], Z_inferred[,,ct], margin, normalization)
  }
  return(as.data.frame(result_matrix))
}

########### functions related to Z dimension changes #########
get_3d_bulk = function(bulk_expr,cell_types){
  bulk_expr = as.matrix(bulk_expr)
  Z = array(NA,dim = c(nrow(bulk_expr),ncol(bulk_expr),length(cell_types)))
  
  for(i in 1:length(cell_types)){
    Z[,,i] = bulk_expr
  }
  
  dimnames(Z)[[1]] = rownames(bulk_expr)
  dimnames(Z)[[2]] = colnames(bulk_expr)
  dimnames(Z)[[3]] = cell_types
  return(Z)
}

get_Z_flattened = function(Z){
  g = dim(Z)[[1]]
  n = dim(Z)[[2]]
  k = dim(Z)[[3]]
  
  Z_flattened = matrix(NA,nrow = g,ncol = n*k)
  
  column_names = c()
  for(i in 1:k){
    Z_flattened[,(1+n*(i-1)):(i*n)] = Z[,,i]
    column_names = c(column_names,paste0(dimnames(Z)[[3]][i],'_',dimnames(Z)[[2]]))
  }
  rownames(Z_flattened) = dimnames(Z)[[1]]
  colnames(Z_flattened) = column_names
  
  return(Z_flattened)
}

get_group_levels_from_Z_flattened = function(Z_flattened,unique_sampleIDs){
  pattern = paste0("_(", paste(unique_sampleIDs, collapse = "|"), ")")
  group_levels = gsub(pattern, "", colnames(Z_flattened))
  return(group_levels)
}


############# functions related to Z intrinsic statistics (gene by cell type features) ############
compute_Z_mean = function(Z){
  Z_flattened = get_Z_flattened(Z)
  group_levels = get_group_levels_from_Z_flattened(Z_flattened,dimnames(Z)[[2]])
  mean_df = as.data.frame(get_mu_by_group(Z_flattened,group_levels))
  return(mean_df)
}

get_Z_limma_statistics = function(Z,scale = 'linear'){
  Z_flattened = get_Z_flattened(Z)
  ct_labels = get_group_levels_from_Z_flattened(Z_flattened,dimnames(Z)[[2]])
  
  if(scale == 'log2'){
    # no log-transformation needed, will directly run limma fit on Z_flattened
    limma_top_genes = compute_limma_statistics(Z_flattened,ct_labels)
  }else if(scale == 'linear'){
    # compensate for negative values (applicable to Unico only) before log transformation
    if(min(Z_flattened) < 0){
      shifts <- pmax(0, -matrixStats::rowMins(Z_flattened))
      Z_flattened = Z_flattened + shifts  
    }
    limma_top_genes = compute_limma_statistics(log2(Z_flattened + 1),ct_labels)
  }
  
  return(limma_top_genes)
}

compute_Z_mean_fc = function(Z_means,scale = 'linear'){
  
  fc_top_genes <- sapply(1:ncol(Z_means), function(j) {
    other_cols <- Z_means[, -j, drop = FALSE]
    max_other_cols <- apply(other_cols, 1, max)
    
    if(scale == 'linear'){
      row_min <- apply(Z_means, 1, min, na.rm = TRUE)
      shift_value <- ifelse(row_min < 0, -row_min, 0)
      Z_means_shifted <- Z_means + shift_value 
      
      log2(pmax(Z_means_shifted[, j], 1e-10) / pmax(max_other_cols + shift_value, 1e-10))
    } else if(scale == 'log2'){
      Z_means[, j] - max_other_cols
      
    }
  })
  
  colnames(fc_top_genes) <- colnames(Z_means)
  return(fc_top_genes)
}

compute_partial_r2 = function(Z,total_expr){
  
  require(rsq)
  
  total_expr = as.matrix(total_expr)
  total_expr = total_expr[rowSums(total_expr > 0) > 1 ,] # keep genes expressed in at least 2 samples
  
  cms = intersect(rownames(total_expr),dimnames(Z)[[1]])
  
  cell_types = dimnames(Z)[[3]]
  
  nGenes <- length(cms)
  if(nGenes<10){
    stop('too few intersect genes between Z and total_expr, please check gene names')
  }
  
  partial_r2_df = data.frame(matrix(NA,nrow = length(cms),ncol = length(cell_types)))
  rownames(partial_r2_df) = cms
  colnames(partial_r2_df) = cell_types
  
  checkpoint <- floor(seq(0.1, 1, by = 0.1) * nGenes)
  
  for (i in seq_along(cms)) {
    g <- cms[i]
    
    for (cell_type in cell_types) {
      assign(cell_type, Z[g, , cell_type])
    }
    
    formula <- as.formula(paste("total_expr[g,] ~", paste(cell_types, collapse = " + ")))
    model <- lm(formula)
    partial_r2 <- rsq.partial(model)
    partial_r2_df[g,] <- round(partial_r2$partial.rsq, 3)
    
    # Print progress at every 10%
    if (i %in% checkpoint) {
      message(paste("Partial_r2 Calculation Progress:", round(i / nGenes * 100), "% (", i, "/", nGenes, "genes)", sep = " "))
    }
  }
  
  return(partial_r2_df)
}

######### functions related to specificity calculation #######
# get covariance structure for the g*k*k covariance values from Z
get_offDiagCor = function(Z, cor_method = 'pearson'){
  offDiagCor = double(dim(Z)[1])
  
  for(i in 1:dim(Z)[[1]]){
    cmat=cor(Z[i,,], Z[i,,],method = cor_method)
    offDiagCor[i]=mean(cmat[upper.tri(cmat)],na.rm = T)
  }
  names(offDiagCor) = dimnames(Z)[[1]]
  
  return(offDiagCor)
}

# for covariance calculation, we only consider Pearson correlation!
get_variation_specificity = function(Z,cellType_stats_matrix,indep_mapping = NULL,top_n = 10000){
  # by default, will use Z-truth derived limma statistics as cellType_stats_matrix!
  # for covariance calculation, we only consider Pearson correlation!
  marker_list = get_marker_list(cellType_stats_matrix,select_method = 'top_n',n = top_n)
  marker_list <- lapply(marker_list, function(gene_list) {
    gene_list[gene_list %in% rownames(Z)]
  })
  
  if(is.null(indep_mapping)){
    stopifnot(all(dimnames(Z)[[3]] %in% names(marker_list)))
  }else{
    stopifnot(all(indep_mapping$maxCorName %in% dimnames(Z)[[3]]))
    Z = Z[,,indep_mapping$maxCorName]
    dimnames(Z)[[3]] = indep_mapping$cell_type
  }
  
  var_specificity_res = list()
  
  compute_var_specificity <- function(cell_type) {
    if (length(marker_list[[cell_type]]) == 0) {
      return(NULL) 
    }
    
    # please note, for cbsx, which can have the same value for all other cell-types, it will return NA 
    f <- do.call(rbind, lapply(marker_list[[cell_type]], function(x)
      cor(Z[x, , cell_type], Z[x, , ], use = 'pairwise.complete.obs'))) %>% as.data.frame()
    
    other_cts = setdiff(names(f), cell_type)

    f$other_ct_avg_cor <- rowMeans(f[, other_cts,drop = F],na.rm = T)
    f$other_ct_max_cor <- apply(f[, other_cts,drop = F],1, max, na.rm = T)
    
    f$varDissimilarity_avg = 1 - abs(f$other_ct_avg_cor)
    f$varDissimilarity_max = 1 - abs(f$other_ct_max_cor)
    
    f$gene <- marker_list[[cell_type]]
    f$specificity_ct <- cell_type
    
    statistics_vec <- setNames(cellType_stats_matrix[, cell_type], rownames(cellType_stats_matrix))
    
    f$specificity_level1 <- sapply(f$gene, find_gene_specificity_group,
                                   statistics_vec = statistics_vec, group_type = "top_n")
    f$specificity_level2 <- sapply(f$gene, find_gene_specificity_group, 
                                   statistics_vec = statistics_vec, group_type = "cutoff")
    return(f)
  }
  
  var_specificity_table = do.call(rbind,lapply(dimnames(Z)[[3]],compute_var_specificity))
  rownames(var_specificity_table) = NULL
  return(var_specificity_table)
}


get_expr_specificity = function(Z,cellType_stats_matrix,indep_mapping = NULL, top_n = 10000, 
                                Z_scale = 'linear'){
  # please note that the expression specificity (Z_means) will be present in log scales!
  # the interpretation of Z-means will be slightly different for Z in linear scale vs log2 scale
  # for Z in linear scales, the difference between the columns (log(mean_ct1) - log(mean_ct2)) can be directly interpreted as logFC
  # but for Z in log2 scales, each column represents mean(log_ct1), which can not be considered as log fold change (but mean differences in log space)
  # considering the majority of the methods we tested are in linear scale, we will not apply further adjustment for Z in log scale (such as unlog)
  
  # please make sure cellType_stats_matrix only contains cell-types from Z-truth! not indep reference!
  marker_list = get_marker_list(cellType_stats_matrix,select_method = 'top_n',n = top_n)
  marker_list <- lapply(marker_list, function(gene_list) {
    gene_list[gene_list %in% rownames(Z)]
  })
  
  Z_means_raw = compute_Z_mean(Z)
  
  if(Z_scale == 'linear'){
    row_min <- apply(Z_means_raw, 1, min, na.rm = TRUE)
    shift_value <- ifelse(row_min < 0, -row_min, 0)
    Z_means_shifted <- Z_means_raw + shift_value 
    
    Z_means = log2(Z_means_shifted + 1)
  }else if(Z_scale == 'log2'){
    Z_means = Z_means_raw
  }
  
  if(is.null(indep_mapping)){
    stopifnot(all(dimnames(Z)[[3]] %in% names(marker_list)))
  }else{
    stopifnot(all(indep_mapping$maxCorName %in% dimnames(Z)[[3]]))
    Z_means = Z_means[,indep_mapping$maxCorName]
    colnames(Z_means) = indep_mapping$cell_type
  }
  
  expr_specificity_res = list()
  
  compute_expr_specificity = function(cell_type){
    if (length(marker_list[[cell_type]]) == 0) {
      return(NULL) 
    }
    
    f = Z_means[marker_list[[cell_type]],]
    other_cts = setdiff(names(f), cell_type)
    
    other_ct_avg_expr = rowMeans(f[, other_cts],na.rm = T)
    other_ct_max_expr = apply(f[, other_cts, drop = FALSE], 1, max,na.rm = T)
    
    f$other_ct_avg_expr <- other_ct_avg_expr
    f$other_ct_max_expr <- other_ct_max_expr
    
    f$logFC_avg <- f[, cell_type] - f$other_ct_avg_expr
    f$logFC_max <- f[, cell_type] - f$other_ct_max_expr
    
    f$gene <- marker_list[[cell_type]]
    f$specificity_ct <- cell_type
    
    statistics_vec <- setNames(cellType_stats_matrix[, cell_type], rownames(cellType_stats_matrix))
    
    f$specificity_level1 <- sapply(f$gene, find_gene_specificity_group, 
                                   statistics_vec = statistics_vec, group_type = "top_n")
    f$specificity_level2 <- sapply(f$gene, find_gene_specificity_group, 
                                   statistics_vec = statistics_vec, group_type = "cutoff")
    
    
    return(f)
  } 
  
  expr_specificity_table = do.call(rbind,lapply(colnames(Z_means),compute_expr_specificity))
  rownames(expr_specificity_table) = NULL
  return(expr_specificity_table)
}


########## functions for cor_res_list processing ############
add_missing_columns <- function(df, full_colnames) {
  missing_cols <- setdiff(full_colnames, names(df))
  for (col in missing_cols) {
    df[[col]] <- NA
  }
  df <- df[, full_colnames]
  return(df)
}

check_same_columns_in_list <- function(df_list) {
  first_df_columns <- colnames(df_list[[1]])
  same_columns <- all(sapply(df_list, function(x) identical(colnames(x), first_df_columns)))
  return(same_columns)
}

cor_by_ct_summary = function(list_of_matrices){
  require(dplyr)
  require(tidyverse)
  
  full_col_method =  names(which.max(sapply(list_of_matrices, ncol)))
  list_of_matrices <- lapply(list_of_matrices, add_missing_columns, 
                             full_colnames = colnames(list_of_matrices[[full_col_method]]))
  stopifnot(check_same_columns_in_list(list_of_matrices))
  
  # Convert each matrix to a data frame and add row names as a column
  list_of_dfs <- lapply(list_of_matrices, function(x) {
    df <- as.data.frame(x)  # Convert matrix to data frame
    df$gene <- row.names(x)  # Add row names as a column
    return(df)
  })
  
  result_list <- list()
  column_names <- colnames(list_of_dfs[[1]])[-length(colnames(list_of_dfs[[1]]))]  # Exclude 'gene'
  
  # Column-wise binding
  for (col_name in column_names) {
    temp_list <- list()
    
    for (i in seq_along(list_of_dfs)) {
      temp_df <- list_of_dfs[[i]][, c("gene", col_name), drop = FALSE]  # Ensure it's a data frame
      colnames(temp_df) <- c("gene", names(list_of_matrices)[i])  # Rename columns correctly
      temp_list[[i]] <- temp_df
    }
    
    # Perform a full join by row names for all data frames in the temp list
    full_joined_df <- Reduce(function(x, y) full_join(x, y, by = "gene"), temp_list)
    
    # Add the fully joined data frame to the result list
    result_list[[col_name]] <- full_joined_df
  }
  return(result_list)
}

quick_summary = function(cor_res_list,marker_list = NULL){
  Z_summary = cor_by_ct_summary(cor_res_list)
  
  if(!is.null(marker_list)){
    
    stopifnot(all(names(marker_list) %in% names(Z_summary)))
    Z_summary = Z_summary[names(marker_list)]
    
    for(ct in names(marker_list)){
      df = Z_summary[[ct]]
      df = df[df$gene %in% marker_list[[ct]],]
      Z_summary[[ct]] = df
    }
  }
  
  Z_summary_df <- do.call(rbind, lapply(names(Z_summary), function(name) {
    Z_summary[[name]] %>%
      gather(key = "method", value = "cor", -1) %>% 
      mutate(cell_type = name)
  }))
  Z_summary_df = Z_summary_df[!is.na(Z_summary_df$cor),]
  return(Z_summary_df)
}

get_by_specificity_summary = function(cor_res_list, 
                                      cellType_stats_matrix, # gene by cell_type matrix, filled with numeric values showing specificity
                                      select_method = 'top_n',
                                      n_vec = c(10,30,100,300,1000,3000,10000),
                                      cutoff_vec = c(0,1,2,3,4),
                                      add_all_genes = F){

  # a unified top gene list will be applied to all method
  if(select_method == 'top_n'){
    cor_by_specificity <- do.call(rbind, lapply(n_vec, function(x) {
      marker_list <- get_marker_list(cellType_stats_matrix,'top_n',x) 
      quick_summary(cor_res_list, marker_list) %>%     
        mutate(group = paste0('top_',x))                               
    }))    

  }else if(select_method == 'cutoff'){
    cor_by_specificity <- do.call(rbind, lapply(cutoff_vec, function(x) {
      marker_list <- get_marker_list(cellType_stats_matrix,'cutoff',NULL,x) 
      quick_summary(cor_res_list, marker_list) %>%     
        mutate(group = paste0('cutoff_',x))                               
    }))    
  }
  
  if(add_all_genes){
    cor_by_specificity = rbind(cor_by_specificity,quick_summary(cor_res_list, NULL) %>%     
                                 mutate(group = 'all_genes') )
  }
  
  cor_by_specificity_summary <- cor_by_specificity %>%
    group_by(method, cell_type, group) %>%
    summarise(
      avg_cor = mean(cor, na.rm = TRUE), 
      q25 = quantile(cor, 0.25, na.rm = TRUE), 
      q75 = quantile(cor, 0.75, na.rm = TRUE), 
      n_genes = n(),
      .groups = "drop" 
    )
  return(cor_by_specificity_summary)
}
