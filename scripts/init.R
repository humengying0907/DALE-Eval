
train_test_split <- function(obs,
                             training_ratio = 0.2,
                             colname_of_cellType = 'cell_type',
                             colname_of_sample = 'sample',
                             min_nCells_per_ct_reference = 100,
                             min_nSamples_per_cell_type = 10,
                             max_iterations = 20) {
  require(dplyr)  
  
  unique_sampleIDs <- unique(obs[[colname_of_sample]])
  iterations <- 0
  
  while (iterations < max_iterations) {
    iterations <- iterations + 1
    
    train_sampleIDs <- sample(unique_sampleIDs, 
                              size = floor(length(unique_sampleIDs) * training_ratio), 
                              replace = FALSE)
    
    training_obs <- obs[obs[[colname_of_sample]] %in% train_sampleIDs, ]
    test_obs <- obs[!obs[[colname_of_sample]] %in% train_sampleIDs, ]
    
    # Count cells for each cell type in the training set
    n <- training_obs %>% 
      group_by_at(colname_of_cellType) %>% 
      summarize(n = n(), .groups = "drop") %>% 
      pull(n)
    
    # Count number of samples for each cell type in the test set
    nSamples <- test_obs %>% 
      group_by_at(colname_of_cellType) %>% 
      summarize(n_samples = n_distinct(!!sym(colname_of_sample)), .groups = "drop") %>% 
      pull(n_samples)
    
    # Check conditions
    if (all(n >= min_nCells_per_ct_reference, na.rm = TRUE) && 
        all(nSamples >= min_nSamples_per_cell_type, na.rm = TRUE)) {
      test_sampleIDs <- setdiff(unique_sampleIDs, train_sampleIDs)
      
      return(list(train_sampleIDs = as.character(train_sampleIDs), 
                  test_sampleIDs = as.character(test_sampleIDs)))
    }
  }
  
  # If the maximum number of iterations is reached and no valid split found
  stop("Unable to satisfy conditions for training and test splits. ",
       "Please check if the cell types and sample sizes meet the requirements.")
}


create_pseudobulk_obj = function(scExpr,scMeta,colnames_of_sample = 'sample',
                                 colnames_of_cellType = 'cell_type',
                                 unit = 'UMI'){
  require(deconvBenchmarking)
  
  stopifnot(unit %in% c('UMI','cpm'))
  stopifnot(all.equal(rownames(scMeta),colnames(scExpr)))
  
  if(max(scExpr)<50){
    stop('single cell input seems to be log transformed, please provide UMI or CPM')
  }
  
  cell_usage = scMeta[,colnames_of_sample,drop = F] %>% rownames_to_column('cell')
  colnames(cell_usage)[2] = 'bulk_id'
  sample_identifiers = scMeta[,colnames_of_sample]
  sample_identifiers = as.character(sample_identifiers)
  
  # use cell count as approximation to fractions
  Meta = data.frame(row.names = rownames(scMeta),
                    cell_type = scMeta[,colnames_of_cellType],
                    sampleID = scMeta[,colnames_of_sample])
  
  cell_usage$cell_type = Meta$cell_type[match(cell_usage$cell,rownames(Meta))]
  
  cell_type_table = Meta %>%
    group_by(sampleID, cell_type) %>%
    summarise(n=n()) %>%
    mutate(frac=n/sum(n)) %>% as.data.frame()
  frac = spread(cell_type_table[,-3], cell_type, frac) %>% column_to_rownames('sampleID') %>% as.matrix()
  frac[is.na(frac)] = 0
  
  if(unit == 'cpm'){
    # rowMeans of cells from the same sample
    pseudobulk = get_mu_by_group(scExpr,sample_identifiers)
    pseudobulk = pseudobulk[,rownames(frac)]
    
    pseudobulk_obj = list(simulated_bulk = pseudobulk,
                          simulated_frac = frac,
                          cell_usage = cell_usage)
    
  }else if(unit == 'UMI'){
    col_sums <- Matrix::colSums(scExpr)
    scaling_factors <- 10^6 / col_sums
    cpm <- scExpr %*% Matrix::Diagonal(x = scaling_factors)
    colnames(cpm) = colnames(scExpr)
    pseudobulk = get_mu_by_group(cpm,sample_identifiers)
    pseudobulk = pseudobulk[,rownames(frac)]
    
    pseudobulk_obj= list(simulated_bulk = pseudobulk,
                         simulated_frac = frac,
                         cell_usage = cell_usage)
    
  }
  
  # raise warning if medians of simulated_bulk is close to 0: indicating not enough cells to create the pseudobulk objects
  col_medians = apply(pseudobulk_obj$simulated_bulk,2,median)
  if(any(col_medians<0.1)){
    poor_quality_count = sum(col_medians < 0.1)
    warning(paste(poor_quality_count,'samples have expression median less than 0.1'))
  }
  
  return(pseudobulk_obj)
}

# this is a versatile method for calculating Z_truth, considering that not all cells in scExpr will be used for aggregation
Z_ground_truth = function(scExpr, 
                          simulated_bulk,
                          cell_usage,
                          unit = 'cpm',
                          n.core = 16){
  
  # for cell types with frac = 0, Z truth is also 0!
  
  scExpr = scExpr[rownames(simulated_bulk),]
  bulk_ids = colnames(simulated_bulk)
  bulk_ct_components = unique(cell_usage$cell_type)
  
  Z_truth = array(NA,dim = c(nrow(simulated_bulk),
                             ncol(simulated_bulk),
                             length(bulk_ct_components)))
  
  dimnames(Z_truth)[[1]] = rownames(simulated_bulk)
  dimnames(Z_truth)[[2]] = colnames(simulated_bulk)
  dimnames(Z_truth)[[3]] = bulk_ct_components
  
  
  if(unit == 'UMI'){
    col_sums <- Matrix::colSums(scExpr)
    scaling_factors <- 10^6 / col_sums
    cpm <- scExpr %*% Matrix::Diagonal(x = scaling_factors)
    colnames(cpm) = colnames(scExpr)
    
    scExpr = cpm
  }
  
  Z_ct_single = function(bulk_id,scExpr,cell_usage_sub){
    cells = cell_usage_sub$cell[cell_usage_sub$bulk_id == bulk_id]
    if(length(cells)==0){
      avg_expr = rep(0,nrow(scExpr))
      names(avg_expr) = rownames(scExpr)
      return(avg_expr)
    }else{
      avg_expr = Matrix::rowMeans(scExpr[,cells,drop=F])
      return(avg_expr)
    }
  }
  
  Z_ct_multi = function(bulk_ids,scExpr,cell_usage_sub,n.core){
    
    mat = do.call(cbind,pblapply(bulk_ids,Z_ct_single,scExpr,cell_usage_sub,cl = n.core))
    colnames(mat) = bulk_ids
    return(mat)
  }
  
  for(i in 1:length(bulk_ct_components)){
    ct = bulk_ct_components[i]
    print(ct)
    cell_usage_sub = cell_usage[cell_usage$cell_type == ct,]
    
    Z_truth[,,i] = Z_ct_multi(bulk_ids,scExpr,cell_usage_sub,n.core)
  }
  
  return(Z_truth)
  
}

# an efficient way to calculate Z_truth, which expects all cells in scExpr to be used in aggregation
Z_ground_truth_scalable = function(scExpr, 
                           cell_table_labels,
                           sample_labels,
                           simulated_bulk,
                           unit = 'cpm'){
  
  scExpr = scExpr[rownames(simulated_bulk),]
  bulk_ids = colnames(simulated_bulk)
  bulk_ct_components = unique(cell_table_labels)
  
  Z_truth = array(NA,dim = c(nrow(simulated_bulk),
                             ncol(simulated_bulk),
                             length(bulk_ct_components)))
  
  dimnames(Z_truth)[[1]] = rownames(simulated_bulk)
  dimnames(Z_truth)[[2]] = colnames(simulated_bulk)
  dimnames(Z_truth)[[3]] = bulk_ct_components
  
  if(unit == 'UMI'){
    col_sums <- Matrix::colSums(scExpr)
    scaling_factors <- 10^6 / col_sums
    cpm <- scExpr %*% Matrix::Diagonal(x = scaling_factors)
    colnames(cpm) = colnames(scExpr)
    
    scExpr = cpm
  }
  
  for(ct in bulk_ct_components){
    print(ct)
    ct_indices = which(cell_table_labels == ct)
    scExpr_sub = scExpr[, ct_indices]
    
    pseudobulk_ct = get_mu_by_group(scExpr_sub,sample_labels[ct_indices])
    
    for(bulk in colnames(pseudobulk_ct)){
      Z_truth[,bulk,ct] = pseudobulk_ct[,bulk]
    }
  }
  
  # for cell types with frac = 0, Z truth is also 0!
  Z_truth[is.na(Z_truth)] = 0
  
  return(Z_truth)
}


ct_recode = function(cell_types, mapping_list) {
  matched <- sapply(cell_types, function(ct) {
    matched_name <- names(mapping_list)[sapply(mapping_list, function(x) ct %in% x)]
    if (length(matched_name) == 0) NA else matched_name
  })
  return(unname(unlist(matched)))
}
