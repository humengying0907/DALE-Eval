
library(deconvBenchmarking)
source('../../scripts/init.R')
source('../../scripts/helpers.R')
source('../../scripts/evalu.R')

SC_DATA_PATH="/path/to/your/scRNA_seq/data.h5ad" 
ad = anndata::read_h5ad(SC_DATA_PATH)

########## train/test splitting ##########
obs = ad$obs
splitting = train_test_split(obs)
splitting_df = data.frame(group = c('train','test'),
                          sampleIDs = c(paste(splitting$train_sampleIDs,collapse = ','),
                                        paste(splitting$test_sampleIDs,collapse = ',')))
write.table(splitting_df,file = 'init/sample_split.csv',sep = ',',row.names = F)

########### create pseudobulk profiles ###########
sample_split = read.delim('init/sample_split.csv',sep = ',')
testing_sampleIDs = strsplit(sample_split$sampleIDs[2], ",")[[1]]

scExpr = Matrix::t(ad$X)
scMeta = ad$obs

keep_id = which(scMeta$sample %in% testing_sampleIDs)

scExpr_test = scExpr[,keep_id]
scMeta_test = scMeta[keep_id,]

# after splitting, some genes may have all zeros in training data, should filter them out (minCells = 10)
# this will ensure consistency between genes in pseudobulk and reference!
training_sampleIDs = strsplit(sample_split$sampleIDs[1], ",")[[1]]
scExpr_train = scExpr[,scMeta$sample %in% training_sampleIDs]
valid_genes = rownames(scExpr_train)[rowSums(scExpr_train >0)>=10]

benchmarking_obj = create_pseudobulk_obj(scExpr_test[valid_genes,],scMeta_test,unit = 'UMI')

write.table(benchmarking_obj$simulated_frac,file = 'init/truthFrac.csv',sep = ',',row.names = T)
write.table(benchmarking_obj$simulated_bulk,file = 'init/bulk_expr.csv',sep = ',',row.names = T)

Z_truth = Z_ground_truth_scalable(scExpr_test,
                                  scMeta_test$cell_type,
                                  scMeta_test$sample,
                                  benchmarking_obj$simulated_bulk,
                                  unit = 'UMI')

X_reconstructed = get_reconstructed_X_from_weighted_Z(Z_truth,benchmarking_obj$simulated_frac)
all.equal(benchmarking_obj$simulated_bulk,X_reconstructed) # TRUE
saveRDS(Z_truth,file = 'init/Z_truth.RDS')

# get Z_truth derived specificity (DE analysis between cell types in flattened Z-truth)
Z_truth_limma_top_genes = Z_limma_statistic(Z_truth)
write.table(Z_truth_limma_top_genes,file = 'init/Z_truth_limma_top_genes.csv',sep = ',')

